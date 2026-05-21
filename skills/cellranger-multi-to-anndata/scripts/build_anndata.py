"""build_anndata.py — parameterised CellRanger Multi → AnnData ingestion.

Ported from <ref-scrna>'s 02_Analysis/00_build_anndata.py +
01_Scripts/Python_scripts/anndata_utils.py. Decision Pauses (metadata source,
join key, species) live in the SKILL.md body; this script is the pure
mechanical part — given resolved parameters, build the .h5ad.

Usage from a notebook or a numbered script:

    from build_anndata import build
    adata = build(
        cellranger_root=Path("00_data/raw/CellRanger"),
        metadata_tsv=Path("00_data/raw/SamplesMetadata.txt"),
        metadata_join_key="Sample",
        species="mouse",
        out_path=Path("03_results/checkpoints/00_raw.h5ad"),
    )

CLI (for batch / cron use):

    python build_anndata.py \\
        --cellranger-root 00_data/raw/CellRanger \\
        --metadata-tsv 00_data/raw/SamplesMetadata.txt \\
        --metadata-join-key Sample \\
        --species mouse \\
        --out 03_results/checkpoints/00_raw.h5ad
"""

from __future__ import annotations

import argparse
import glob
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


# --------------------------------------------------------------------------- #
# species defaults — Decision Pause 3 has resolved this before build() runs
# --------------------------------------------------------------------------- #

SPECIES_DEFAULTS = {
    "mouse": {"organism": "mmusculus", "mt": "mt-",  "ribo": ("Rpl", "Rps")},
    "human": {"organism": "hsapiens",  "mt": "MT-",  "ribo": ("RPL", "RPS")},
}


# --------------------------------------------------------------------------- #
# step 2 — read robust per sample
# --------------------------------------------------------------------------- #

def read_10x_h5_robust(
    fp: str | Path,
    sample_id: str,
    pool_id: str,
    manif_row: pd.DataFrame | None = None,
    verbose: bool = True,
) -> ad.AnnData:
    """Read one 10x .h5, tolerate column-name drift, prefix obs_names."""
    a = sc.read_10x_h5(str(fp), gex_only=True)

    # Detect ID and symbol column names across CellRanger versions
    vc = {c.lower(): c for c in a.var.columns}
    id_col = next(
        (vc[c] for c in ("gene_ids", "gene_id", "id", "feature_id") if c in vc), None,
    )
    name_col = next(
        (vc[c] for c in ("gene_symbols", "gene_symbol", "gene_names", "gene_name", "name", "feature_name") if c in vc),
        None,
    )

    a.var["gene_id"] = a.var[id_col].astype(str) if id_col else pd.Index(a.var_names).astype(str)
    a.var["gene_name"] = a.var[name_col].astype(str) if name_col else a.var["gene_id"]
    a.var_names = a.var["gene_id"]
    a.var_names_make_unique()

    # cell-level provenance
    a.obs["barcode"] = a.obs_names.astype(str)
    a.obs["sample_id"] = sample_id
    a.obs["pool_id"] = pool_id
    a.obs_names = [f"{sample_id}-{bc}" for bc in a.obs["barcode"]]

    # optional manifest row
    if manif_row is not None and not manif_row.empty:
        for col in set(manif_row.columns) - {"sample_id", "h5"}:
            a.obs[col] = manif_row.iloc[0][col]

    if verbose:
        symbols_differ = int((a.var["gene_name"] != a.var["gene_id"]).sum())
        print(f"  {sample_id}: {a.n_obs} cells × {a.n_vars} genes "
              f"(symbols resolved for {symbols_differ}/{a.n_vars})")
    return a


# --------------------------------------------------------------------------- #
# step 3 — gene-union concat
# --------------------------------------------------------------------------- #

def robust_concat_with_gene_union(
    adatas: list[ad.AnnData],
    keys: list[str],
    label: str = "orig_ident",
) -> ad.AnnData:
    """Union-of-genes concat with zero-fill and preserved gene metadata.

    See references/gene-union-concat.md for the full pattern and the
    sparse-vs-dense memory caveat.
    """
    if not adatas:
        raise ValueError("No AnnData objects to concatenate")

    # 1. union of genes + first-occurrence metadata
    all_genes: set[str] = set()
    gene_metadata: dict[str, dict] = {}
    for a in adatas:
        for g in a.var_names:
            all_genes.add(g)
            gene_metadata.setdefault(g, a.var.loc[g].to_dict())
    all_genes_sorted = sorted(all_genes)

    # 2. union var table
    union_var = pd.DataFrame(index=all_genes_sorted)
    union_var["gene_id"] = union_var.index.astype(str)
    sym_cols = ("gene_symbol", "gene_symbols", "gene_name", "feature_name")
    union_var["gene_name"] = [
        next(
            (str(gene_metadata[g][c]) for c in sym_cols
             if c in gene_metadata[g] and gene_metadata[g][c]
             and not pd.isna(gene_metadata[g][c])
             and str(gene_metadata[g][c]) != str(g)),
            g,
        )
        for g in all_genes_sorted
    ]

    # 3. re-index each input, zero-fill missing genes
    standardized: list[ad.AnnData] = []
    for a in adatas:
        missing = [g for g in all_genes_sorted if g not in a.var_names]
        if missing:
            n_cells, n_missing = a.n_obs, len(missing)
            if sparse.isspmatrix(a.X):
                zero = sparse.csr_matrix((n_cells, n_missing))
                new_X = sparse.hstack([a.X, zero])
            else:
                zero = np.zeros((n_cells, n_missing))
                new_X = np.hstack([a.X, zero])
            new_var = pd.concat([a.var, union_var.loc[missing]])
            new = ad.AnnData(
                X=new_X, obs=a.obs.copy(), var=new_var,
                obsm=dict(a.obsm), uns=dict(a.uns),
            )
            new.var_names = list(a.var_names) + missing
        else:
            new = a.copy()
        new = new[:, all_genes_sorted].copy()
        new.var = union_var.copy()
        standardized.append(new)

    # 4. concat
    return ad.concat(
        standardized,
        axis=0,
        join="outer",
        fill_value=0,
        merge="first",
        uns_merge="unique",
        label=label,
        keys=keys,
        index_unique=None,
    )


# --------------------------------------------------------------------------- #
# step 4 — biomart symbol annotation
# --------------------------------------------------------------------------- #

def annotate_genes_biomart(
    adata: ad.AnnData,
    organism: str = "mmusculus",
    use_cache: bool = True,
    verbose: bool = True,
) -> ad.AnnData:
    """Map Ensembl IDs to symbols + biotype + MT flag via scanpy.queries.biomart."""
    try:
        annotations = sc.queries.biomart_annotations(
            organism,
            ["ensembl_gene_id", "external_gene_name", "gene_biotype"],
            use_cache=use_cache,
        ).set_index("ensembl_gene_id")

        name_map = annotations["external_gene_name"].to_dict()
        adata.var["gene_name"] = [name_map.get(g, g) for g in adata.var_names]

        if "gene_biotype" in annotations.columns:
            biotype_map = annotations["gene_biotype"].to_dict()
            adata.var["gene_biotype"] = [biotype_map.get(g, "unknown") for g in adata.var_names]

        try:
            mt = sc.queries.mitochondrial_genes(organism, attrname="ensembl_gene_id")
            adata.var["mt_biomart"] = adata.var_names.isin(set(mt["ensembl_gene_id"]))
        except Exception as e:
            if verbose:
                print(f"  MT-genes biomart call failed: {e}")
            adata.var["mt_biomart"] = False

        if verbose:
            mapped = int((adata.var["gene_name"] != adata.var_names).sum())
            print(f"biomart: {mapped}/{adata.n_vars} symbols resolved "
                  f"(mt_biomart count: {int(adata.var['mt_biomart'].sum())})")
    except Exception as e:
        if verbose:
            print(f"biomart annotation failed: {e}; keeping Ensembl IDs as gene_name")
        adata.var["gene_name"] = adata.var_names.astype(str)
        adata.var["mt_biomart"] = False
    return adata


# --------------------------------------------------------------------------- #
# step 5 — metadata join (Decision Pauses 1+2 resolved upstream)
# --------------------------------------------------------------------------- #

def join_metadata(
    adata: ad.AnnData,
    metadata_path: Path,
    join_key: str,
    sep: str | None = None,
    verbose: bool = True,
) -> ad.AnnData:
    """Left-join metadata into obs after validation."""
    if sep is None:
        sep = "\t" if metadata_path.suffix in (".tsv", ".txt") else ","
    metadata = pd.read_csv(metadata_path, sep=sep)

    if join_key not in metadata.columns:
        raise KeyError(f"Join key '{join_key}' not in metadata columns: {list(metadata.columns)}")
    key_values = metadata[join_key].astype(str)
    if not key_values.is_unique:
        raise ValueError(f"Join key '{join_key}' has duplicates in {metadata_path}")

    adata_keys = set(adata.obs["sample_id"].astype(str).unique())
    meta_keys = set(key_values)
    unmatched_obs = adata_keys - meta_keys
    unmatched_meta = meta_keys - adata_keys
    if unmatched_obs:
        print(f"WARNING: {len(unmatched_obs)} samples in AnnData not in metadata: "
              f"{sorted(unmatched_obs)[:5]}...")
    if unmatched_meta:
        print(f"WARNING: {len(unmatched_meta)} metadata rows not in AnnData: "
              f"{sorted(unmatched_meta)[:5]}...")

    obs_idx = adata.obs.index
    merged = adata.obs.merge(
        metadata, left_on="sample_id", right_on=join_key,
        how="left", validate="many_to_one",
    )
    merged.index = obs_idx
    adata.obs = merged
    if verbose:
        new_cols = [c for c in metadata.columns if c != join_key]
        print(f"metadata join: added obs columns {new_cols}")
    return adata


# --------------------------------------------------------------------------- #
# step 1 — discover, then orchestrate
# --------------------------------------------------------------------------- #

@dataclass
class BuildResult:
    adata: ad.AnnData
    n_pools: int
    n_samples: int
    n_h5_files: int


def discover_h5s(
    cellranger_root: Path,
    pool_segment: str | None = None,
    filtered_or_raw: str = "filtered",
) -> list[tuple[str, str, Path]]:
    """Glob per_sample_outs and return [(sample_id, pool_id, path)] tuples."""
    matrix_name = (
        "sample_filtered_feature_bc_matrix.h5" if filtered_or_raw == "filtered"
        else "raw_feature_bc_matrix.h5"
    )
    pattern = str(cellranger_root / "**" / "per_sample_outs" / "*" / "*" / matrix_name)
    h5s = sorted(glob.glob(pattern, recursive=True))

    pool_seg = pool_segment or cellranger_root.name
    out: list[tuple[str, str, Path]] = []
    for fp in h5s:
        sid_match = re.search(r"per_sample_outs/([^/]+)/", fp)
        pool_match = re.search(rf"{re.escape(pool_seg)}/([^/]+)/", fp)
        if not sid_match or not pool_match:
            print(f"WARNING: skipping path that did not match regexes: {fp}")
            continue
        out.append((sid_match.group(1), pool_match.group(1), Path(fp)))
    return out


def build(
    cellranger_root: Path,
    out_path: Path,
    metadata_tsv: Path | None = None,
    metadata_join_key: str | None = None,
    species: str = "mouse",
    species_overrides: dict | None = None,
    pool_segment: str | None = None,
    filtered_or_raw: str = "filtered",
    biomart_use_cache: bool = True,
    seed: int = 42,
) -> ad.AnnData:
    """Run the full ingestion pipeline. Decision Pauses resolved upstream."""
    np.random.seed(seed)
    sc.settings.verbosity = 2

    cellranger_root = Path(cellranger_root)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # 1. discover
    triples = discover_h5s(cellranger_root, pool_segment=pool_segment,
                           filtered_or_raw=filtered_or_raw)
    if not triples:
        raise FileNotFoundError(
            f"No matrices under {cellranger_root}. Check pool_segment "
            f"(currently '{pool_segment or cellranger_root.name}') and the glob in "
            "references/pool-detection.md."
        )
    print(f"Discovered {len(triples)} matrices "
          f"across {len(set(t[1] for t in triples))} pools, "
          f"{len(set(t[0] for t in triples))} samples")

    # optional manifest pre-load (passed per-row into read_10x_h5_robust)
    manif: pd.DataFrame | None = None
    if metadata_tsv and metadata_join_key:
        sep = "\t" if Path(metadata_tsv).suffix in (".tsv", ".txt") else ","
        manif = pd.read_csv(metadata_tsv, sep=sep)
        # Only carry the join column as `sample_id` for per-sample seeding;
        # full join happens after concat to keep concat metadata clean.
        manif = manif.rename(columns={metadata_join_key: "sample_id"})

    # 2. read robust per sample
    adatas, keys = [], []
    for sid, pool, fp in triples:
        row = manif.loc[manif["sample_id"] == sid] if manif is not None else None
        adatas.append(read_10x_h5_robust(fp, sid, pool, row))
        keys.append(sid)

    # 3. gene-union concat
    adata = robust_concat_with_gene_union(adatas, keys=keys, label="orig_ident")
    print(f"Concatenated: {adata.n_obs} cells × {adata.n_vars} genes")

    # 4. biomart symbols
    species_cfg = (species_overrides or {}) if species == "other" else SPECIES_DEFAULTS.get(species)
    if not species_cfg:
        raise ValueError(f"Unknown species '{species}'. Use 'mouse', 'human', or 'other' with species_overrides.")
    adata = annotate_genes_biomart(adata, organism=species_cfg["organism"], use_cache=biomart_use_cache)

    # 5. metadata join (full)
    if metadata_tsv and metadata_join_key:
        adata = join_metadata(adata, Path(metadata_tsv), metadata_join_key)

    # 6. write checkpoint
    adata.write_h5ad(out_path, compression="gzip")
    print(f"Wrote {out_path} ({adata.n_obs} × {adata.n_vars})")
    return adata


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def _parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build a single AnnData from CellRanger Multi outputs.")
    p.add_argument("--cellranger-root", type=Path, required=True)
    p.add_argument("--out", dest="out_path", type=Path, required=True)
    p.add_argument("--metadata-tsv", type=Path, default=None)
    p.add_argument("--metadata-join-key", type=str, default=None)
    p.add_argument("--species", choices=("mouse", "human", "other"), default="mouse")
    p.add_argument("--pool-segment", type=str, default=None,
                   help="Path segment used to extract pool_id (default: basename of --cellranger-root)")
    p.add_argument("--filtered-or-raw", choices=("filtered", "raw"), default="filtered")
    p.add_argument("--biomart-no-cache", action="store_true",
                   help="Force-refresh biomart cache (use after suspected drift)")
    return p.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> int:
    args = _parse_args(argv)
    if (args.metadata_tsv is None) ^ (args.metadata_join_key is None):
        raise SystemExit("--metadata-tsv and --metadata-join-key must be provided together")
    build(
        cellranger_root=args.cellranger_root,
        out_path=args.out_path,
        metadata_tsv=args.metadata_tsv,
        metadata_join_key=args.metadata_join_key,
        species=args.species,
        pool_segment=args.pool_segment,
        filtered_or_raw=args.filtered_or_raw,
        biomart_use_cache=not args.biomart_no_cache,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
