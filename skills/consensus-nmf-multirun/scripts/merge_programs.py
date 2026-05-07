"""merge_programs.py — Stage 5, port of 13403-YD's 07_merge_programs.py.

Hierarchical clustering at r > threshold (default 0.7), rank aggregation of
top-50 genes per source program -> top-100 consensus, classification by gene
panel composition (Biological / Technical / CellCycle / Ribosomal /
Mitochondrial / ImmediateEarly).

Bug preserved: source prefixes are sorted by length descending before any
startswith() match (so cNMF_fullQC is matched before cNMF_full).
"""

from __future__ import annotations

import argparse
import pickle
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform


# --------------------------------------------------------------------------- #
# Gene panels — mouse (default). Set --species human to swap.
# --------------------------------------------------------------------------- #

PANELS_MOUSE = {
    "ribo_pattern": r"^Rp[sl]\d*[a-z]?$|^Rplp\d$|^Rpsa$",
    "mt_pattern":   r"^mt-",
    "cell_cycle":   {
        "Top2a", "Mki67", "Cdk1", "Ccna2", "Ccnb1", "Ccnb2", "Cdc20",
        "Birc5", "Ube2c", "Nusap1", "Cenpf", "Tpx2", "Plk1", "Aurka",
        "Aurkb", "Bub1", "Bub1b", "Esco2", "Ncapd2", "Ncapg", "Smc2", "Smc4",
    },
    "ieg": {
        "Fos", "Fosb", "Jun", "Junb", "Jund", "Egr1", "Egr2", "Egr3",
        "Nr4a1", "Nr4a2", "Nr4a3", "Ier2", "Ier3", "Atf3", "Btg2",
    },
    "housekeeping": {
        "Gapdh", "Actb", "Actg1", "Tuba1b", "Tubb5", "B2m", "Hsp90ab1",
        "Eef1a1", "Eef2", "Aldoa", "Pgk1", "Ldha", "Eno1", "Tpi1",
    },
}

PANELS_HUMAN = {
    "ribo_pattern": r"^RP[SL]\d*[A-Z]?$|^RPLP\d$|^RPSA$",
    "mt_pattern":   r"^MT-",
    "cell_cycle":   {g.upper() for g in PANELS_MOUSE["cell_cycle"]},
    "ieg":          {g.upper() for g in PANELS_MOUSE["ieg"]},
    "housekeeping": {g.upper() for g in PANELS_MOUSE["housekeeping"]},
}


def get_panels(species: str) -> dict:
    if species.lower() in ("mouse", "mmusculus"):
        return PANELS_MOUSE
    if species.lower() in ("human", "hsapiens"):
        return PANELS_HUMAN
    raise ValueError(f"Unknown species {species!r}; pass 'mouse' or 'human'.")


# --------------------------------------------------------------------------- #
# Source prefix priority (length-descending — bug-fix)
# --------------------------------------------------------------------------- #

def make_source_priority(variant_names: list[str]) -> list[str]:
    """Build prefix list `cNMF_<variant>` sorted by length descending."""
    prefixes = [f"cNMF_{v}" for v in variant_names] + ["cNMF_"]
    return sorted(set(prefixes), key=len, reverse=True)


def get_source_from_program(program_name: str, sorted_prefixes: list[str]) -> str:
    for p in sorted_prefixes:
        if program_name.startswith(p):
            return p
    return "unknown"


# --------------------------------------------------------------------------- #
# Classification
# --------------------------------------------------------------------------- #

def classify_program(genes: list[str], panels: dict, n_top: int = 50) -> tuple[str, dict]:
    """Return (class, fractions) for a program's top genes."""
    if not genes:
        return "Unknown", {}
    top = genes[:n_top]
    n = len(top)

    ribo  = sum(1 for g in top if re.match(panels["ribo_pattern"], g, re.IGNORECASE))
    mito  = sum(1 for g in top if re.match(panels["mt_pattern"], g, re.IGNORECASE))
    cc    = sum(1 for g in top if g in panels["cell_cycle"])
    ieg   = sum(1 for g in top if g in panels["ieg"])
    hk    = sum(1 for g in top if g in panels["housekeeping"])

    fracs = {
        "pct_ribosomal":   round(ribo / n * 100, 1),
        "pct_mito":        round(mito / n * 100, 1),
        "pct_cell_cycle":  round(cc   / n * 100, 1),
        "pct_ieg":         round(ieg  / n * 100, 1),
        "pct_housekeeping":round(hk   / n * 100, 1),
        "pct_technical":   round((ribo + mito + hk) / n * 100, 1),
    }

    if cc / n > 0.20:                                         return "CellCycle", fracs
    if ribo / n > 0.20:                                       return "Ribosomal", fracs
    if mito / n > 0.20:                                       return "Mitochondrial", fracs
    if (ribo + mito + hk) / n > 0.30:                         return "Technical", fracs
    if ieg / n > 0.15:                                        return "ImmediateEarly", fracs
    return "Biological", fracs


# --------------------------------------------------------------------------- #
# Rank aggregation
# --------------------------------------------------------------------------- #

def compute_consensus_genes(
    programs: list[str],
    program_data: dict[str, dict],
    top_genes_to_keep: int = 100,
    top_genes_per_source: int = 50,
) -> tuple[list[dict], str]:
    """Return ([{gene, consensus_score, n_sources, avg_rank}, ...], merge_method)."""
    if len(programs) == 1:
        prog = programs[0]
        if prog not in program_data:
            return [], "singleton"
        return [
            {"gene": g, "consensus_score": 1.0, "n_sources": 1, "avg_rank": rank + 1}
            for rank, g in enumerate(program_data[prog]["genes"][:top_genes_to_keep])
        ], "singleton"

    n_progs = len(programs)
    max_rank = top_genes_per_source
    stats: dict[str, dict] = defaultdict(lambda: {"ranks": [], "sources": []})

    for prog in programs:
        if prog not in program_data:
            continue
        for rank, gene in enumerate(program_data[prog]["genes"][:max_rank]):
            stats[gene]["ranks"].append(rank + 1)
            stats[gene]["sources"].append(prog)

    out: list[dict] = []
    for gene, st in stats.items():
        n_src = len(st["sources"])
        avg_rank = float(np.mean(st["ranks"]))
        rank_score = (max_rank - avg_rank + 1) / max_rank
        coverage_score = n_src / n_progs
        score = rank_score * (1 + coverage_score)
        out.append({"gene": gene, "consensus_score": score,
                    "n_sources": n_src, "avg_rank": avg_rank})
    out.sort(key=lambda x: x["consensus_score"], reverse=True)
    return out[:top_genes_to_keep], "consensus"


def select_canonical(programs: list[str], sorted_prefixes: list[str]) -> str:
    for prefix in sorted_prefixes:
        for prog in programs:
            if prog.startswith(prefix):
                return prog
    return programs[0]


def determine_confidence(n_unique_runs: int) -> str:
    if n_unique_runs >= 4: return "High"
    if n_unique_runs >= 2: return "Medium"
    return "Low"


# --------------------------------------------------------------------------- #
# Main entry point
# --------------------------------------------------------------------------- #

@dataclass
class MergeResult:
    merged_programs_path: Path
    program_annotations_path: Path
    gene_consensus_path: Path
    n_clusters: int


def load_program_data(
    sources_pkls: dict[str, str | Path],
    sorted_prefixes: list[str],
    top_genes_per_source: int = 50,
) -> dict[str, dict]:
    """Build {program_name: {genes, gene_scores, source, run}} from per-variant pkls."""
    program_data: dict[str, dict] = {}
    for variant_name, pkl_path in sources_pkls.items():
        pkl_path = Path(pkl_path)
        if not pkl_path.exists():
            print(f"  skipping {variant_name}: {pkl_path} not found")
            continue
        with open(pkl_path, "rb") as f:
            res = pickle.load(f)
        top_genes_df: pd.DataFrame = res["top_genes"]
        spectra = res.get("spectra_scores") or res.get("spectra_tpm")
        for i, col in enumerate(top_genes_df.columns):
            program_name = f"cNMF_{variant_name}_P{i+1}"
            genes = top_genes_df[col].head(top_genes_per_source).tolist()
            gene_scores: dict[str, dict] = {}
            if spectra is not None:
                try:
                    s = spectra.loc[:, col].sort_values(ascending=False)
                    for rank, (g, sc) in enumerate(s.head(top_genes_per_source).items()):
                        gene_scores[g] = {"score": float(sc), "rank": rank + 1}
                except Exception:
                    pass
            program_data[program_name] = {
                "genes": genes,
                "gene_scores": gene_scores,
                "source": get_source_from_program(program_name, sorted_prefixes),
                "run": variant_name,
            }
    return program_data


def merge(
    correlation_csv: str | Path,
    sources_pkls: dict[str, str | Path],
    out_dir: str | Path,
    correlation_threshold: float = 0.7,
    top_genes_to_keep: int = 100,
    top_genes_per_source: int = 50,
    species: str = "mouse",
) -> MergeResult:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    correlation_csv = Path(correlation_csv)
    panels = get_panels(species)

    corr_df = pd.read_csv(correlation_csv, index_col=0)
    program_names: list[str] = list(corr_df.index)
    sorted_prefixes = make_source_priority(list(sources_pkls.keys()))
    program_data = load_program_data(sources_pkls, sorted_prefixes, top_genes_per_source)

    # Hierarchical clustering on 1 - |r|
    corr = corr_df.values.copy()
    nan_mask = np.isnan(corr)
    if nan_mask.any():
        print(f"  WARNING: {int(nan_mask.sum())} NaN entries in correlation matrix; filling with 0")
        corr[nan_mask] = 0.0
    dist = 1 - np.abs(corr)
    np.fill_diagonal(dist, 0)
    dist = (dist + dist.T) / 2
    condensed = squareform(dist, checks=False)
    Z = linkage(condensed, method="ward")
    cutoff = 1 - correlation_threshold
    clusters = fcluster(Z, t=cutoff, criterion="distance")
    n_clusters = int(len(set(clusters)))
    print(f"  {n_clusters} clusters at r > {correlation_threshold}")

    # Group programs by cluster
    cluster_groups: dict[int, list[str]] = defaultdict(list)
    for prog, c in zip(program_names, clusters):
        cluster_groups[c].append(prog)

    merged_rows: list[dict] = []
    annotation_rows: list[dict] = []
    gene_rows: list[dict] = []

    for cluster_id in sorted(cluster_groups):
        progs = cluster_groups[cluster_id]
        runs = {program_data[p]["run"] for p in progs if p in program_data}
        canonical = select_canonical(progs, sorted_prefixes)
        confidence = determine_confidence(len(runs))
        gene_info, method = compute_consensus_genes(
            progs, program_data,
            top_genes_to_keep=top_genes_to_keep,
            top_genes_per_source=top_genes_per_source,
        )
        top_genes_list = [g["gene"] for g in gene_info]
        cls, fracs = classify_program(top_genes_list, panels, n_top=top_genes_per_source)

        canonical_name = f"merged_{cluster_id:02d}_{canonical}"
        merged_rows.append({
            "cluster_id": cluster_id,
            "canonical_program": canonical,
            "merged_name": canonical_name,
            "n_programs": len(progs),
            "n_unique_runs": len(runs),
            "confidence": confidence,
            "class": cls,
            "merge_method": method,
            "source_programs": ",".join(progs),
            "top_genes": ",".join(top_genes_list),
        })
        annotation_rows.append({"merged_name": canonical_name, "class": cls, **fracs})
        for g in gene_info:
            gene_rows.append({"merged_name": canonical_name, **g})

    merged_df = pd.DataFrame(merged_rows)
    annot_df = pd.DataFrame(annotation_rows)
    gene_df = pd.DataFrame(gene_rows)

    merged_path = out_dir / "merged_programs_v2.csv"
    annot_path  = out_dir / "program_annotations.csv"
    gene_path   = out_dir / "gene_consensus_scores.csv"
    merged_df.to_csv(merged_path, index=False)
    annot_df.to_csv(annot_path, index=False)
    gene_df.to_csv(gene_path, index=False)

    print(f"  wrote {merged_path}, {annot_path}, {gene_path}")
    print(f"  confidence distribution: {merged_df['confidence'].value_counts().to_dict()}")
    print(f"  class distribution:      {merged_df['class'].value_counts().to_dict()}")

    return MergeResult(
        merged_programs_path=merged_path,
        program_annotations_path=annot_path,
        gene_consensus_path=gene_path,
        n_clusters=n_clusters,
    )


def _parse() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Merge cross-variant cNMF programs.")
    p.add_argument("--correlation-csv", type=Path, required=True)
    p.add_argument("--source", action="append", required=True, metavar="NAME=PATH")
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--threshold", type=float, default=0.7)
    p.add_argument("--top-genes-to-keep", type=int, default=100)
    p.add_argument("--top-genes-per-source", type=int, default=50)
    p.add_argument("--species", choices=("mouse", "human"), default="mouse")
    return p.parse_args()


def main() -> int:
    a = _parse()
    sources: dict[str, str] = {}
    for entry in a.source:
        name, path = entry.split("=", 1)
        sources[name] = path
    merge(
        correlation_csv=a.correlation_csv,
        sources_pkls=sources,
        out_dir=a.out_dir,
        correlation_threshold=a.threshold,
        top_genes_to_keep=a.top_genes_to_keep,
        top_genes_per_source=a.top_genes_per_source,
        species=a.species,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
