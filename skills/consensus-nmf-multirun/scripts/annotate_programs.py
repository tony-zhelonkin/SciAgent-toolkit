"""annotate_programs.py — Stage 6, port of <ref-scrna>'s 07_annotate_programs.py.

g:Profiler GO/KEGG/Reactome enrichment on top-50 genes per merged program.
Per-program CSV + combined `program_annotations.csv` (overwrites the merge
step's output of the same name; the merge file is kept as `_classes.csv`).

Rate-limit handling: sleep 1s on Exception; retry once.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import pandas as pd


def annotate(
    merged_csv: str | Path,
    organism: str = "mmusculus",
    out_csv: str | Path = "program_annotations_gprofiler.csv",
    out_dir: str | Path | None = None,
    top_n_genes: int = 50,
    sources: list[str] | None = None,
    sleep_on_error: float = 1.0,
) -> Path:
    """Run g:Profiler per merged program, write per-program + combined CSVs."""
    try:
        from gprofiler import GProfiler
    except ImportError as e:
        raise ImportError(
            "gprofiler-official is not installed. `pip install gprofiler-official`."
        ) from e

    merged_csv = Path(merged_csv)
    out_csv = Path(out_csv)
    out_dir = Path(out_dir) if out_dir else out_csv.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    sources = sources or ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"]

    merged = pd.read_csv(merged_csv)
    if "merged_name" not in merged.columns or "top_genes" not in merged.columns:
        raise ValueError(f"merged_csv {merged_csv} must have 'merged_name' and 'top_genes' columns")

    gp = GProfiler(return_dataframe=True)
    all_results: list[pd.DataFrame] = []

    for _, row in merged.iterrows():
        program_name = row["merged_name"]
        gene_list = [g.strip() for g in str(row["top_genes"]).split(",") if g.strip()][:top_n_genes]
        if not gene_list:
            print(f"  {program_name}: empty gene list; skipping")
            continue

        print(f"  {program_name}: querying g:Profiler ({len(gene_list)} genes)")

        result = None
        for attempt in range(2):
            try:
                result = gp.profile(
                    organism=organism,
                    query=gene_list,
                    sources=sources,
                    significance_threshold_method="fdr",
                )
                break
            except Exception as e:
                print(f"    error (attempt {attempt + 1}): {e}")
                if attempt == 0:
                    time.sleep(sleep_on_error)

        if result is None or result.empty:
            print("    no significant terms")
            continue

        result = result.sort_values("p_value")
        per_path = out_dir / f"program_{program_name}_GO.csv"
        result.to_csv(per_path, index=False)
        result = result.assign(program=program_name)
        all_results.append(result)
        top = result.iloc[0]
        print(f"    top: {top['name']} (p={top['p_value']:.2e})")

    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        combined.to_csv(out_csv, index=False)
        print(f"  wrote {out_csv} ({len(combined)} rows)")
    else:
        print("  no annotations produced for any program")
        # Still write an empty file with the expected header so downstream stages don't error
        pd.DataFrame(columns=[
            "program", "name", "source", "p_value", "intersection_size",
            "term_size", "query_size", "precision", "recall", "native",
        ]).to_csv(out_csv, index=False)
    return out_csv


def _parse() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="g:Profiler annotation of merged cNMF programs.")
    p.add_argument("--merged-csv", type=Path, required=True)
    p.add_argument("--organism", type=str, default="mmusculus",
                   help="Source from analysis_config.yaml::decisions::cellranger-multi-to-anndata::species")
    p.add_argument("--out-csv", type=Path, required=True)
    p.add_argument("--top-n-genes", type=int, default=50)
    return p.parse_args()


def main() -> int:
    a = _parse()
    annotate(
        merged_csv=a.merged_csv,
        organism=a.organism,
        out_csv=a.out_csv,
        top_n_genes=a.top_n_genes,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
