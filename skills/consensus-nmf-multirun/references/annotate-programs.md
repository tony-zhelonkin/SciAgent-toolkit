# Annotate programs via g:Profiler

Stage 6 takes the merged programs' top-50 gene lists and queries g:Profiler for GO / KEGG / Reactome enrichment. Output is one CSV per program plus a combined `program_annotations.csv`.

## The query

```python
from gprofiler import GProfiler
gp = GProfiler(return_dataframe=True)

result = gp.profile(
    organism="mmusculus",                                            # 'hsapiens' for human
    query=top_50_genes,                                              # list of gene symbols
    sources=["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"],
    significance_threshold_method="fdr",                             # BH-FDR over all sources
)
```

`return_dataframe=True` makes `gp.profile()` return a `pd.DataFrame` directly; without it, you get a list of dicts. Always set it.

`organism` should be consumed from `analysis_config.yaml::decisions::cellranger-multi-to-anndata::species` (the species was resolved in Stage 0). Hard-coding here couples the skill to one species; pulling from config keeps it portable.

## What comes back

A DataFrame with columns including `name` (term name), `p_value` (FDR-corrected), `intersection_size`, `term_size`, `query_size`, `source` (the database). One row per significant term across all five `sources`.

## Per-program output

```python
result_sorted = result.sort_values("p_value")
result_sorted.to_csv(out_dir / f"program_{program_name}_GO.csv", index=False)
print(f"  Top: {result_sorted.iloc[0]['name']} (p={result_sorted.iloc[0]['p_value']:.2e})")
```

Expected pattern: a real biological program returns dozens of significant terms; the top is interpretable (e.g., `T cell activation`, `G2-M cell cycle transition`). A program that returns zero significant terms is suspicious — likely the gene list is dominated by ribosomal / MT genes that g:Profiler rejects as background.

## Rate limits

g:Profiler's public API has soft rate limits (~10 queries/sec; bursts allowed). For ≤ 20 merged programs, the limit is invisible. For larger queries:

```python
import time

for program_name, genes in merged_programs.items():
    try:
        result = gp.profile(...)
    except Exception as e:
        print(f"  Error on {program_name}: {e}")
        time.sleep(1)                # back off
        try:
            result = gp.profile(...)  # retry once
        except Exception as e2:
            print(f"  Retry also failed: {e2}; skipping")
            continue
```

The retry is what the <ref-scrna> reference does. For batch jobs running >100 programs, consider switching to the asynchronous endpoint or self-hosting `gprofiler-official`.

## Combined output

```python
combined = pd.concat([
    df.assign(program=name) for name, df in all_annotations.items()
])
combined.to_csv(out_dir / "program_annotations.csv", index=False)
```

`combined` is the single file downstream skills (`bulk-rnaseq-pathway-explorer`) read. The schema:

| Column | Type | Notes |
|--------|------|-------|
| program | str | merged program canonical name |
| name | str | term name |
| source | str | one of `GO:BP`, `GO:MF`, `GO:CC`, `KEGG`, `REAC` |
| p_value | float | FDR-corrected |
| intersection_size | int | overlap count |
| term_size | int | term gene count |
| query_size | int | always 50 here |
| precision | float | overlap / query_size |
| recall | float | overlap / term_size |
| native | str | native term ID (e.g., `GO:0042110`) |

## Falling back when g:Profiler is unreachable

Long deploys (CI, air-gapped) may not reach `biit.cs.ut.ee`. Two options:

- **Self-host gprofiler-official** (Docker image available; documented at https://biit.cs.ut.ee/gprofiler/page/docs).
- **Use decoupler-py with the same gene panels** — `decoupler.enrich()` with MSigDB / Reactome from `omnipath`. Different stats; same shape of output.

The shipped script logs a warning when g:Profiler is unreachable and falls through to write empty CSVs (so downstream stages do not error on missing files).

## Cross-reference

g:Profiler's category prefixes:

- `GO:BP` — biological process
- `GO:MF` — molecular function
- `GO:CC` — cellular component
- `KEGG` — KEGG pathways
- `REAC` — Reactome pathways

Anton's reference picks these five and excludes WikiPathways, miRBase, Tissues, etc. The skill follows that — adding more sources increases noise without much information gain on already-well-annotated programs.
