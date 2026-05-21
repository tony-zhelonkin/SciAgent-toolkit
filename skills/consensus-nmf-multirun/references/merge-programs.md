# Merge programs — hierarchical clustering, rank aggregation, classification

Stage 5 turns the cross-source correlation matrix into a merged-program table. Three logical steps: cluster, aggregate, classify.

## Step 1 — Hierarchical clustering on `1 − |r|` distance

```python
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import numpy as np

corr = pd.read_csv(correlation_csv, index_col=0).values
nan_mask = np.isnan(corr)
if nan_mask.any():
    corr[nan_mask] = 0   # NaN -> 0 distance contribution; programs with no shared cells are 'far'

dist = 1 - np.abs(corr)
np.fill_diagonal(dist, 0)
dist = (dist + dist.T) / 2          # symmetric (numerical safety)

condensed = squareform(dist, checks=False)
Z = linkage(condensed, method="ward")

# Cut at threshold. Reference uses a heuristic: cutoff = (1 - r_threshold) * Z[-1, 2] * 0.5
# That heuristic is a project-tunable knob; the shipped script exposes both 'absolute' and
# 'relative' cut modes. Default 'absolute' cuts at 1 - r_threshold (= 0.3 for r=0.7).
cutoff = 1 - correlation_threshold      # 0.3 for r=0.7
clusters = fcluster(Z, t=cutoff, criterion="distance")
```

Each `clusters[i]` is the cluster id for `program_names[i]`. Programs in the same cluster are merged.

**Why `1 - |r|` and not `1 - r`.** A program and its negative are biologically the same factor; absolute correlation captures that. Pure cNMF gives non-negative spectra so anti-correlation is rare in practice, but the absolute-value form is harmless when present.

**Why ward linkage.** Reference uses ward; centroid and complete also work. Ward minimises within-cluster variance — good for tight, biologically-coherent clusters.

## Step 2 — Rank aggregation per cluster

Per cluster:

- **Singleton** (one program in the cluster) → directly inherit top-100 genes; `consensus_score=1.0` for each gene.
- **Multi-program** (≥ 2 programs) → rank-aggregate top-50 genes per source.

Rank aggregation logic (port of <ref-scrna> `compute_consensus_genes`):

```python
from collections import defaultdict

n_programs = len(programs)
max_rank = TOP_GENES_PER_SOURCE   # 50

gene_stats = defaultdict(lambda: {"ranks": [], "sources": []})
for prog in programs:
    for rank, gene in enumerate(program_data[prog]["genes"][:max_rank]):
        gene_stats[gene]["ranks"].append(rank + 1)
        gene_stats[gene]["sources"].append(prog)

scored = []
for gene, st in gene_stats.items():
    n_sources = len(st["sources"])
    avg_rank = float(np.mean(st["ranks"]))
    rank_score     = (max_rank - avg_rank + 1) / max_rank        # higher rank in source -> closer to 1
    coverage_score = n_sources / n_programs                       # fraction of source programs hit
    consensus_score = rank_score * (1 + coverage_score)           # 1 + coverage so coverage doubles weight
    scored.append({"gene": gene, "consensus_score": consensus_score,
                   "n_sources": n_sources, "avg_rank": avg_rank})

scored.sort(key=lambda x: x["consensus_score"], reverse=True)
top_100 = scored[:TOP_GENES_TO_KEEP]                              # 100
```

The `(1 + coverage_score)` factor doubles the weight of genes that appear in *every* source program of the cluster — a robustness signal.

## Step 3 — Confidence tier and classification

### Confidence tier

```python
n_unique_runs = len({program_data[p]["run"] for p in programs})    # distinct runs, not programs
if   n_unique_runs >= 4: confidence = "High"
elif n_unique_runs >= 2: confidence = "Medium"
else:                    confidence = "Low"
```

`run` is the variant name (`full`, `full_qc`, `T_cells`, `T_cells_qc`, …). High-confidence programs survived ≥4 of the ≥4 variants in the design space.

### Classification (mouse panels — adapt for human)

Top 50 genes per merged program are scanned against fixed panels:

```python
# Mouse panels — all symbols are mouse-cased.
RIBOSOMAL_PATTERN = r"^Rp[sl]\d*[a-z]?$|^Rplp\d$|^Rpsa$"
MITO_PATTERN       = r"^mt-"
CELL_CYCLE_GENES = {
    "Top2a", "Mki67", "Cdk1", "Ccna2", "Ccnb1", "Ccnb2", "Cdc20",
    "Birc5", "Ube2c", "Nusap1", "Cenpf", "Tpx2", "Plk1", "Aurka",
    "Aurkb", "Bub1", "Bub1b", "Esco2", "Ncapd2", "Ncapg", "Smc2", "Smc4",
}
IMMEDIATE_EARLY_GENES = {
    "Fos", "Fosb", "Jun", "Junb", "Jund", "Egr1", "Egr2", "Egr3",
    "Nr4a1", "Nr4a2", "Nr4a3", "Ier2", "Ier3", "Atf3", "Btg2",
}
HOUSEKEEPING_GENES = {
    "Gapdh", "Actb", "Actg1", "Tuba1b", "Tubb5", "B2m", "Hsp90ab1",
    "Eef1a1", "Eef2", "Aldoa", "Pgk1", "Ldha", "Eno1", "Tpi1",
}
```

Classification rules (in order; first match wins):

| Condition (top-50 fraction) | Class |
|------------------------------|-------|
| `cell_cycle > 0.20` | `CellCycle` |
| `ribo > 0.20` | `Ribosomal` |
| `mt > 0.20` | `Mitochondrial` |
| `(ribo + mt + housekeeping) > 0.30` | `Technical` |
| `IEG > 0.15` | `ImmediateEarly` |
| else | `Biological` |

For human projects, replace symbols with capitalised versions (`MKI67`, `TOP2A`, `RPL\\d`, `MT-`, …). The shipped script accepts a `species` parameter that picks the right panel.

## Output

```
03_results/tables/programs/
├── merged_programs_v2.csv        # one row per merged program: name, confidence, class, top_100_genes, source_programs
├── program_annotations.csv       # canonical name + classification details (% ribo, % MT, % cell-cycle, ...)
└── gene_consensus_scores.csv     # one row per (program, gene): consensus_score, n_sources, avg_rank
```

The dendrogram from `Z` is saved to `03_results/plots/cnmf/program_dendrogram.png` (with the cluster cut line drawn). It is the single most useful diagnostic if the merge looks wrong.

## Sanity checks before downstream

- `merged_programs_v2['confidence'].value_counts()` — expect a mix of High/Medium/Low; if all Low, the multi-run design did not converge on shared programs (threshold too strict, or runs too disparate).
- `merged_programs_v2['class'].value_counts()` — expect Biological as the majority; a project where most programs are Technical/Ribosomal indicates QC-filter failure.
- Inspect `program_dendrogram.png` — clusters should look tight (short within-cluster height) and well-separated (long inter-cluster height).
