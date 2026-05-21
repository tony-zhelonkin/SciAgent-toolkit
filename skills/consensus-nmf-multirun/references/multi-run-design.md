# Multi-run cNMF design — full ± QC, subset ± QC

The single-run cNMF result has three known failure modes. (1) Programs sensitive to QC: a noisy ribosomal program shows up because dying cells were not filtered. (2) Programs dominated by majority-celltype variance: the T-cell-specific state never separates because the all-cells embedding is dominated by myeloid-vs-lymphoid axis. (3) Programs unstable to K: choosing K=10 vs K=11 gives slightly different gene lists — which is the "real" answer?

The multi-run design addresses all three by running cNMF in *parallel* across orthogonal slices and merging programs that recur. Programs robust to all three perturbations land in the High-confidence tier.

## The 2×2 grid

| | Raw counts | QC-filtered counts |
|------|------------|---------------------|
| Full dataset | `cNMF_full` | `cNMF_fullQC` |
| Subset `<S>` | `cNMF_<S>` | `cNMF_<S>QC` |

For N subsets, the grid expands to 2×(N+1) variants (the +1 is "no subset"). The defaults the user can clear in one back-and-forth:

- 1 subset (none) × 1 QC variant (raw) = **1 run** (whole-dataset only, raw)
- N subsets × 2 QC variants = **2N runs** (full Anton-style consensus design)

The <ref-scrna> reference uses 6 runs: full ± QC, Th1 ± QC, Th17 ± QC.

## Why pair raw + QC

QC filtering removes cells that look broken (high MT%, low gene counts, doublets). The cells retained are higher quality; the programs cNMF finds in QC-filtered data are cleaner. *But* QC also removes biology — a stressed-cell population is exactly the kind of thing that can look like dying cells. Running both raw and QC and merging captures the high-quality programs *and* keeps a record of programs that disappear under filtering (those are the candidates for "is this real biology or QC artefact?").

In the merge stage, a program present in `cNMF_full` but absent from `cNMF_fullQC` is the canary. It might be:

- Real stress-response biology that QC over-filtered → expand QC thresholds
- Dying-cell contamination → keep filter, drop the program

The skill flags such asymmetries; the call is the user's.

## Why pair full + subset

A "broad" program from the full run (e.g., "all cycling cells") is often a coarse decomposition that the subset run can refine into "early-S cycling" + "G2-M cycling" within just the celltype of interest. After merging, the full-dataset program clusters with both subset programs and inherits a broader gene list; the subset programs stand alone in their high-confidence tier with focused gene lists.

For a single-celltype project (e.g., the <ref-scrna> T-cell reference), the full-vs-subset distinction is less important — the "subset" *is* the T cells. For a multicellular project (e.g., spleen), full + per-celltype is informative.

## Compute budget

cNMF on 100k cells with K=8..12 and n_iter=100 takes ~30–60 minutes on a 16-core machine. Six runs is 3–6 hours. The parameterised script `run_cnmf_subset.py` is a single call per variant; orchestrate via a shell loop or a Snakemake/nf-core wrapper if the variant count is large. (The skill itself does not assume a workflow runner.)

`TEST_MODE` (subsample 1000 cells, n_iter=10) brings a single run to ~2 minutes — useful for sanity-checking the pipeline before committing the real budget. **Never use TEST_MODE results downstream**; they only test the wiring.

## Why not "single run with K-selected post-hoc"

A common alternative: run cNMF once with K large (e.g., K=20), then post-hoc cluster the resulting programs to find the "real" K. This conflates two separate uncertainties — the K-stability one and the QC-stability one. The multi-run design treats them as orthogonal, which is what they are.

## Naming convention

Variant names are short, kebab-or-underscore, and prefix the program columns:

```
cNMF_full_P1        # full dataset raw, K=10, program 1
cNMF_full_P2
...
cNMF_fullQC_P1      # full dataset QC-filtered, K=10, program 1
...
cNMF_T_cells_P1     # subset 'T_cells' raw, K=8, program 1
cNMF_T_cellsQC_P1   # subset 'T_cells' QC-filtered, K=8, program 1
```

The `QC` suffix is on the *variant name*, immediately before `_P<i>`. Stage 5's prefix-matching sorts variant names by length descending so `cNMF_fullQC` is matched before `cNMF_full`.
