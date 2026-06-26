# External validation: the PRIMARY correlation test and CellRanger overlap

The unpaired atlas earns trust through an external-correlation test that doubles
as the Tier-0c decision, plus a CellRanger overlap sanity check and a
strategy-comparison view. These run alongside the framework's biology battery
(FRiP retention, marker promoter retention, embedding —
`peak-atlas-framework/references/validation-battery.md`). Provenance:
`S6.2a_validation_primary.R`, `S6.2a_validate_cellranger_overlap.R`,
`S4_validation_strategy_comparison.R`.

## PRIMARY: pseudobulk accessibility vs external datasets

The idea: if the atlas captures real biology, a cell type's pseudobulk
accessibility profile should correlate with external datasets that share that
biology. The test makes that quantitative.

### 1. ATAC pseudobulk per group

For each ATAC cell type (or cell type × condition, with `min_cells = 20` in the
condition-separated mode), sum fragments over the atlas peaks and normalize
TPM-style:

```r
counts <- FeatureMatrix(fragments = Fragments(atac),
                        features = atlas, cells = cells_in_group)
pb     <- Matrix::rowSums(counts)
pb     <- (pb / sum(pb)) * 1e6        # TPM-like
```

### 2. External binary accessibility

For each external dataset, mark which atlas peaks it overlaps (1) or not (0):

```r
overlap <- findOverlaps(atlas, external_peaks, ignore.strand = TRUE)
acc <- rep(0, length(atlas)); acc[unique(queryHits(overlap))] <- 1
```

### 3. Correlate

Build a matrix of all ATAC profiles + all external profiles and take
`cor(..., method = "pearson")`. The full matrix is saved for figure
regeneration; the per-pair values go to a CSV sorted by `r`.

### 4. The named hypothesis gate

The decisive number is a **specific, pre-registered hypothesis**, not the whole
heatmap. In the source: a putative BATF3-independent cDC1 subtype (cDC1B) should
track an IRF8-KO dataset, because IRF8-KO mice lack the cDC1 lineage. The gate:

```
cDC1B  ~  IRF8-KO   Pearson r >= 0.4   -> the subtype's peaks are real biology
```

This is what decides Tier 0c (`references/tier-stratification.md`):

```
r >= 0.5   keep all Tier 0
r 0.4-0.5  keep 0a + 0b, monitor 0c
r <  0.4   keep only 0a; the risky single-strategy peaks were likely noise
```

Pick the hypothesis that your rare biology *predicts* — a dataset whose
perturbation should leave a measurable accessibility signature in the cell type
you are trying to validate.

## CellRanger overlap (baseline sanity)

CellRanger's default peak set is an aggregated pan-cell pseudobulk — broad but
not cell-type resolved. Checking how much of the de-novo atlas overlaps it tells
you whether you dropped baseline accessibility:

```
overlap = 100 * (# atlas peaks overlapping CellRanger) / (# atlas peaks)
>= 80%   CellRanger is redundant; excluding it was appropriate (no action)
50-80%   partial coverage; consider adding CellRanger-only peaks to fill gaps
<  50%   CRITICAL: missing major regions; add CellRanger as an input and re-assemble
```

High overlap is the desired outcome — it means the multi-strategy atlas already
contains CellRanger's landscape and then adds cell-type-specific peaks on top.

## Strategy comparison (what each lens contributed)

An UpSet plot or 3-way Venn of A vs B vs C peaks shows the division of labor:
the A∩B∩C core is high-confidence shared accessibility; C-only peaks are the
label-free hedge catching signal A/B missed; B-only peaks are candidate rare
biology to be gated by the PRIMARY test. It is descriptive, not a pass/fail
gate, but it is the quickest way to see whether C (the hedge) is actually
adding value.

## Pitfalls

| Symptom | Cause | Fix |
|---|---|---|
| PRIMARY r low for the right reason vs wrong reason | external dataset processed for a different build / chromosome style | harmonize external peaks first (`process_peaks_to_501bp`); confirm overlaps are non-zero |
| Everything correlates highly | comparing total accessibility, not a specific contrast | use the named hypothesis (the perturbation that should specifically affect the target type) |
| Atlas-CellRanger overlap < 50% | the atlas genuinely lost baseline regions, or peaks are not on the same build | verify build, then add CellRanger peaks as a fourth input and re-assemble |

## See also

- `references/tier-stratification.md` — how the PRIMARY r decides Tier 0c.
- `peak-atlas-framework/references/validation-battery.md` — FRiP / marker /
  embedding gates that run alongside this test.
- `chromvar-motif-accessibility` — TF-activity validation of the final atlas.
