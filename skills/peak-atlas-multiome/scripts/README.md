# scripts/

Multiome-specific R helpers, reproduced faithfully from the source pipeline
(see each file's provenance header). Shared primitives are SOURCED from
`peak-atlas-framework/scripts/` (iterative_overlap, support_voting,
normalize_width, blacklist, frip), not re-defined here.

| Script | Functions | What it does |
|---|---|---|
| `call_peaks_multistrategy.R` | `call_peaks_with_pseudoreps`, `call_peaks_by_strategy`, `reconcile_strategies` (+ `create_pseudobulk_bedpe`, `call_macs3_peaks`, `split_pseudoreplicates`, `filter_by_pseudoreps`, `read_narrowpeak`) | Per-group pseudobulk + 50/50 pseudo-rep split + MACS3 + reproducibility filter, then the four-strategy combine, support-before-merge, and final iterative-overlap merge (CALL -> MERGE -> TIER). |
| `apply_primary_rescue_filter.R` | `apply_primary_rescue_filter`, `calculate_adaptive_thresholds` | The adaptive Primary+Rescue prune (FILTER Phase 4): per-cell-type clamp(0.02*n, 15, 200) PRIMARY plus the 4-strategy RESCUE. |

Genome / MACS path / blacklist / chromosome set are function arguments with
mouse mm39 defaults. For human pass `macs_genome = "hs"`,
`std_chroms = paste0("chr", c(1:22, "X", "Y"))`, and an hg38 blacklist BED.
`call_peaks_multistrategy.R` resolves the framework primitives **relative to its
own file**, so the sibling-skill default (`../../peak-atlas-framework/scripts`)
works from any working directory. Set `PEAK_ATLAS_FRAMEWORK_SCRIPTS` (env) only
when the framework lives somewhere else. A script it cannot find is a **fatal
error**, not a warning — the alternative is a much later "could not find function
clusterGRanges" that points nowhere near the cause. If you have already sourced
the primitives by hand they are detected and the miss is downgraded to a message.
