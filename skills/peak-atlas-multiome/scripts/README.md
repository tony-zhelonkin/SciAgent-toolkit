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
Set `PEAK_ATLAS_FRAMEWORK_SCRIPTS` (env) to the framework scripts dir before
sourcing `call_peaks_multistrategy.R`.
