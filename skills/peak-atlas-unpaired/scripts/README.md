# scripts/

Unpaired-specific R helpers, reproduced faithfully from the DC_Dictionary
foundation layer (see each file's provenance header). Shared primitives
(iterative-overlap merge, normalize_width, support_voting, frip, blacklist) live
in `peak-atlas-framework/scripts/` — source those rather than reimplementing.

| Script | Functions | What it does |
|---|---|---|
| `label_transfer_strategy.R` | `add_gene_activity_bridge`, `transfer_labels`, `process_group_peaks`, `call_strategy_peaks` | GeneActivity bridge, staged R1/R2 CCA label transfer, and per-strategy (A/B/C) Signac peak calling + the shared 501bp/QC post-process. |
| `build_consensus.R` | `process_peaks_to_501bp`, `merge_replicates`, `build_external_consensus` | Harmonize external datasets to 501bp, merge replicates (>=50% support), union scaffold, per-peak dataset support, 0.25 consensus threshold, and the bias check. |
| `tier_stratification.R` | `assign_tiers`, `stratify_tier0`, `keep_tier0c` | Priority-order tier assignment over the A/B/C union and Tier 0 sub-stratification into 0a/0b/0c. |

Genome / path specifics (mouse build, `effective.genome.size`, macs3 path,
blacklist) are function arguments with sensible defaults — override for your
build. Source: `source("scripts/label_transfer_strategy.R")` etc.
