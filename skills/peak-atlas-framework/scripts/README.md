# scripts/

Shared R helpers for the peak-atlas-* family, reproduced faithfully from the
source pipeline (see each file's provenance header). Children link here rather
than reimplementing.

| Script | Functions | What it does |
|---|---|---|
| `iterative_overlap.R` | `clusterGRanges`, `convergeClusterGRanges` | Summit-faithful iterative-overlap winner selection (Corces & Granja 2018). |
| `normalize_width.R` | `normalize_to_501bp` | Re-center peaks on the MACS3 summit and resize to a fixed width (501bp at extend=250). |
| `support_voting.R` | `calculate_strategy_support`, `add_adjusted_score` | Cross-strategy `n_strategies` support and the `adjusted_score` boost. Compute support BEFORE the final merge. |
| `frip.R` | `calculate_frip` | Fraction of reads in peaks, per cell, via Signac FeatureMatrix. |
| `blacklist.R` | `load_blacklist`, `remove_blacklist_peaks` | Load an ENCODE blacklist (any build via `bed_path`) and invert-overlap to drop contaminated peaks. |

Genome / blacklist specifics are function arguments with sensible defaults, not
hard-coded paths. Source: `source("scripts/iterative_overlap.R")` etc.
