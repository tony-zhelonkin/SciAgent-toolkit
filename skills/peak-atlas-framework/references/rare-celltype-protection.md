# Rare cell-type protection — the first principle

The whole point of a multi-strategy consensus atlas is to *not* lose the peaks that define rare populations (HSCs, pDCs, cycling progenitors). A naive global filter — "keep peaks present in at least N cells" — silently deletes rare biology, because a rare type simply does not have N cells. This file is the design that prevents that. Provenance: `peak_utils.R` (`calculate_adaptive_thresholds`, `apply_primary_rescue_filter`) and `02_analysis/config/peak_filtration_config.yaml`.

## The principle

> Don't chase a target peak count if it means losing biology. Chase FRiP and marker coverage. Give a rare type a low absolute bar and an abundant type a higher relative bar, and keep a peak if it clears the bar in ANY one cell type.

## Adaptive per-cluster threshold

For each cluster (cell type), the minimum number of accessible cells required is:

```
threshold(cluster) = clamp(prevalence_rate * n,  floor,  cap)
                   = max(15,  min(0.02 * n,  200))
```

with `floor = 15`, `prevalence_rate = 0.02` (2%), `cap = 200`.

| Cluster | n cells | threshold | as % |
|---|---|---|---|
| HSC | 147 | 15 | 10.2% |
| pDC | 180 | 15 | 8.3% |
| Neutrophil | 4500 | 90 | 2.0% |

The **floor** (15) keeps a minimum statistical footing for rare types (standard error ≈ 0.13) without demanding a percentage they can never reach. The **cap** (200) stops abundant clusters from imposing an absurd absolute bar that would over-filter common peaks. A peak is kept by the PRIMARY filter if it clears its cluster's threshold in **any one** cluster. (`calculate_adaptive_thresholds`; `apply_primary_rescue_filter` PRIMARY arm.)

## Stratified cell sampling

Quantification (the cell-by-peak count) is expensive, so the source subsamples cells — but stratified by type so rare types are not sampled away:

```
cells_per_type = max(min_cells_per_type, sample_fraction * n_type)
# min_cells_per_type = 100,  sample_fraction = 0.10
```

Abundant types contribute ~10%; rare types contribute at least `min_cells_per_type` (capped at their actual count). This keeps rare-type signal representative in the quantified matrix. (`cell_sampling` in the config.)

## Rescue mechanism

A peak that *fails* the adaptive PRIMARY filter is still rescued if it is high-confidence replicated rarity:

```
rescue_keep = (n_strategies == 4) & (global_presence >= 10 cells)
final_keep  = primary_keep | rescue_keep
```

Found by ALL FOUR strategies and present in at least 10 cells globally → keep it even though it never cleared a per-cluster bar. This protects peaks that sit right at the threshold boundary but are corroborated by every independent grouping. (`apply_primary_rescue_filter` RESCUE arm; `sparsity_filter_v2.rescue_mechanism` in the config: `require_n_strategies: 4`, `min_global_cells: 10`.)

## Metadata pre-filter

Before any quantification, cheaply drop obvious noise using only the peak metadata:

```
keep if  n_strategies >= 2
     OR  score > quantile(score, 0.85)   # high-scoring singleton rescue
```

Multi-strategy consensus passes automatically; lone single-strategy peaks must be in the top 15% by score to survive the pre-filter. (`metadata_filter.score_percentile: 0.85`.)

## Putting it together (FILTER stage)

1. **Metadata pre-filter** — `n_strategies >= 2` OR top-15% singleton.
2. **Stratified sampling** — `max(100, 0.10*n)` cells per type, then quantify (FeatureMatrix).
3. **Adaptive PRIMARY** — keep if cell count clears `max(15, min(0.02*n, 200))` in ANY cluster.
4. **RESCUE** — re-add 4-strategy peaks with `>= 10` global cells.
5. Validate retention (FRiP, marker coverage) — `references/validation-battery.md`.

Tunables live in `assets/peak_filtration_config.yaml`. The whole battery is parameterized so you can relax (`relaxed_min_cells`) or tighten (`strict_min_cells`) without touching code.

## See also

- `references/support-voting.md` — where `n_strategies` comes from.
- `references/validation-battery.md` — confirming biology survived.
- `assets/peak_filtration_config.yaml` — the tunable template.
