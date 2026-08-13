# Cross-strategy support voting

`n_strategies` is the framework's confidence signal: for each candidate peak, how many *independent* groupings of the same cells discovered an overlapping peak. A peak found by RNA clusters AND ATAC clusters AND WNN AND cell-type annotation is far more likely to be real chromatin than a peak seen by only one. This file covers how support is computed, the one ordering rule that must not be broken, the score boost, and how the FILTER/TIER stages key off support.

## The strategies

The multiome regime uses four complementary groupings (see `peak-atlas-multiome`):

| Strategy | Grouping basis |
|---|---|
| RNA | clusters from RNA expression |
| ATAC | clusters from ATAC accessibility |
| WNN | clusters from weighted-nearest-neighbor (joint) |
| CellType | expert biological annotation |

Each strategy is peak-called and merged independently, producing one merged GRanges per strategy. The unpaired regime substitutes its own groupings (e.g. transferred-identity clusters) but the voting machinery is identical.

## Computing support — `calculate_strategy_support`

```r
combined <- do.call(c, unname(strategy_merged))   # all 4 merged sets stacked
combined <- calculate_strategy_support(combined, strategy_merged)
table(mcols(combined)$n_strategies)   # distribution over 1..4
```

It builds a peaks x strategies boolean overlap matrix and sets `n_strategies = rowSums(matrix)` plus a `strategies` label string. (Provenance: `peak_utils.R::calculate_strategy_support`.)

## CRITICAL: compute support BEFORE the final merge

This is the single most important ordering rule in the pipeline.

- **Right:** compute `n_strategies` on the per-strategy *merged* sets, THEN run the final iterative-overlap merge. Each peak's support reflects how many strategies independently found it.
- **Wrong:** run the final merge first, then compute support on the survivors. The final merge has already collapsed overlapping peaks into single winners that overlap peaks from *every* strategy — so every survivor "bridges" all four strategies and `n_strategies` reads 4 for essentially everything. The signal is destroyed (the "bridging" artifact).

Symptom of getting this wrong: the support distribution is almost entirely `n_strategies == 4`. If you see that, you computed support after the merge. (This was the v1→v2 fix in the source pipeline; see `README_peak_calling.md`.)

## The score boost

```r
adjusted_score = score * (1 + 0.5 * (n_strategies - 1))
```

| n_strategies | multiplier |
|---|---|
| 1 | x1.0 |
| 2 | x1.5 |
| 3 | x2.0 |
| 4 | x2.5 |

`adjusted_score` becomes the rank key for the final `convergeClusterGRanges` merge, so when a high-`score` single-strategy peak overlaps a slightly-lower-`score` four-strategy peak, the consensus peak wins. (Provenance: `S2_6b_atac_peak_recalling_multistrategies.R`; `add_adjusted_score` in `scripts/support_voting.R`.)

## How FILTER and TIER key off support

- **Metadata pre-filter (FILTER stage):** keep a peak if `n_strategies >= 2` OR it is a high-scoring singleton (`score > quantile(score, 0.85)`). Multi-strategy consensus passes automatically; lone peaks must be exceptionally strong. (See `references/rare-celltype-protection.md`.)
- **Rescue mechanism:** a peak that would be pruned for sparsity is rescued if it has FULL support (`n_strategies == 4`) AND global presence `>= 10` cells — replicated rare biology near the threshold boundary.
- **TIER:** the children assign confidence tiers directly from `n_strategies`. Treat it as the primary tier axis; raw score is a tiebreaker, not the headline.

## See also

- `scripts/support_voting.R` — `calculate_strategy_support`, `add_adjusted_score`.
- `references/iterative-overlap-merge.md` — where `adjusted_score` is consumed.
- `references/rare-celltype-protection.md` — the pre-filter and rescue in full.
