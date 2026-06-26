# Iterative-overlap merge: the algorithm and its embedding here

The MERGE stage resolves a pile of overlapping fixed-width peaks into a non-overlapping consensus set using **iterative-overlap winner selection** (Corces & Granja, Science 2018; ArchR). This file describes the algorithm in prose, its relationship to the standalone `iterative-peak-merging` skill, and the SPM/rule reproducibility variant.

## The algorithm (summit-faithful winner-take-all)

Given many overlapping fixed-width peaks, each with a score:

1. **Reduce** the peaks to define overlap *clusters* — maximal stretches of mutually-overlapping peaks (`GenomicRanges::reduce(..., min.gapwidth=0L)`). Reduce is used only to LABEL clusters, never to emit coordinates.
2. **Rank** peaks within each cluster by score (descending).
3. **Take the best peak** in each cluster — the actual highest-score peak, with its real summit-centered coordinates.
4. **Remove** that winner and everything still overlapping it from the pool.
5. **Repeat** on the remaining pool until no overlaps survive.

This is implemented as `clusterGRanges` (one pass: label clusters, optionally keep the best per cluster) wrapped by `convergeClusterGRanges` (loop until convergence) in `scripts/iterative_overlap.R`.

### Why not `reduce()` / `bedtools merge`?

A plain reduce/merge would emit the *union interval* of overlapping peaks — a coordinate no individual peak actually has, effectively a midpoint smear. That destroys summit fidelity and breaks fixed-width assumptions. Iterative overlap instead keeps the **actual best peak** and discards its overlappers, so every survivor is a real 501bp summit-centered call.

It also avoids **daisy-chaining**: with a transitive merge, peak A (overlaps B) and peak C (overlaps B) get fused into one giant region even though A and C never overlap. Winner-take-all removal cannot daisy-chain.

## Support-weighted embedding here

Standalone, you rank by `score`. In this framework the rank key for the FINAL merge is the support-boosted `adjusted_score = score * (1 + 0.5*(n_strategies-1))` (see `references/support-voting.md`), so a peak corroborated by multiple independent strategies beats a single-strategy peak in a tie. The within-strategy merges still rank by raw `score` — support is only known once each strategy has been merged.

Canonical order:
```r
# 1. merge within each strategy (rank by score)
strategy_merged[[s]] <- convergeClusterGRanges(strat_peaks[[s]], by = "score")
# 2. combine, then compute support BEFORE the final merge
combined <- do.call(c, unname(strategy_merged))
combined <- calculate_strategy_support(combined, strategy_merged)
combined <- add_adjusted_score(combined)
# 3. final winner-take-all on the boosted score
final <- convergeClusterGRanges(combined, by = "adjusted_score")
```

## Relation to `iterative-peak-merging`

| | `iterative-peak-merging` (sibling) | Here (embedded) |
|---|---|---|
| Input | MACS `_summits.bed` files in a directory + metadata | Per-strategy merged GRanges in memory |
| Design | Within-study replicate/group design | Cross-strategy single-cell consensus |
| Rank key | `score` (CPM), reproducibility rule | support-boosted `adjusted_score` |
| Role | Standalone primitive | One stage of CALL→MERGE→FILTER→TIER→VALIDATE |

They are the **same core algorithm**. Do not reimplement the primitive in a child skill — for a from-BED within-study merge route to `iterative-peak-merging`; for the in-pipeline support-aware merge use `scripts/iterative_overlap.R` here.

## SPM / reproducibility-rule variant

`createIterativeOverlapPeakSet.R` (`extendedPeakSet`) is the bulk/pseudobulk reproducibility path. After extending summits to fixed width and removing cliffed/blacklist peaks, it:
- merges within group with `convergeClusterGRanges(by="score")`,
- rescales score to **score-per-million** via `edgeR::cpm`,
- keeps peaks above `scorePerMillion` (default `spm = 5`),
- applies a **reproducibility rule** across samples (default `rule = "2"`; `"(n+1)/2"` = majority): a peak must clear SPM in at least `minSamples` samples,
- re-normalizes SPM and does a final cross-group `convergeClusterGRanges`.

Defaults: `extend = 250` (→ 501bp), `spm = 5`, `rule` evaluated as `(n+1)/2`. Use this path when you have explicit replicates and want a hard reproducibility gate; use the support-voting path (above) for single-cell multi-strategy consensus.

## See also

- `scripts/iterative_overlap.R` — `clusterGRanges`, `convergeClusterGRanges`.
- `references/support-voting.md` — the boost and the compute-before-merge rule.
