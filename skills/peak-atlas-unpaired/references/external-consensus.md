# External-dataset consensus

In the unpaired regime, published ATAC datasets from the same biological system
are a **peak source**, not just an annotation layer. Folding them into a
consensus adds reproducible coordinates the de-novo strategies might have
missed, and gives you the binary accessibility profiles the PRIMARY validation
test needs. Provenance: `S4.1a_process_external_datasets.R` (harmonize),
`config.R::merge_replicates`, `S4.1b_build_consensus.R` (consensus + bias).

## 1. Harmonize each dataset to 501bp

External peaks come in many widths, chromosome styles, and summit conventions.
Before they can be compared, every dataset is put through the same CALL-stage
pipeline as the de-novo peaks:

1. coerce chromosome names to UCSC style (`chr1`, `chrM`, …),
2. center on the summit if present (`summit` or narrowPeak `peak` column),
   else on the peak center,
3. `resize(width = 501, fix = "center")`,
4. keep standard chromosomes only,
5. drop blacklist overlaps,
6. drop N-rich peaks (N fraction `>= 0.1`).

This makes per-peak overlap comparisons fair (equal widths) and removes the same
artifacts removed from de-novo peaks. Code: `build_consensus.R::process_peaks_to_501bp`.

## 2. Merge replicates first (>= 50% support)

Counting biological/technical replicates as independent datasets inflates the
apparent support of whatever those replicates happen to contain
(pseudo-replication). Collapse each replicate group first: build a union
scaffold (`reduce(min.gapwidth = 50)`) and keep peaks present in at least half
the replicates:

```r
min_support <- ceiling(length(rep_ids) * 0.5)   # 2 reps -> 100% agreement
consensus   <- scaffold[n_support >= min_support]
```

In the source this turned 35 replicate-level datasets into ~22 independent ones.
Code: `build_consensus.R::merge_replicates`.

## 3. Union scaffold across datasets

Stack all independent datasets and reduce with a wider gap so nearby peaks join:

```r
union_peaks <- GenomicRanges::reduce(do.call(c, unname(external_peaks)),
                                     min.gapwidth = 100)
```

**Two gap widths — do not conflate them.** `gapwidth = 100` is used HERE, for the
union *across independent datasets*, because different pipelines place peaks with
looser positional agreement. The tighter `gapwidth = 50` is used for the other,
WITHIN-study merges: collapsing replicate groups (step 2 above) and merging peaks
across groups *within a single A/B/C strategy* (see `references/abc-strategies.md`).
Using 50 for the cross-dataset union under-merges near-duplicate peaks from
different studies; using 100 within a strategy over-merges genuinely distinct
adjacent peaks.

## 4. Per-peak support matrix

For each union peak, record which datasets overlap it:

```r
support_matrix <- sapply(external_peaks, function(p) overlapsAny(union_peaks, p))
n_per_peak     <- rowSums(support_matrix)
```

## 5. Consensus threshold 0.25

Keep peaks present in at least 25% of independent datasets:

```r
min_datasets    <- ceiling(length(external_peaks) * consensus_threshold)  # 0.25
consensus_peaks <- union_peaks[n_per_peak >= min_datasets]
```

25% balances breadth (don't demand near-universal agreement, which would lose
lineage-specific regions) against noise (don't keep a peak seen by one dataset).
Tune up for a stricter atlas, down for broader coverage.

## 6. Bias check (the dataset can't drive the consensus)

A consensus dominated by one big atlas is not a consensus. Compute each
dataset's contribution — the fraction of consensus peaks it overlaps — and gate:

```r
contribution(d) = 100 * (# consensus peaks overlapping d) / (# consensus peaks)
warn if max(contribution) > 40%
FAIL if max(contribution) > 50%
```

A `> 50%` breach usually means replicates were not merged (one study's many
replicates dominate) or one atlas is simply much larger than the rest. The
runnable gate is `checks/check_consensus_bias.R`. Code:
`build_consensus.R::build_external_consensus`.

## Annotation-only external sets

Some external sets cannot contribute coordinates — e.g. an aging atlas processed
for a different genome build, or peak tables without usable positions for your
assembly. Keep those OUT of the coordinate union (they would add noise or
mislifted intervals) but still use them at validation: they can flag whether a
de-novo peak overlaps known biology even when they cannot donate coordinates.

## Pitfalls

| Symptom | Cause | Fix |
|---|---|---|
| One dataset contributes > 50% | replicates counted as independent | run `merge_replicates` (>=50% support) before consensus |
| Consensus peaks not 501bp | external peaks not harmonized first | always `process_peaks_to_501bp` each dataset |
| Almost no consensus peaks | threshold too high for heterogeneous datasets | lower `consensus_threshold` or check chromosome-style mismatch (UCSC vs Ensembl) dropped overlaps |
| Near-duplicate peaks from different studies kept separate | used the within-study `gapwidth = 50` for the cross-dataset union | use `gapwidth = 100` for the cross-dataset union; reserve `50` for replicate and within-strategy merges |

## See also

- `references/external-validation.md` — using these datasets for the PRIMARY
  correlation test.
- `checks/check_consensus_bias.R`, `scripts/build_consensus.R`.
