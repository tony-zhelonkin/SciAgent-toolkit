# support_voting.R - Cross-strategy support and support-weighted score boost
#
# Provenance:
#   calculate_strategy_support() reproduced faithfully from:
#     /data2/users/JCRLab/JBader/JBader_scHFD/01_scripts/R/peak_utils.R
#   The adjusted_score boost is reproduced from the pipeline:
#     /data2/users/JCRLab/JBader/JBader_scHFD/02_analysis/
#       S2_6b_atac_peak_recalling_multistrategies.R
#       (mcols$adjusted_score <- score * (1 + 0.5 * (n_strategies - 1)))
#
# CRITICAL: support must be computed on the per-strategy MERGED sets BEFORE the
# final iterative merge. Computing it after the final merge "bridges" every
# survivor onto all strategies and makes everything look fully supported.

suppressPackageStartupMessages({
  library(GenomicRanges)
})

#' Calculate cross-strategy support for peaks
#'
#' For each peak, count how many strategies produced an overlapping peak, and
#' record which strategies. Call this BEFORE the final merge.
#'
#' @param peaks GRanges object with peaks to annotate.
#' @param strategy_peaks Named list of GRanges (one per strategy, e.g.
#'   list(RNA=..., ATAC=..., WNN=..., CellType=...)). Names become the labels.
#' @return GRanges with `n_strategies` (integer) and `strategies` (comma string)
#'   columns added.
calculate_strategy_support <- function(peaks, strategy_peaks) {
  # Build support matrix: rows = peaks, cols = strategies
  support_matrix <- sapply(strategy_peaks, function(strat_peaks) {
    overlapsAny(peaks, strat_peaks)
  })

  mcols(peaks)$n_strategies <- rowSums(support_matrix)
  mcols(peaks)$strategies <- apply(support_matrix, 1, function(row) {
    paste(names(strategy_peaks)[row], collapse = ",")
  })

  return(peaks)
}

#' Apply the support-weighted score boost
#'
#' adjusted_score = score * (1 + 0.5 * (n_strategies - 1))
#'   1 strategy  -> x1.0
#'   2 strategies-> x1.5
#'   3 strategies-> x2.0
#'   4 strategies-> x2.5
#' The boosted score is the rank key for the final convergeClusterGRanges merge,
#' so multi-strategy consensus peaks win ties against single-strategy peaks.
#'
#' @param peaks GRanges with `score` and `n_strategies` columns.
#' @param score_col Character. Base score column (default "score").
#' @param boost Numeric. Per-extra-strategy boost factor (default 0.5).
#' @return GRanges with `adjusted_score` column added.
add_adjusted_score <- function(peaks, score_col = "score", boost = 0.5) {
  stopifnot(score_col %in% colnames(mcols(peaks)))
  stopifnot("n_strategies" %in% colnames(mcols(peaks)))
  mcols(peaks)$adjusted_score <- mcols(peaks)[[score_col]] *
    (1 + boost * (mcols(peaks)$n_strategies - 1))
  return(peaks)
}

# ----------------------------------------------------------------------------
# Canonical order (see references/support-voting.md):
#   1. Merge WITHIN each strategy:    strategy_merged[[s]] <- convergeClusterGRanges(...)
#   2. combined <- do.call(c, unname(strategy_merged))
#   3. combined <- calculate_strategy_support(combined, strategy_merged)  # BEFORE merge
#   4. combined <- add_adjusted_score(combined)
#   5. final <- convergeClusterGRanges(combined, by = "adjusted_score")
# ----------------------------------------------------------------------------
