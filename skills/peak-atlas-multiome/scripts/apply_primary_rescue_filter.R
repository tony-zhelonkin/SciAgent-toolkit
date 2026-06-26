# apply_primary_rescue_filter.R
#   The adaptive "Primary + Rescue" peak filter (FILTER stage, Phase 4).
#
# Provenance (reproduced faithfully):
#   apply_primary_rescue_filter() and calculate_adaptive_thresholds() from
#     /data2/users/JCRLab/JBader/JBader_scHFD/01_scripts/R/peak_utils.R
#   Wiring (config-driven call) from
#     /data2/users/JCRLab/JBader/JBader_scHFD/02_analysis/S2_6b_filter_peak_atlas.R
#
# Design (see references/primary-rescue-filter.md and the framework's
#   peak-atlas-framework/references/rare-celltype-protection.md):
#   PRIMARY: per cell-type adaptive threshold max(floor, min(rate*n, cap)).
#            Keep a peak if it clears the bar in ANY ONE cluster, so a rare
#            type faces a low absolute bar and an abundant type a higher
#            relative bar.
#   RESCUE : keep peaks with n_strategies == 4 AND >= rescue_global_min cells
#            globally (replicated rare biology near the boundary).
#   final_keep = primary_keep | rescue_keep.
#
# Default parameters are the source values (floor=15, rate=0.02, cap=200,
# rescue_global_min=10); they are arguments so you can relax/tighten per dataset.

suppressPackageStartupMessages({
  library(Matrix)
})

#' Adaptive per-cluster cell-count threshold: clamp(rate * n, floor, cap).
#'
#' @param cluster_sizes Named integer vector of cells per cluster (e.g.
#'   table(cluster_labels)).
#' @param floor Absolute minimum cells (protects rare types). Default 15.
#' @param rate  Within-cluster prevalence for abundant types. Default 0.02 (2%).
#' @param cap   Maximum cells (stops abundant types over-filtering). Default 200.
#' @return Named numeric vector of thresholds, aligned to names(cluster_sizes).
#'
#' @examples
#'   calculate_adaptive_thresholds(c(HSC = 147, pDC = 180, Neutrophil = 4500))
#'   # HSC: 15 (10.2%), pDC: 15 (8.3%), Neutrophil: 90 (2.0%)
calculate_adaptive_thresholds <- function(cluster_sizes,
                                          floor = 15, rate = 0.02, cap = 200) {
  thresholds <- pmax(floor, pmin(rate * as.numeric(cluster_sizes), cap))
  names(thresholds) <- names(cluster_sizes)
  thresholds
}

#' Primary + Rescue peak filter.
#'
#' @param peaks GRanges with an `n_strategies` metadata column. Row order MUST
#'   match the rows of `counts`.
#' @param counts Sparse peaks-by-cells count matrix (e.g. a subsampled
#'   Signac::FeatureMatrix). nrow(counts) == length(peaks).
#' @param cluster_labels Per-cell cluster/cell-type labels, length ncol(counts).
#' @param floor,rate,cap PRIMARY adaptive-threshold parameters (see above).
#' @param rescue_global_min Minimum global cells for the RESCUE arm. Default 10.
#' @param verbose Print progress. Default TRUE.
#' @return Logical vector (length == length(peaks)) of peaks to KEEP.
apply_primary_rescue_filter <- function(peaks, counts, cluster_labels,
                                        floor = 15, rate = 0.02, cap = 200,
                                        rescue_global_min = 10, verbose = TRUE) {
  stopifnot(nrow(counts) == length(peaks))
  stopifnot(ncol(counts) == length(cluster_labels))
  stopifnot("n_strategies" %in% colnames(GenomicRanges::mcols(peaks)))

  if (verbose) {
    message("\n=== Primary + Rescue Filtration ===")
    message("Parameters: floor=", floor, ", rate=", rate,
            ", cap=", cap, ", rescue>=", rescue_global_min)
  }

  # Adaptive thresholds per cluster.
  cluster_sizes <- table(cluster_labels)
  thresholds <- calculate_adaptive_thresholds(cluster_sizes, floor, rate, cap)
  if (verbose) {
    message("\nAdaptive thresholds per cluster:")
    for (ct in names(thresholds)) {
      pct <- 100 * thresholds[ct] / cluster_sizes[ct]
      message(sprintf("  %s (n=%d): %d cells (%.1f%%)",
                      ct, cluster_sizes[ct], round(thresholds[ct]), pct))
    }
  }

  # PRIMARY: per-cluster count of cells with the peak accessible (count > 0).
  clusters <- unique(cluster_labels)
  cluster_counts <- matrix(0, nrow = nrow(counts), ncol = length(clusters))
  colnames(cluster_counts) <- clusters
  if (verbose) message("\nCalculating per-cluster accessibility...")
  for (i in seq_along(clusters)) {
    ct_idx <- which(cluster_labels == clusters[i])
    cluster_counts[, i] <- Matrix::rowSums(counts[, ct_idx, drop = FALSE] > 0)
  }

  # PRIMARY KEEP: clears its cluster's threshold in ANY ONE cluster.
  primary_keep <- apply(cluster_counts, 1, function(row) {
    any(row >= thresholds[colnames(cluster_counts)])
  })
  if (verbose) message("  PRIMARY filter: ",
                       format(sum(primary_keep), big.mark = ","), " peaks retained")

  # RESCUE: full 4-strategy support + minimal global presence.
  peak_presence_global <- Matrix::rowSums(counts > 0)
  rescue_keep <- (peaks$n_strategies == 4) & (peak_presence_global >= rescue_global_min)
  if (verbose) {
    message("  RESCUE mechanism: ", format(sum(rescue_keep), big.mark = ","),
            " peaks (4-strategy + >=", rescue_global_min, " cells)")
    message("    Added by RESCUE: ",
            format(sum(rescue_keep & !primary_keep), big.mark = ","), " not in PRIMARY")
  }

  final_keep <- primary_keep | rescue_keep
  if (verbose) {
    message("\nFINAL: ", format(sum(final_keep), big.mark = ","), " peaks retained (",
            format(100 * sum(final_keep) / length(peaks), digits = 2), "%)")
  }
  final_keep
}

# ------------------------------------------------------------------------------
# Usage (after the subsampled FeatureMatrix in Phase 3):
#   keep <- apply_primary_rescue_filter(
#             peaks = peaks_prefiltered,        # GRanges with n_strategies
#             counts = counts_sampled,          # peaks x sampled-cells
#             cluster_labels = obj$refined_cell_type[colnames(counts_sampled)])
#   peaks_filtered <- peaks_prefiltered[keep]
# ------------------------------------------------------------------------------
