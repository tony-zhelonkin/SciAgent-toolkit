# iterative_overlap.R - Iterative-overlap winner selection (Corces & Granja 2018)
#
# Provenance:
#   Reproduced faithfully from:
#     /data2/users/JCRLab/JBader/JBader_scHFD/01_scripts/R/peak_utils.R
#       (clusterGRanges, convergeClusterGRanges)
#   The SPM/rule reproducibility variant is documented in:
#     /data2/users/JCRLab/JBader/JBader_scHFD/01_scripts/R/createIterativeOverlapPeakSet.R
#       (extendedPeakSet: defaults extend=250, spm=5, rule="2" => (n+1)/2)
#   Originally from Corces MR, Granja JM, et al. Science 2018;362(6413):eaav1898.
#
# Why this method (vs reduce()/bedtools merge):
#   It is SUMMIT-FAITHFUL. It picks the ACTUAL highest-score peak per overlap
#   cluster and removes everything that overlaps it, rather than collapsing
#   overlapping peaks to a midpoint. Embedded in the peak-atlas pipeline it is
#   also SUPPORT-WEIGHTED: pass by = "adjusted_score" to rank on the
#   cross-strategy support boost (see support_voting.R).

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(GenomeInfoDb)
})

#' Cluster overlapping peaks and optionally filter to best per cluster
#'
#' Finds overlapping peaks, assigns them to clusters (via reduce), and
#' optionally selects the best peak per cluster by a score column.
#'
#' @param gr GRanges object with peaks
#' @param filter Logical. If TRUE, return only best peak per cluster
#' @param by Character. Column name to use for scoring (default "score")
#' @param decreasing Logical. If TRUE, higher scores are better (default TRUE)
#' @param verbose Logical. Print progress (default FALSE)
#' @return GRanges object (filtered if filter=TRUE)
clusterGRanges <- function(gr, filter = TRUE, by = "score",
                           decreasing = TRUE, verbose = FALSE) {
  # Sort peaks
  gr <- sort(sortSeqlevels(gr))

  # Find overlapping regions (reduce to get cluster boundaries)
  r <- GenomicRanges::reduce(gr, min.gapwidth = 0L, ignore.strand = TRUE)

  # Assign cluster IDs based on overlap
  o <- findOverlaps(gr, r)
  mcols(gr)$cluster <- subjectHits(o)

  if (verbose) {
    n_overlaps <- length(gr) - max(subjectHits(o))
    message(sprintf("  Found %d overlapping peaks in %d clusters",
                    n_overlaps, max(subjectHits(o))))
  }

  if (filter) {
    if (by %in% colnames(mcols(gr))) {
      if (verbose) message(sprintf("  Filtering by %s...", by))
      # Sort by score and take best per cluster
      gr <- gr[order(mcols(gr)[[by]], decreasing = decreasing), ]
      gr <- gr[!duplicated(mcols(gr)$cluster), ]
      gr <- sort(sortSeqlevels(gr))
    } else {
      # Fall back to order-based deduplication
      if (verbose) message("  Filtering by order...")
      gr <- gr[!duplicated(mcols(gr)$cluster), ]
    }
    mcols(gr)$cluster <- NULL
  }

  return(gr)
}

#' Iteratively resolve overlapping peaks by picking the winner per cluster
#'
#' Repeatedly applies clusterGRanges until no overlaps remain. Preserves
#' summit fidelity - picks the ACTUAL best peak, not a midpoint.
#'
#' @param gr GRanges object with peaks (must have the `by` score column)
#' @param by Character. Column name to use for scoring (default "score")
#' @param decreasing Logical. If TRUE, higher scores are better (default TRUE)
#' @param verbose Logical. Print progress (default FALSE)
#' @return GRanges object with no overlapping peaks
convergeClusterGRanges <- function(gr, by = "score",
                                   decreasing = TRUE, verbose = FALSE) {
  stopifnot(by %in% colnames(mcols(gr)))

  i <- 0
  gr_initial <- gr

  if (verbose) message("  Starting iterative convergence...")

  while (length(gr_initial) > 0) {
    i <- i + 1
    if (verbose && i %% 10 == 1) message("    Iteration ", i, "...")

    # Get best peak from each overlap cluster
    gr_clustered <- clusterGRanges(
      gr = gr_initial,
      filter = TRUE,
      by = by,
      decreasing = decreasing,
      verbose = FALSE
    )

    # Remove selected peaks (and everything overlapping them) from the pool
    gr_initial <- subsetByOverlaps(gr_initial, gr_clustered, invert = TRUE)

    # Accumulate results
    if (i == 1) {
      gr_all <- gr_clustered
    } else {
      gr_all <- c(gr_all, gr_clustered)
    }
  }

  if (verbose) message("  Converged after ", i, " iterations")
  return(sort(sortSeqlevels(gr_all)))
}

# ----------------------------------------------------------------------------
# Usage (embedded, support-weighted final merge):
#   source("scripts/support_voting.R")  # adds adjusted_score / n_strategies
#   final <- convergeClusterGRanges(combined, by = "adjusted_score",
#                                   decreasing = TRUE, verbose = TRUE)
#
# Usage (within-strategy merge before support is known):
#   strat_merged <- convergeClusterGRanges(strat_peaks, by = "score")
# ----------------------------------------------------------------------------
