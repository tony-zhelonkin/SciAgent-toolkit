# build_consensus.R - External-dataset consensus peak source for the unpaired atlas.
#
# Provenance (reproduced faithfully from the DC_Dictionary foundation layer):
#   process_peaks_to_501bp() -> foundation/S4.1a_process_external_datasets.R
#   merge_replicates()       -> 00_main/pipeline/config.R (>=50% replicate support)
#   union / support matrix / 0.25 threshold / bias check
#                            -> foundation/S4.1b_build_consensus.R
#
# External ATAC datasets from the same biological system contribute COORDINATES
# (a peak source), not merely annotation overlap. This builds a consensus from
# them: harmonize to 501bp -> optionally merge replicates -> union scaffold ->
# per-peak dataset support -> keep peaks in >= 25% of datasets -> bias check.
#
# Genome/path specifics (mouse mm10, blacklist) are function arguments.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(IRanges)
})

# ---------------------------------------------------------------------------
# Harmonize one external dataset to uniform 501bp summit-centered peaks
# (provenance: S4.1a_process_external_datasets.R::process_peaks_to_501bp)
# ---------------------------------------------------------------------------

#' Harmonize external peaks to 501bp with the standard CALL-stage filters
#'
#' @param peaks GRanges of external peaks (summit in mcols$summit or mcols$peak).
#' @param dataset_name Label stored in mcols$dataset.
#' @param blacklist GRanges blacklist (empty = skip).
#' @param bsgenome BSgenome for the N-content filter.
#' @param std_chr Standard chromosomes (default mouse).
#' @param use_summit If TRUE and a summit column exists, center on the summit.
#' @return GRanges of 501bp peaks on standard chromosomes, blacklist/N filtered.
process_peaks_to_501bp <- function(peaks, dataset_name, blacklist, bsgenome,
                                    std_chr = paste0("chr", c(1:19, "X", "Y")),
                                    use_summit = TRUE) {
  # Coerce to UCSC chromosome style if naked numerics/X/Y/MT
  if (all(grepl("^\\d+$|^X$|^Y$|^MT$", seqnames(peaks)))) {
    seqlevels(peaks) <- paste0("chr", seqlevels(peaks))
    seqlevels(peaks)[seqlevels(peaks) == "chrMT"] <- "chrM"
  }
  GenomeInfoDb::seqlevelsStyle(peaks) <- "UCSC"

  if (use_summit && "summit" %in% names(mcols(peaks))) {
    s <- start(peaks) + mcols(peaks)$summit
    peaks <- GRanges(seqnames(peaks), IRanges(start = s, width = 1), strand = "*")
  } else if (use_summit && "peak" %in% names(mcols(peaks))) {
    s <- start(peaks) + mcols(peaks)$peak
    peaks <- GRanges(seqnames(peaks), IRanges(start = s, width = 1), strand = "*")
  }

  peaks <- resize(peaks, width = 501, fix = "center")
  present_chr <- intersect(seqlevels(peaks), std_chr)
  peaks <- keepSeqlevels(peaks, present_chr, pruning.mode = "coarse")
  if (length(blacklist) > 0) peaks <- peaks[!overlapsAny(peaks, blacklist)]

  peak_seqs <- Biostrings::getSeq(bsgenome, peaks)
  n_content <- Biostrings::letterFrequency(peak_seqs, "N") / 501
  peaks <- peaks[n_content < 0.1]

  mcols(peaks)$dataset <- dataset_name
  peaks
}

# ---------------------------------------------------------------------------
# Merge biological/technical replicates within a group (>= 50% support)
# (provenance: config.R::merge_replicates)
# ---------------------------------------------------------------------------

#' Collapse replicate datasets into one consensus per group
#'
#' Corrects pseudo-replication BEFORE counting datasets: within each replicate
#' group, build a union scaffold then keep peaks present in >= 50% of replicates.
#'
#' @param peaks_list Named list of GRanges (one per replicate-level dataset).
#' @param replicate_groups Named list mapping group -> character vector of
#'   replicate ids in peaks_list.
#' @return list(merged_peaks = named list of group-level GRanges, merge_stats).
merge_replicates <- function(peaks_list, replicate_groups) {
  merged <- list()
  stats <- list()
  for (group_name in names(replicate_groups)) {
    rep_ids <- replicate_groups[[group_name]]
    group_peaks <- peaks_list[rep_ids]
    scaffold <- GenomicRanges::reduce(do.call(c, unname(group_peaks)),
                                      min.gapwidth = 50)
    support <- rowSums(sapply(group_peaks, function(pk) overlapsAny(scaffold, pk)))
    min_support <- ceiling(length(rep_ids) * 0.5)   # >= 50% of replicates
    merged[[group_name]] <- scaffold[support >= min_support]
    stats[[group_name]] <- c(n_replicates = length(rep_ids),
                             n_union = length(scaffold),
                             n_consensus = length(merged[[group_name]]))
  }
  list(merged_peaks = merged, merge_stats = do.call(rbind, stats))
}

# ---------------------------------------------------------------------------
# Build the external consensus (union -> support -> 0.25 threshold -> bias)
# (provenance: S4.1b_build_consensus.R)
# ---------------------------------------------------------------------------

#' Build a consensus peak set from a list of (replicate-merged) external datasets
#'
#' @param external_peaks Named list of GRanges (one per INDEPENDENT dataset;
#'   run merge_replicates() first if you have replicates).
#' @param consensus_threshold Fraction of datasets a peak must appear in
#'   (default 0.25 -> keep peaks present in >= 25% of datasets).
#' @param gap_width reduce() min.gapwidth for the union scaffold (default 100).
#' @param warn_frac,fail_frac Bias thresholds on max single-dataset contribution
#'   (default warn 0.40, fail 0.50). fail_frac breach raises an error.
#' @return list(consensus_peaks, support_matrix, dataset_contribution,
#'   max_contribution, min_datasets).
build_external_consensus <- function(external_peaks, consensus_threshold = 0.25,
                                      gap_width = 100,
                                      warn_frac = 0.40, fail_frac = 0.50) {
  external_peaks <- external_peaks[sapply(external_peaks, length) > 0]
  n_datasets <- length(external_peaks)
  if (n_datasets == 0) stop("No non-empty external datasets supplied")

  # Union scaffold across all datasets
  union_peaks <- GenomicRanges::reduce(
    do.call(c, unname(external_peaks)), min.gapwidth = gap_width)

  # Per-peak dataset support matrix
  support_matrix <- sapply(external_peaks, function(p) overlapsAny(union_peaks, p))
  n_per_peak <- rowSums(support_matrix)
  mcols(union_peaks)$n_datasets <- n_per_peak
  mcols(union_peaks)$pct_datasets <- round(100 * n_per_peak / n_datasets, 1)

  # Consensus threshold: >= 25% of datasets
  min_datasets <- ceiling(n_datasets * consensus_threshold)
  consensus_peaks <- union_peaks[n_per_peak >= min_datasets]

  # Per-dataset contribution to the consensus, and bias check
  contribution <- sapply(external_peaks, function(p)
    100 * sum(overlapsAny(consensus_peaks, p)) / length(consensus_peaks))
  max_contribution <- max(contribution)
  if (max_contribution > fail_frac * 100) {
    stop(sprintf("Consensus biased: top dataset contributes %.1f%% (> %.0f%%)",
                 max_contribution, fail_frac * 100))
  } else if (max_contribution > warn_frac * 100) {
    warning(sprintf("Top dataset contributes %.1f%% of consensus (> %.0f%%) - monitor",
                    max_contribution, warn_frac * 100))
  }

  list(consensus_peaks = consensus_peaks,
       support_matrix = support_matrix,
       dataset_contribution = sort(contribution, decreasing = TRUE),
       max_contribution = max_contribution,
       min_datasets = min_datasets)
}
