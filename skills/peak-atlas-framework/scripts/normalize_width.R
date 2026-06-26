# normalize_width.R - Normalize MACS3 peaks to fixed-width summit-centered peaks
#
# Provenance:
#   Reproduced faithfully from normalize_to_501bp() in:
#     /data2/users/JCRLab/JBader/JBader_scHFD/01_scripts/R/peak_utils.R
#
# Why 501bp summit-centered:
#   Fixed-width peaks are the chromVAR / TFBSTools standard - equal-width
#   windows make per-peak signal comparable and motif scanning unbiased by
#   width. 501bp = summit +/- 250 + the 1bp summit. The summit (MACS3 `peak`
#   column = offset from peak start) is the Tn5 insertion maximum, so centering
#   on it keeps the highest-signal base at the window center.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
})

#' Normalize peaks to (2*extend + 1) bp centered on the summit
#'
#' Takes MACS3 narrowPeak-derived GRanges and creates fixed-width peaks
#' centered on the summit position. With extend=250 this yields 501bp peaks.
#'
#' @param peaks GRanges object with MACS3 peaks. Summit is read from the first
#'   available of: `peak` (MACS3 offset-from-start), `summit`, `summit_pos`.
#'   Falls back to the peak midpoint with a warning if none is present.
#' @param extend Integer. bp to extend on each side of the summit (default 250).
#' @return GRanges with width = 2*extend + 1, centered on the summit.
normalize_to_501bp <- function(peaks, extend = 250) {
  stopifnot(extend > 0)

  # Determine summit position
  if ("peak" %in% colnames(mcols(peaks))) {
    # MACS3 'peak' column contains summit offset from start
    centers <- start(peaks) + mcols(peaks)$peak
  } else if ("summit" %in% colnames(mcols(peaks))) {
    centers <- mcols(peaks)$summit
  } else if ("summit_pos" %in% colnames(mcols(peaks))) {
    centers <- mcols(peaks)$summit_pos
  } else {
    warning("No summit column found, using peak midpoint")
    centers <- start(peaks) + floor(width(peaks) / 2)
  }

  # Create 1bp summit GRanges
  summits_1bp <- GRanges(
    seqnames = seqnames(peaks),
    ranges = IRanges(start = centers, width = 1),
    strand = "*"
  )

  # Copy metadata, record the summit position used
  mcols(summits_1bp) <- mcols(peaks)
  mcols(summits_1bp)$summit_pos <- centers

  # Extend to fixed width: extend on each side + 1bp summit = 2*extend + 1
  peaks_fixed <- resize(summits_1bp, width = 2 * extend + 1, fix = "center")

  return(peaks_fixed)
}
