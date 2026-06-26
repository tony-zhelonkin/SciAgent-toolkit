# frip.R - Fraction of Reads in Peaks (FRiP) per cell
#
# Provenance:
#   Reproduced faithfully from calculate_frip() in:
#     /data2/users/JCRLab/JBader/JBader_scHFD/01_scripts/R/peak_utils.R
#
# FRiP = (fragments overlapping the peak set) / (total fragments) per cell.
# It is the headline retention metric for the FILTER and VALIDATE stages:
# a filtered atlas should keep >=0.90 of the original per-cell median FRiP
# (see checks/check_frip_retention.R and references/validation-battery.md).

suppressPackageStartupMessages({
  library(Matrix)
})

#' Calculate Fraction of Reads in Peaks (FRiP) per cell
#'
#' Counts fragments overlapping the peak set via Signac::FeatureMatrix and
#' divides by the per-cell total fragments.
#'
#' @param obj Seurat object with an ATAC assay carrying Fragments().
#' @param peaks GRanges object with the peak set to score.
#' @param assay Character. Assay holding fragment data (default "ATAC").
#' @param total_frags_col Character. Per-cell total-fragments metadata column
#'   (default "atac_fragments").
#' @return Numeric vector of FRiP per cell (named by barcode).
calculate_frip <- function(obj, peaks, assay = "ATAC",
                           total_frags_col = "atac_fragments") {
  # Fragment file objects from the assay
  assay_obj <- obj[[assay]]
  frags <- Signac::Fragments(assay_obj)

  # Count fragments overlapping peaks: returns sparse peaks x cells matrix
  peak_counts <- Signac::FeatureMatrix(
    fragments = frags,
    features = peaks,
    cells = colnames(obj)
  )

  # Sum counts across all peaks for each cell
  counts_in_peaks <- Matrix::colSums(peak_counts)

  # Per-cell total fragments from metadata
  total_frags <- obj[[total_frags_col, drop = TRUE]]

  # FRiP per cell
  frip_per_cell <- counts_in_peaks / total_frags

  return(frip_per_cell)
}
