# blacklist.R - Load an ENCODE blacklist and remove overlapping peaks
#
# Provenance:
#   Reproduced faithfully from load_blacklist_mm39() and
#   remove_blacklist_peaks() in:
#     /data2/users/JCRLab/JBader/JBader_scHFD/01_scripts/R/peak_utils.R
#
#   The source project works in mm39 and lifts the ENCODE mm10 blacklist v2
#   over to mm39. That project-specific genome/path is PARAMETERIZED here:
#   pass your own `bed_path` (and `cache_path`) for any genome build.
#
# Blacklists are ENCODE/Boyle-Lab anomalous high-signal regions. Removing peaks
# that overlap them (subsetByOverlaps invert=TRUE) drops artifactual signal.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
})

#' Load an ENCODE blacklist BED into GRanges (with optional RDS cache)
#'
#' @param bed_path Path to an uncompressed blacklist BED. Build-agnostic: supply
#'   e.g. an mm39 liftover of mm10-blacklist.v2, hg38-blacklist.v2.bed, etc.
#'   If NULL or missing, returns an empty GRanges with a warning (peaks may then
#'   overlap problematic regions).
#' @param cache_path Optional path to cache the parsed GRanges as RDS to speed
#'   up repeated loads. If it exists it is read directly.
#' @param verbose Logical. Print progress (default TRUE).
#' @return GRanges of blacklist regions (possibly empty).
load_blacklist <- function(bed_path = NULL, cache_path = NULL, verbose = TRUE) {
  # Use cache if present
  if (!is.null(cache_path) && file.exists(cache_path)) {
    if (verbose) message("Loading cached blacklist from ", cache_path)
    return(readRDS(cache_path))
  }

  if (is.null(bed_path) || !file.exists(bed_path)) {
    warning("Blacklist BED not found at: ", if (is.null(bed_path)) "NULL" else bed_path)
    warning("Using empty GRanges - peaks may overlap problematic regions")
    return(GRanges())
  }

  if (verbose) message("Loading blacklist from ", bed_path)

  tryCatch({
    bl <- rtracklayer::import.bed(bed_path)

    if (!is.null(cache_path)) {
      dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
      saveRDS(bl, cache_path)
    }

    if (verbose) message(sprintf("Loaded %d blacklist regions", length(bl)))
    return(bl)
  }, error = function(e) {
    warning("Could not load blacklist: ", e$message)
    return(GRanges())
  })
}

#' Remove peaks overlapping blacklist regions
#'
#' @param gr GRanges with peaks.
#' @param blacklist GRanges with blacklist regions (empty = no-op).
#' @return GRanges with blacklist-overlapping peaks removed.
remove_blacklist_peaks <- function(gr, blacklist) {
  if (length(blacklist) == 0) {
    return(gr)
  }
  return(subsetByOverlaps(gr, blacklist, invert = TRUE))
}
