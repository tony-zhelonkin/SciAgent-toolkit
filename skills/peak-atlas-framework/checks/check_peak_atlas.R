#!/usr/bin/env Rscript
# check_peak_atlas.R - Validate a fixed-width peak atlas, exit non-zero on failure
#
# Provenance:
#   Reproduces validate_peak_atlas() from:
#     /data2/users/JCRLab/JBader/JBader_scHFD/01_scripts/R/peak_utils.R
#
# Assertions (any failure => exit 1):
#   - No overlapping peaks in the atlas
#   - Every peak exactly the fixed width (default 501bp)
#   - Zero blacklist overlaps (when a blacklist is provided)
#   - Valid coordinates (start >= 1, start <= end)
#
# Usage:
#   Rscript check_peak_atlas.R atlas.rds [blacklist.bed] [expected_width]
#   # atlas.rds may be an .rds (GRanges) or a BED file.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript check_peak_atlas.R <atlas.rds|atlas.bed> [blacklist.bed] [expected_width]")
}

atlas_path <- args[[1]]
blacklist_path <- if (length(args) >= 2 && nzchar(args[[2]])) args[[2]] else NULL
expected_width <- if (length(args) >= 3) as.integer(args[[3]]) else 501L

load_peaks <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    obj <- readRDS(path)
    if (inherits(obj, "GRanges")) return(obj)
    if (inherits(obj, "RangedSummarizedExperiment")) return(rowRanges(obj))
    stop("RDS does not contain a GRanges / RangedSummarizedExperiment")
  }
  rtracklayer::import(path)
}

peaks <- load_peaks(atlas_path)
blacklist <- if (!is.null(blacklist_path)) rtracklayer::import.bed(blacklist_path) else NULL

issues <- character()

# Check 1: overlapping peaks
overlaps <- findOverlaps(peaks, drop.self = TRUE, drop.redundant = TRUE)
n_overlaps <- length(overlaps)
if (n_overlaps > 0) {
  issues <- c(issues, sprintf("%d overlapping peak pairs", n_overlaps))
}

# Check 2: width uniformity
widths <- width(peaks)
if (!all(widths == expected_width)) {
  n_wrong <- sum(widths != expected_width)
  issues <- c(issues, sprintf("%d peaks not %dbp (range %d-%d)",
                              n_wrong, expected_width, min(widths), max(widths)))
}

# Check 3: blacklist contamination
if (!is.null(blacklist) && length(blacklist) > 0) {
  n_bl <- sum(overlapsAny(peaks, blacklist))
  if (n_bl > 0) issues <- c(issues, sprintf("%d peaks overlap blacklist", n_bl))
}

# Check 4: valid coordinates
bad_coords <- sum(start(peaks) < 1 | start(peaks) > end(peaks))
if (bad_coords > 0) {
  issues <- c(issues, sprintf("%d peaks with invalid coordinates", bad_coords))
}

message("=== Peak Atlas Validation ===")
message(sprintf("Total peaks: %s", format(length(peaks), big.mark = ",")))

if (length(issues) == 0) {
  message("[PASS] No overlaps; all peaks ", expected_width, "bp; coordinates valid",
          if (!is.null(blacklist)) "; no blacklist contamination" else "")
  message("\nValidation: PASSED")
  quit(status = 0)
} else {
  message("[FAIL] ", paste(issues, collapse = "; "))
  message("\nValidation: FAILED")
  quit(status = 1)
}
