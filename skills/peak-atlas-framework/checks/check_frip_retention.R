#!/usr/bin/env Rscript
# check_frip_retention.R - Assert FRiP retention of a filtered atlas, exit non-zero on failure
#
# Provenance:
#   Threshold from validation.frip_retention_min in:
#     /data2/users/JCRLab/JBader/JBader_scHFD/02_analysis/config/peak_filtration_config.yaml
#   FRiP computed via calculate_frip() (scripts/frip.R), reproduced from
#     /data2/users/JCRLab/JBader/JBader_scHFD/01_scripts/R/peak_utils.R
#
# Gate: median FRiP of the FILTERED atlas must be >= threshold (default 0.90)
# of the median FRiP of the ORIGINAL atlas. Losing more than 10% of FRiP means
# the prune dropped too much real signal.
#
# Inputs are precomputed per-cell FRiP vectors (recommended) so the check stays
# fast and dependency-light. Each vector may be a plain RDS numeric vector or a
# one-column CSV. (Compute them with scripts/frip.R::calculate_frip.)
#
# Usage:
#   Rscript check_frip_retention.R original_frip.rds filtered_frip.rds [min_ratio]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript check_frip_retention.R <original_frip> <filtered_frip> [min_ratio]")
}

orig_path <- args[[1]]
filt_path <- args[[2]]
min_ratio <- if (length(args) >= 3) as.numeric(args[[3]]) else 0.90

load_frip <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    v <- readRDS(path)
  } else {
    df <- read.csv(path, header = TRUE)
    v <- df[[ncol(df)]]
  }
  v <- as.numeric(v)
  v[is.finite(v)]
}

orig <- load_frip(orig_path)
filt <- load_frip(filt_path)

if (length(orig) == 0 || length(filt) == 0) {
  message("[FAIL] Empty FRiP vector(s) - cannot assess retention")
  quit(status = 1)
}

med_orig <- median(orig)
med_filt <- median(filt)
ratio <- med_filt / med_orig

message("=== FRiP Retention ===")
message(sprintf("Original median FRiP: %.4f", med_orig))
message(sprintf("Filtered median FRiP: %.4f", med_filt))
message(sprintf("Retention ratio: %.4f (threshold >= %.2f)", ratio, min_ratio))

if (ratio >= min_ratio) {
  message("[PASS] FRiP retention OK")
  quit(status = 0)
} else {
  message(sprintf("[FAIL] FRiP retention %.4f below threshold %.2f", ratio, min_ratio))
  quit(status = 1)
}
