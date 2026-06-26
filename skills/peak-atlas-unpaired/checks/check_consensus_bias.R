#!/usr/bin/env Rscript
# check_consensus_bias.R - Assert no single external dataset dominates the
#                          consensus. Exit non-zero on failure.
#
# Provenance:
#   Bias thresholds from foundation/S4.1b_build_consensus.R
#     ("Single dataset contributes >50% - consensus may be biased!" / >40% note)
#   Reproduces scripts/build_consensus.R::build_external_consensus bias logic.
#
# Gate: the maximum single-dataset contribution to the consensus peak set must
# stay below the fail threshold. A high contribution means one dataset (or an
# un-merged replicate set) is driving the consensus rather than agreement.
#   contribution(d) = 100 * (# consensus peaks overlapping dataset d) / (# consensus peaks)
#   warn if max contribution > warn_pct (default 40)
#   FAIL if max contribution > fail_pct (default 50)
#
# Input options (one of):
#   A) An RDS with a numeric vector of per-dataset contribution percentages
#      (e.g. build_external_consensus()$dataset_contribution).
#   B) A two-column CSV: dataset, pct_contribution.
#
# Usage:
#   Rscript check_consensus_bias.R <contribution.rds|csv> [warn_pct] [fail_pct]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript check_consensus_bias.R <contribution.rds|csv> [warn_pct] [fail_pct]")
}

path     <- args[[1]]
warn_pct <- if (length(args) >= 2) as.numeric(args[[2]]) else 40
fail_pct <- if (length(args) >= 3) as.numeric(args[[3]]) else 50

load_contribution <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    x <- readRDS(path)
    if (is.list(x) && "dataset_contribution" %in% names(x)) x <- x$dataset_contribution
    v <- as.numeric(x)
    names(v) <- names(x)
    v
  } else {
    df <- read.csv(path, header = TRUE, stringsAsFactors = FALSE)
    pct_col <- if ("pct_contribution" %in% names(df)) "pct_contribution" else names(df)[ncol(df)]
    name_col <- if ("dataset" %in% names(df)) "dataset" else names(df)[1]
    v <- as.numeric(df[[pct_col]])
    names(v) <- df[[name_col]]
    v
  }
}

contribution <- load_contribution(path)
contribution <- contribution[is.finite(contribution)]

if (length(contribution) == 0) {
  message("[FAIL] No dataset contributions found - cannot assess bias")
  quit(status = 1)
}

max_contribution <- max(contribution)
top_dataset <- names(contribution)[which.max(contribution)]
if (is.null(top_dataset) || is.na(top_dataset)) top_dataset <- "<unnamed>"

message("=== Consensus Dataset-Bias Check ===")
message(sprintf("Datasets: %d", length(contribution)))
message(sprintf("Top contributor: %s = %.1f%%", top_dataset, max_contribution))
message(sprintf("Thresholds: warn > %.0f%%, fail > %.0f%%", warn_pct, fail_pct))

if (max_contribution > fail_pct) {
  message(sprintf("[FAIL] %s contributes %.1f%% (> %.0f%%) - consensus is biased; merge replicates or add datasets",
                  top_dataset, max_contribution, fail_pct))
  quit(status = 1)
} else if (max_contribution > warn_pct) {
  message(sprintf("[WARN] %s contributes %.1f%% (> %.0f%%) - acceptable but monitor",
                  top_dataset, max_contribution, warn_pct))
  quit(status = 0)
} else {
  message("[PASS] No single dataset dominates the consensus")
  quit(status = 0)
}
