# helper-fixtures.R - fixtures + sourcing for peak-atlas-multiome test scaffolds.
# See peak-atlas-framework/tests/README.md for the testing philosophy (invariants
# / properties on synthetic inputs; skip_scaffold() until wired in an R env).

suppressWarnings(suppressMessages({
  have <- requireNamespace("GenomicRanges", quietly = TRUE)
  if (have) { library(testthat); library(GenomicRanges); library(IRanges) } else { library(testthat) }
  if (requireNamespace("Matrix", quietly = TRUE)) library(Matrix)
}))

skip_scaffold <- function(what = "wire the fixture to the real function and verify in an R env") {
  testthat::skip(paste("SCAFFOLD -", what))
}

skill_root    <- function() normalizePath(file.path(testthat::test_path(), "..", ".."))
source_script <- function(name) {
  try(sys.source(file.path(skill_root(), "scripts", name), envir = globalenv()), silent = TRUE)
}
check_path <- function(name) file.path(skill_root(), "checks", name)
run_check  <- function(name, args = character()) {
  system2("Rscript", c(check_path(name), args), stdout = FALSE, stderr = FALSE)
}
gr_key <- function(gr) paste(as.character(seqnames(gr)), start(gr), end(gr))

# call_peaks_multistrategy.R sources the framework primitives; point it at them.
Sys.setenv(PEAK_ATLAS_FRAMEWORK_SCRIPTS =
  suppressWarnings(normalizePath(
    file.path(skill_root(), "..", "peak-atlas-framework", "scripts"), mustWork = FALSE)))

## --- synthetic fixtures ------------------------------------------------------

# A peaks(3) x cells accessibility matrix with a known rare cluster:
#  peak 1 "rare_marker" accessible in `rare_hits` of the rare cluster's cells only;
#  peak 2 "ubiquitous"  accessible broadly; peak 3 "noise" in 2 cells.
make_counts_fixture <- function(n_rare = 20, n_abund = 1000, rare_hits = 15, seed = 1) {
  set.seed(seed)
  n_cells <- n_rare + n_abund
  labels  <- c(rep("Rare", n_rare), rep("Abundant", n_abund))
  m <- Matrix::Matrix(0, nrow = 3, ncol = n_cells, sparse = TRUE)
  m[1, sample(seq_len(n_rare), rare_hits)] <- 1
  m[2, sample(seq_len(n_cells), floor(0.5 * n_cells))] <- 1
  m[3, sample(seq_len(n_cells), 2)] <- 1
  list(counts = m, labels = labels)
}

# A 3-peak GRanges carrying the n_strategies column the filter keys off.
make_peaks_n_strategies <- function(n_strategies = c(1L, 4L, 1L)) {
  gr <- GRanges("chr1", IRanges(start = c(1, 2000, 4000), width = 501))
  mcols(gr)$n_strategies <- n_strategies
  gr
}

# A tiny MACS3 narrowPeak file (0-based) for read_narrowpeak().
write_narrowpeak_tmp <- function() {
  p  <- tempfile(fileext = ".narrowPeak")
  df <- data.frame(chr = "chr1", start = c(1000, 5000), end = c(1500, 5500),
                   name = c("p1", "p2"), score = c(100, 50), strand = ".",
                   sv = c(5, 3), pv = c(9, 8), qv = c(7, 6), peak = c(250, 250))
  write.table(df, p, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  p
}

# TODO(R-env): stub Signac::CallPeaks to return a known GRanges per group, so
# reconcile_strategies / call_strategy_peaks are testable without MACS3.
mock_call_peaks <- function(...) stop("TODO(R-env): stub Signac::CallPeaks with a known return")
