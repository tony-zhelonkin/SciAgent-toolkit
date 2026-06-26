# helper-fixtures.R - fixtures + sourcing for peak-atlas-framework test scaffolds.
#
# Testing philosophy (see ../README.md): assert INVARIANTS / CONTRACTS /
# PROPERTIES on SYNTHETIC inputs, never byte-for-byte golden outputs on frozen
# data. The GRanges builders below are concrete and deterministic; heavy
# fixtures (Seurat + fragments) are TODO stubs. Every test opens with
# skip_scaffold() so the suite stays green until implemented + verified in an
# R + Bioconductor env.

suppressWarnings(suppressMessages({
  has_gr <- requireNamespace("GenomicRanges", quietly = TRUE)
  if (has_gr) {
    library(testthat); library(GenomicRanges); library(IRanges)
  } else {
    library(testthat)
  }
}))

skip_scaffold <- function(what = "wire the fixture to the real function and verify in an R env") {
  testthat::skip(paste("SCAFFOLD -", what))
}

# Resolve <skill>/{scripts,checks}/ from tests/testthat/.
skill_root    <- function() normalizePath(file.path(testthat::test_path(), "..", ".."))
# Source a script into the global env; tolerant of missing Bioconductor deps
# (tests skip_scaffold() first, so undefined functions are never reached).
source_script <- function(name) {
  try(sys.source(file.path(skill_root(), "scripts", name), envir = globalenv()),
      silent = TRUE)
}
check_path <- function(name) file.path(skill_root(), "checks", name)

# Run a CLI check script and return its exit status.
run_check <- function(name, args = character()) {
  system2("Rscript", c(check_path(name), args), stdout = FALSE, stderr = FALSE)
}

# Stable "is this exact interval present" key (used for summit-fidelity checks).
gr_key <- function(gr) paste(as.character(seqnames(gr)), start(gr), end(gr))

## --- synthetic peak fixtures (concrete, non-brittle) ------------------------

# n disjoint fixed-width peaks with a unique integer `score`. No two overlap.
make_disjoint_peaks <- function(n = 20, width = 501, gap = 1000, seed = 1) {
  set.seed(seed)
  gr <- GRanges("chr1", IRanges(start = seq(1, by = gap, length.out = n), width = width))
  mcols(gr)$score <- sample(seq_len(n))
  gr
}

# k overlap clusters of `per` mutually-overlapping peaks; a correct merge keeps
# exactly k peaks (the max-score peak of each cluster). true_cluster labels them.
make_overlapping_clusters <- function(k = 5, per = 4, width = 501, seed = 1) {
  set.seed(seed)
  parts <- lapply(seq_len(k), function(i) {
    base <- (i - 1L) * 5000L + 1L
    GRanges("chr1", IRanges(start = base + sample(0:150, per), width = width))
  })
  gr <- do.call(c, parts)
  mcols(gr)$score <- sample(length(gr))            # unique scores -> unambiguous winner
  mcols(gr)$true_cluster <- rep(seq_len(k), each = per)
  gr
}

# Peaks whose true summit is known, given as a MACS-style `peak` offset-from-start.
make_summit_peaks <- function(summits = c(1000, 8000, 20000), in_width = 137) {
  starts <- summits - (in_width %/% 2L)
  GRanges("chr1", IRanges(start = starts, width = in_width),
          peak = summits - starts,                 # MACS3 offset-from-start
          score = seq_along(summits))
}

# TODO(R-env): a minimal Signac ATAC object with a Fragments()-backed assay and
# an `atac_fragments` metadata column, for calculate_frip(). Keep it SMALL and
# SYNTHETIC (a hand-written fragments.tsv.gz over a few cells + peaks).
make_toy_atac_seurat <- function(...) {
  stop("TODO(R-env): build a minimal Signac object with a fragments file + peaks")
}
