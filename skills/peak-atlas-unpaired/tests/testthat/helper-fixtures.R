# helper-fixtures.R - fixtures + sourcing for peak-atlas-unpaired test scaffolds.
# See peak-atlas-framework/tests/README.md for the testing philosophy.

suppressWarnings(suppressMessages({
  have <- requireNamespace("GenomicRanges", quietly = TRUE)
  if (have) {
    library(testthat); library(GenomicRanges); library(IRanges); library(S4Vectors)
  } else {
    library(testthat)
  }
}))

skip_scaffold <- function(what = "wire the fixture to the real function and verify in an R env") {
  testthat::skip(paste("SCAFFOLD -", what))
}

skill_root    <- function() normalizePath(file.path(testthat::test_path(), "..", ".."))
source_script <- function(name) {
  try(sys.source(file.path(skill_root(), "scripts", name), envir = globalenv()), silent = TRUE)
}
check_path    <- function(name) file.path(skill_root(), "checks", name)
run_check     <- function(name, args = character()) {
  system2("Rscript", c(check_path(name), args), stdout = FALSE, stderr = FALSE)
}
gr_key       <- function(gr) paste(as.character(seqnames(gr)), start(gr), end(gr))
write_rds_tmp <- function(x) { p <- tempfile(fileext = ".rds"); saveRDS(x, p); p }

# build_external_consensus with the bias gate disabled, for non-bias tests
# (small synthetic inputs trivially trip the 50% bias guard otherwise).
bc_nobias <- function(external_peaks, consensus_threshold = 0.25, gap_width = 100) {
  build_external_consensus(external_peaks, consensus_threshold = consensus_threshold,
                           gap_width = gap_width, warn_frac = 2, fail_frac = 2)
}

## --- synthetic fixtures ------------------------------------------------------

# A union of A/B/C peaks with the per-peak metadata assign_tiers() expects.
#  row 1: rare + B(score>15, n_groups>0) AND 3-strategy -> Tier 0 (priority test)
#  row 2: 3-strategy, not rare                          -> Tier 1
#  row 3: rare but no B support                          -> falls through (not Tier 0)
#  row 4: single strategy, nothing special              -> NotReproducible
make_tier_union <- function() {
  gr <- GRanges("chr1", IRanges(start = c(1, 2000, 4000, 6000), width = 501))
  mcols(gr) <- S4Vectors::DataFrame(
    strategy_a = c(TRUE,  TRUE,  TRUE,  FALSE),
    strategy_b = c(TRUE,  TRUE,  FALSE, FALSE),
    strategy_c = c(TRUE,  TRUE,  FALSE, TRUE),
    has_rare_celltype        = c(TRUE,  FALSE, TRUE,  FALSE),
    has_celltype_of_interest = c(TRUE,  FALSE, TRUE,  FALSE),
    a_max_n_cells = c(100, 100, 100, 100),
    b_max_n_cells = c(50,  50,  0,   0),
    b_max_score   = c(30,  5,   0,   0),
    c_max_score   = c(5,   5,   0,   2),
    a_n_groups    = c(1,   1,   1,   0),
    b_n_groups    = c(2,   2,   0,   0),
    a_n_clusters  = c(1,   1,   1,   0),
    b_n_clusters  = c(1,   1,   0,   0),
    b_n_conditions = c(2,  1,   0,   0),
    c_n_conditions = c(1,  1,   0,   1)
  )
  gr
}

# Three independent external datasets sharing one region (1000) plus private
# regions, for consensus / bias tests.
make_external_datasets <- function() {
  list(
    d1 = GRanges("chr1", IRanges(c(1000, 5000),        width = 501)),
    d2 = GRanges("chr1", IRanges(c(1000, 9000),        width = 501)),
    d3 = GRanges("chr1", IRanges(c(1000, 5000, 13000), width = 501))
  )
}
