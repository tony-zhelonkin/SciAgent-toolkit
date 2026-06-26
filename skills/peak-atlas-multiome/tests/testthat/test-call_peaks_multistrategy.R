# Invariants for multi-strategy calling (scripts/call_peaks_multistrategy.R).
# The pure helpers are tested directly; MACS3/Signac paths are mocked.
source_script("call_peaks_multistrategy.R")

test_that("split_pseudoreplicates partitions cells ~50/50, deterministically", {
  skip_scaffold()
  cells <- paste0("c", 1:101)
  s <- split_pseudoreplicates(cells, seed = 42)
  expect_setequal(c(s$rep1, s$rep2), cells)          # union == input
  expect_length(intersect(s$rep1, s$rep2), 0L)       # disjoint
  expect_equal(length(s$rep1), floor(101 / 2))
  expect_identical(split_pseudoreplicates(cells, seed = 42), s)  # deterministic
})

test_that("filter_by_pseudoreps keeps only rep1 peaks overlapping rep2", {
  skip_scaffold()
  rep1 <- GRanges("chr1", IRanges(c(1000, 5000, 9000), width = 501))  # 9000 = rep1-only
  rep2 <- GRanges("chr1", IRanges(c(1000, 5000),       width = 501))
  out  <- filter_by_pseudoreps(rep1, rep2, verbose = FALSE)
  expect_true(all(gr_key(out) %in% gr_key(rep1)))    # subset of rep1
  expect_true(all(overlapsAny(out, rep2)))           # all reproducible
  expect_length(out, 2L)                             # the rep1-only peak dropped
  expect_length(filter_by_pseudoreps(GRanges(), rep2, verbose = FALSE), 0L)  # empty-safe
})

test_that("read_narrowpeak converts 0-based starts to 1-based and keeps the summit offset", {
  skip_scaffold("needs a narrowPeak file (write_narrowpeak_tmp)")
  gr <- read_narrowpeak(write_narrowpeak_tmp())
  expect_equal(start(gr), c(1001, 5001))             # 0-based 1000/5000 -> 1-based
  expect_true("peak" %in% colnames(mcols(gr)))
})

test_that("reconcile_strategies: non-overlapping output with support metadata", {
  skip_scaffold("needs framework primitives sourced + a BSgenome for the boundary filter")
  mk <- function(starts, sc) {
    g <- GRanges("chr1", IRanges(starts, width = 501)); mcols(g)$score <- sc; list(g1 = g)
  }
  strat <- list(RNA = mk(c(1000, 5000), c(10, 8)), ATAC = mk(1000, 9),
                WNN = mk(1000, 7),                 CellType = mk(5000, 6))
  # out <- reconcile_strategies(strat, blacklist = GRanges(), BSgenome = toy_bsgenome)
  # expect_length(findOverlaps(out, drop.self = TRUE, drop.redundant = TRUE), 0L)
  # expect_true("n_strategies" %in% colnames(mcols(out)))
  # expect_true(all(grepl("^chr1:", names(out))))
})

test_that("CONTRACT: per-group calling requires group.by (mock MACS3)", {
  skip_scaffold("mock_call_peaks - without group.by every group returns identical peaks (a real Signac footgun)")
  # With a CallPeaks mock, assert call_strategy_peaks passes group.by/idents per
  # group, skips groups below min_cells, and scores merged peaks by #supporting groups.
})
