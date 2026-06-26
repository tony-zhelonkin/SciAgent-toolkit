# Invariants for blacklist handling (scripts/blacklist.R).
source_script("blacklist.R")

test_that("remove_blacklist_peaks against an empty blacklist is a no-op", {
  skip_scaffold()
  gr <- make_disjoint_peaks(10)
  expect_identical(remove_blacklist_peaks(gr, GRanges()), gr)
})

test_that("output has zero blacklist overlaps and is a subset of the input", {
  skip_scaffold()
  gr  <- make_disjoint_peaks(10, gap = 1000)            # peaks at 1, 1001, 2001, ...
  bl  <- GRanges("chr1", IRanges(start = c(1001, 4001), width = 100))
  out <- remove_blacklist_peaks(gr, bl)
  expect_equal(sum(overlapsAny(out, bl)), 0L)
  expect_true(all(gr_key(out) %in% gr_key(gr)))
  expect_equal(length(out), length(gr) - 2L)
})

test_that("load_blacklist degrades to an empty GRanges for NULL / missing paths", {
  skip_scaffold("a missing blacklist must warn and continue, never error")
  expect_warning(bl1 <- load_blacklist(bed_path = NULL, verbose = FALSE))
  expect_length(bl1, 0L)
  expect_warning(bl2 <- load_blacklist(bed_path = tempfile(fileext = ".bed"), verbose = FALSE))
  expect_length(bl2, 0L)
})
