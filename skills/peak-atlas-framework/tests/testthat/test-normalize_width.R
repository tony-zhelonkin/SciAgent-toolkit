# Invariants for fixed-width summit-centering (scripts/normalize_width.R).
source_script("normalize_width.R")

test_that("all peaks become 2*extend+1 wide", {
  skip_scaffold()
  p <- make_summit_peaks()
  expect_true(all(width(normalize_to_501bp(p, extend = 250)) == 501L))
  expect_true(all(width(normalize_to_501bp(p, extend = 150)) == 301L))
})

test_that("peaks are centered on the summit", {
  skip_scaffold("summit fidelity: window center == MACS 'peak' summit position")
  p   <- make_summit_peaks(summits = c(1000, 8000, 20000))
  out <- normalize_to_501bp(p, extend = 250)
  centers <- start(out) + (width(out) - 1L) %/% 2L
  expect_equal(centers, c(1000, 8000, 20000))
  expect_equal(mcols(out)$summit_pos, c(1000, 8000, 20000))
})

test_that("idempotent in width and center when driven by summit_pos", {
  skip_scaffold("CONTRACT: normalize ONCE. (Re-running with the original MACS 'peak' offset column still present mis-centers, since 'peak' is an offset from the ORIGINAL start. Callers normalize once; here we drive the second pass from summit_pos.)")
  p <- GRanges("chr1", IRanges(start = c(1000, 8000), width = 1), summit_pos = c(1000, 8000))
  once  <- normalize_to_501bp(p)
  twice <- normalize_to_501bp(once)
  expect_true(all(width(twice) == 501L))
  expect_equal(start(twice), start(once))
})

test_that("summit column precedence: peak > summit > summit_pos > midpoint", {
  skip_scaffold("documents the documented fallback order")
  p <- GRanges("chr1", IRanges(start = 901, width = 200), peak = 100, summit = 5)
  expect_equal(mcols(normalize_to_501bp(p))$summit_pos, 1001)   # 'peak' wins: 901 + 100
})

test_that("non-summit metadata is preserved", {
  skip_scaffold()
  p <- make_summit_peaks()
  mcols(p)$tag <- letters[seq_along(p)]
  expect_equal(mcols(normalize_to_501bp(p))$tag, letters[seq_along(p)])
})
