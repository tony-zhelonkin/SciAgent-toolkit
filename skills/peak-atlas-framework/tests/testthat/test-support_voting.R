# Invariants for cross-strategy support + score boost (scripts/support_voting.R).
source_script("support_voting.R")

three_strategies <- function() {
  list(
    RNA  = GRanges("chr1", IRanges(c(1000, 5000), width = 501)),
    ATAC = GRanges("chr1", IRanges(c(1000, 9000), width = 501)),
    WNN  = GRanges("chr1", IRanges(c(1000),       width = 501))
  )
}

test_that("n_strategies counts overlapping strategies and is bounded by their number", {
  skip_scaffold()
  s     <- three_strategies()
  peaks <- GRanges("chr1", IRanges(c(1000, 5000, 9000, 50000), width = 501))
  out   <- calculate_strategy_support(peaks, s)
  expect_equal(mcols(out)$n_strategies, c(3L, 1L, 1L, 0L))   # 1000 in all 3; 50000 in none
  expect_true(all(mcols(out)$n_strategies <= length(s)))
})

test_that("strategies string names exactly the overlapping strategies", {
  skip_scaffold()
  s   <- three_strategies()
  out <- calculate_strategy_support(GRanges("chr1", IRanges(1000, width = 501)), s)
  expect_setequal(strsplit(mcols(out)$strategies, ",")[[1]], names(s))
})

test_that("adjusted_score = score*(1 + boost*(n_strategies-1)) and is monotone in support", {
  skip_scaffold()
  peaks <- GRanges("chr1", IRanges(c(1, 1000, 2000, 3000), width = 501), score = rep(10, 4))
  mcols(peaks)$n_strategies <- c(1L, 2L, 3L, 4L)
  out <- add_adjusted_score(peaks, boost = 0.5)
  expect_equal(mcols(out)$adjusted_score, c(10, 15, 20, 25))
  expect_true(all(diff(mcols(out)$adjusted_score) > 0))
})

test_that("REGRESSION: support reflects the per-strategy sets, not a bridged merge", {
  skip_scaffold("the single most important correctness property (compute support BEFORE the final merge)")
  # Two strategies whose peaks do NOT overlap each other. Support must mark each
  # peak as single-strategy, NOT bridge them to n=2.
  s     <- list(A = GRanges("chr1", IRanges(1000, width = 501)),
                B = GRanges("chr1", IRanges(9000, width = 501)))
  peaks <- GRanges("chr1", IRanges(c(1000, 9000), width = 501))
  out   <- calculate_strategy_support(peaks, s)
  expect_equal(mcols(out)$n_strategies, c(1L, 1L))   # NOT c(2L, 2L)
})
