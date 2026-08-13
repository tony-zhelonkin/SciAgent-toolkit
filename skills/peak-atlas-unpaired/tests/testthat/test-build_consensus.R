# Invariants for external consensus (scripts/build_consensus.R). merge_replicates
# and build_external_consensus are pure GRanges logic - fully testable here.
source_script("build_consensus.R")

test_that("merge_replicates keeps peaks present in >= 50% of replicates (ceiling)", {
  skip_scaffold()
  peaks_list <- list(
    r1 = GRanges("chr1", IRanges(c(1000, 5000), width = 501)),
    r2 = GRanges("chr1", IRanges(c(1000, 9000), width = 501)),
    r3 = GRanges("chr1", IRanges(c(1000),       width = 501))
  )
  merged <- merge_replicates(peaks_list, list(grp = c("r1", "r2", "r3")))$merged_peaks$grp
  # 1000 is in all 3 (kept); 5000 and 9000 are each in 1/3 < ceil(1.5)=2 (dropped)
  expect_true(overlapsAny(GRanges("chr1", IRanges(1000, width = 501)), merged))
  expect_length(merged, 1L)
})

test_that("build_external_consensus keeps peaks in >= 25% of datasets", {
  skip_scaffold()
  res <- bc_nobias(make_external_datasets(), consensus_threshold = 0.25)
  expect_equal(res$min_datasets, ceiling(3 * 0.25))               # = 1
  expect_true(all(mcols(res$consensus_peaks)$n_datasets >= res$min_datasets))
})

test_that("TWO-GAP-WIDTH rule: 70bp-apart peaks merge at gap=100 but NOT at gap=50", {
  skip_scaffold("encodes the distinction in references/external-consensus.md")
  near <- list(
    a = GRanges("chr1", IRanges(1000, width = 501)),                 # 1000-1500
    b = GRanges("chr1", IRanges(1000 + 501 + 70, width = 501))       # 70bp gap after a
  )
  wide   <- bc_nobias(near, consensus_threshold = 0, gap_width = 100)
  narrow <- bc_nobias(near, consensus_threshold = 0, gap_width = 50)
  expect_lt(length(wide$consensus_peaks), length(narrow$consensus_peaks))  # 1 vs 2
})

test_that("bias gate errors when a single dataset contributes > 50% of the consensus", {
  skip_scaffold()
  biased <- list(
    big   = GRanges("chr1", IRanges(seq(1, by = 1000, length.out = 20), width = 501)),
    small = GRanges("chr1", IRanges(1, width = 501))
  )
  expect_error(build_external_consensus(biased, consensus_threshold = 0, fail_frac = 0.50))
})

test_that("monotone: a higher consensus threshold yields no more peaks", {
  skip_scaffold()
  d <- make_external_datasets()
  expect_gte(length(bc_nobias(d, 0.25)$consensus_peaks),
             length(bc_nobias(d, 0.75)$consensus_peaks))
})

test_that("empty datasets are dropped; all-empty errors", {
  skip_scaffold()
  expect_error(build_external_consensus(list(a = GRanges(), b = GRanges())))
})
