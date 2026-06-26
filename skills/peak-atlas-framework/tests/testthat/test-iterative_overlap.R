# Invariants for the iterative-overlap merge (scripts/iterative_overlap.R).
source_script("iterative_overlap.R")

test_that("convergeClusterGRanges returns a strictly non-overlapping set", {
  skip_scaffold("verify GenomicRanges is installed and the script sources cleanly")
  gr  <- make_overlapping_clusters(k = 6, per = 5)
  out <- convergeClusterGRanges(gr, by = "score")
  expect_length(findOverlaps(out, drop.self = TRUE, drop.redundant = TRUE), 0L)
})

test_that("merge is summit-faithful: every output peak is an INPUT peak (no midpoints)", {
  skip_scaffold("the key property that distinguishes this from reduce()/bedtools merge")
  gr  <- make_overlapping_clusters(k = 6, per = 5)
  out <- convergeClusterGRanges(gr, by = "score")
  expect_true(all(gr_key(out) %in% gr_key(gr)))
})

test_that("one winner per overlap cluster, and it is the max-score peak", {
  skip_scaffold("k overlap clusters -> exactly k survivors, each its cluster's argmax score")
  gr  <- make_overlapping_clusters(k = 6, per = 5)
  out <- convergeClusterGRanges(gr, by = "score")
  expect_equal(length(out), 6L)
  winners <- tapply(seq_along(gr), mcols(gr)$true_cluster, function(idx) {
    gr_key(gr[idx[which.max(mcols(gr)$score[idx])]])
  })
  expect_setequal(gr_key(out), as.character(winners))
})

test_that("output count never exceeds input; disjoint input is returned intact", {
  skip_scaffold()
  disj <- make_disjoint_peaks(n = 12)
  expect_setequal(gr_key(convergeClusterGRanges(disj, by = "score")), gr_key(disj))
  cl <- make_overlapping_clusters()
  expect_lte(length(convergeClusterGRanges(cl, by = "score")), length(cl))
})

test_that("idempotent and deterministic", {
  skip_scaffold()
  gr    <- make_overlapping_clusters()
  once  <- convergeClusterGRanges(gr, by = "score")
  twice <- convergeClusterGRanges(once, by = "score")
  expect_setequal(gr_key(once), gr_key(twice))                                   # idempotent
  expect_setequal(gr_key(once), gr_key(convergeClusterGRanges(gr, by = "score")))# deterministic
})

test_that("respects an alternate ranking column (e.g. adjusted_score)", {
  skip_scaffold("the embedded pipeline ranks by the support-weighted adjusted_score")
  gr <- make_overlapping_clusters()
  mcols(gr)$adjusted_score <- mcols(gr)$score * 2
  out <- convergeClusterGRanges(gr, by = "adjusted_score")
  expect_length(findOverlaps(out, drop.self = TRUE, drop.redundant = TRUE), 0L)
})

test_that("edge cases: empty and singleton inputs", {
  skip_scaffold()
  empty <- make_disjoint_peaks(n = 2)[0]
  expect_length(convergeClusterGRanges(empty, by = "score"), 0L)
  one <- make_disjoint_peaks(n = 1)
  expect_setequal(gr_key(convergeClusterGRanges(one, by = "score")), gr_key(one))
})
