# Invariants for the adaptive Primary+Rescue filter (scripts/apply_primary_rescue_filter.R).
source_script("apply_primary_rescue_filter.R")

test_that("calculate_adaptive_thresholds = clamp(rate*n, floor, cap)", {
  skip_scaffold()
  th <- calculate_adaptive_thresholds(c(HSC = 147, pDC = 180, Neutrophil = 4500),
                                      floor = 15, rate = 0.02, cap = 200)
  expect_equal(unname(round(th)), c(15, 15, 90))     # the docstring example
  expect_equal(unname(calculate_adaptive_thresholds(c(Huge = 100000))), 200)  # cap
})

test_that("PRIMARY keeps a peak that clears its cluster threshold in ANY one cluster", {
  skip_scaffold()
  fx    <- make_counts_fixture(n_rare = 20, n_abund = 1000, rare_hits = 15)
  peaks <- make_peaks_n_strategies(c(1L, 1L, 1L))    # no rescue eligibility
  keep  <- apply_primary_rescue_filter(peaks, fx$counts, fx$labels,
                                       floor = 15, rate = 0.02, cap = 200, verbose = FALSE)
  expect_true(keep[1])    # rare_marker: 15 hits in the rare cluster >= floor(15)
  expect_true(keep[2])    # ubiquitous
  expect_false(keep[3])   # noise: 2 cells, below every cluster threshold
})

test_that("RARE-TYPE PROTECTION: adaptive keeps a peak a global threshold would drop", {
  skip_scaffold("the whole reason Primary+Rescue exists")
  fx    <- make_counts_fixture(n_rare = 20, n_abund = 1000, rare_hits = 15)
  peaks <- make_peaks_n_strategies(c(1L, 1L, 1L))
  keep  <- apply_primary_rescue_filter(peaks, fx$counts, fx$labels, verbose = FALSE)
  # 15 / 1020 cells globally (~1.5%) - a naive global min (say 50 cells) drops it.
  expect_true(keep[1])
})

test_that("RESCUE keeps 4-strategy peaks with >= rescue_global_min cells even if PRIMARY fails", {
  skip_scaffold()
  fx    <- make_counts_fixture()
  peaks <- make_peaks_n_strategies(c(1L, 4L, 1L))    # peak 2 is 4-strategy
  keep  <- apply_primary_rescue_filter(peaks, fx$counts, fx$labels,
                                       floor = 1e6,    # impossible PRIMARY bar
                                       rescue_global_min = 10, verbose = FALSE)
  expect_true(keep[2])    # rescued by full cross-strategy support
  expect_false(keep[1])   # not 4-strategy, fails PRIMARY
})

test_that("final_keep == primary | rescue; output is a logical of length(peaks)", {
  skip_scaffold()
  fx    <- make_counts_fixture()
  peaks <- make_peaks_n_strategies(c(1L, 4L, 1L))
  keep  <- apply_primary_rescue_filter(peaks, fx$counts, fx$labels, verbose = FALSE)
  expect_type(keep, "logical")
  expect_length(keep, length(peaks))
})

test_that("monotone: lowering the floor keeps a superset", {
  skip_scaffold()
  fx     <- make_counts_fixture()
  peaks  <- make_peaks_n_strategies(c(1L, 1L, 1L))
  strict <- apply_primary_rescue_filter(peaks, fx$counts, fx$labels, floor = 50, verbose = FALSE)
  relax  <- apply_primary_rescue_filter(peaks, fx$counts, fx$labels, floor = 5,  verbose = FALSE)
  expect_true(all(which(strict) %in% which(relax)))
})
