# Contract for per-cell FRiP (scripts/frip.R). calculate_frip needs a Signac
# object with a Fragments()-backed assay, so these stay skipped until the toy
# fixture make_toy_atac_seurat() exists. The assertions ARE the spec.
source_script("frip.R")

test_that("FRiP is a per-cell fraction in [0, 1]", {
  skip_scaffold("needs make_toy_atac_seurat() - a minimal Signac object with fragments")
  obj       <- make_toy_atac_seurat()
  toy_peaks <- make_disjoint_peaks(50)
  frip      <- calculate_frip(obj, peaks = toy_peaks, total_frags_col = "atac_fragments")
  expect_length(frip, ncol(obj))
  expect_true(all(frip >= 0 & frip <= 1))
})

test_that("FRiP is monotone in the peak set: a superset cannot lower it", {
  skip_scaffold("adding peaks can only add reads-in-peaks")
  obj   <- make_toy_atac_seurat()
  big   <- make_disjoint_peaks(50)
  small <- big[1:5]
  expect_true(all(calculate_frip(obj, big) >= calculate_frip(obj, small)))
})

test_that("cells with zero total fragments yield a defined value (no Inf/NaN crash)", {
  skip_scaffold("decide the empty-cell behaviour (drop or 0), then assert it explicitly")
  # obj <- make_toy_atac_seurat(with_empty_cell = TRUE)
  # expect_false(any(is.nan(calculate_frip(obj, make_disjoint_peaks(5)))))
})
