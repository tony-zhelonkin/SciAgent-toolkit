# Exit-code contract for the CLI gates (checks/check_peak_atlas.R,
# checks/check_frip_retention.R). Each runs via Rscript on a constructed fixture
# and must return 0 on a clean input and 1 on each injected violation.
source_script("normalize_width.R")  # not required, but keeps the helper env consistent

write_rds_tmp <- function(x) { p <- tempfile(fileext = ".rds"); saveRDS(x, p); p }

test_that("check_peak_atlas passes a valid atlas and fails each violation", {
  skip_scaffold("needs Rscript + GenomicRanges/rtracklayer; runs the CLI and checks exit codes")
  good <- make_disjoint_peaks(n = 20, width = 501, gap = 1000)
  expect_equal(run_check("check_peak_atlas.R", write_rds_tmp(good)), 0L)

  overlap <- c(good, GRanges("chr1", IRanges(start(good)[1] + 10, width = 501)))
  expect_equal(run_check("check_peak_atlas.R", write_rds_tmp(overlap)), 1L)   # overlapping peaks

  badw <- good
  badw[1] <- resize(badw[1], 777, fix = "center")
  expect_equal(run_check("check_peak_atlas.R", write_rds_tmp(badw)), 1L)      # wrong width
})

test_that("check_peak_atlas flags blacklist contamination when a blacklist is given", {
  skip_scaffold("needs Rscript; second arg is a blacklist BED")
  good <- make_disjoint_peaks(n = 10, width = 501, gap = 1000)
  bl   <- GRanges("chr1", IRanges(start(good)[2], width = 200))
  bl_path <- tempfile(fileext = ".bed")
  rtracklayer::export.bed(bl, bl_path)
  expect_equal(run_check("check_peak_atlas.R", c(write_rds_tmp(good), bl_path)), 1L)
})

test_that("check_frip_retention passes >= threshold and fails below", {
  skip_scaffold("needs Rscript; FRiP vectors as RDS")
  orig <- rep(0.40, 100)
  ok   <- rep(0.38, 100)   # ratio 0.95 >= 0.90 default
  bad  <- rep(0.30, 100)   # ratio 0.75 <  0.90
  expect_equal(run_check("check_frip_retention.R", c(write_rds_tmp(orig), write_rds_tmp(ok))), 0L)
  expect_equal(run_check("check_frip_retention.R", c(write_rds_tmp(orig), write_rds_tmp(bad))), 1L)
})
