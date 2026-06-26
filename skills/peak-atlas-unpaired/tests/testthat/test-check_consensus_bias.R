# Exit-code contract for checks/check_consensus_bias.R.
test_that("check_consensus_bias: PASS <= warn, WARN in (warn, fail], FAIL > fail", {
  skip_scaffold("needs Rscript; contribution vector as RDS")
  pass <- c(d1 = 30, d2 = 25, d3 = 20)   # max 30 <= 40 -> exit 0 (PASS)
  warn <- c(d1 = 45, d2 = 30, d3 = 10)   # 45 in (40, 50] -> exit 0 (WARN)
  fail <- c(d1 = 60, d2 = 30, d3 = 10)   # 60 > 50 -> exit 1 (FAIL)
  expect_equal(run_check("check_consensus_bias.R", write_rds_tmp(pass)), 0L)
  expect_equal(run_check("check_consensus_bias.R", write_rds_tmp(warn)), 0L)
  expect_equal(run_check("check_consensus_bias.R", write_rds_tmp(fail)), 1L)
})

test_that("accepts a two-column CSV (dataset, pct_contribution)", {
  skip_scaffold("needs Rscript")
  p <- tempfile(fileext = ".csv")
  write.csv(data.frame(dataset = c("a", "b"), pct_contribution = c(70, 30)),
            p, row.names = FALSE)
  expect_equal(run_check("check_consensus_bias.R", p), 1L)   # 70% > 50% fail
})
