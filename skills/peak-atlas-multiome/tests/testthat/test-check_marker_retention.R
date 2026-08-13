# Exit-code contract for checks/check_marker_retention.R. The check needs an
# EnsDb for gene coordinates, so this stays a documented scaffold: build a
# synthetic atlas whose peaks overlap known marker promoters, then verify the
# gate. The assertions below are the spec.
test_that("check_marker_retention passes >= threshold and fails below / on a lost panel", {
  skip_scaffold("needs Rscript + an EnsDb + GRanges atlases whose peaks overlap marker promoters")
  # Fixture plan (in an R env with EnsDb.Mmusculus.v79):
  #   proms <- promoters(genes(ensdb)[symbol %in% KEY_MARKERS], -2000, 500)
  #   orig  <- peaks covering ALL marker promoters
  #   full  <- orig                         # retains everything
  #   lossy <- orig minus > 5% of markers   # drops below 0.95
  #   pdc0  <- orig minus the whole pDC panel
  # Expected exit codes:
  #   run_check("check_marker_retention.R", c(orig_rds, full_rds))  == 0
  #   run_check("check_marker_retention.R", c(orig_rds, lossy_rds)) == 1   # overall < 0.95
  #   run_check("check_marker_retention.R", c(orig_rds, pdc0_rds))  == 1   # a panel hit zero
  #   run_check("check_marker_retention.R", c(empty_rds, full_rds)) == 1   # no markers in ORIGINAL
  succeed()  # placeholder so the scaffold is a valid, skipped test
})
