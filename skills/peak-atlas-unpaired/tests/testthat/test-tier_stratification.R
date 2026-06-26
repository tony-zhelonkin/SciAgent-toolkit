# Invariants for tier assignment (scripts/tier_stratification.R). Pure GRanges
# logic - fully testable with the synthetic union fixture.
source_script("tier_stratification.R")

KNOWN_TIERS <- c("Tier0_BiologicalPriority", "Tier1_CrossStrategy",
                 "Tier2_HighCellFraction", "Tier3a_CrossCondition",
                 "Tier3b_ConditionSpecificBiology", "Tier4_CrossCluster",
                 "NotReproducible")

test_that("every peak gets exactly one known tier (a partition)", {
  skip_scaffold()
  tiers <- mcols(assign_tiers(make_tier_union()))$tier
  expect_true(all(tiers %in% KNOWN_TIERS))
  expect_length(tiers, 4L)
})

test_that("PRIORITY ORDER: a peak meeting Tier 0 AND Tier 1 is labelled Tier 0", {
  skip_scaffold("first match wins - rare biology is protected before the generic cross-strategy rule")
  tiers <- mcols(assign_tiers(make_tier_union()))$tier
  expect_equal(tiers[1], "Tier0_BiologicalPriority")  # rare + B>15 AND also 3-strategy
  expect_equal(tiers[2], "Tier1_CrossStrategy")       # 3-strategy, not rare
})

test_that("cross_strategy_count == strategy_a + strategy_b + strategy_c", {
  skip_scaffold()
  m <- mcols(assign_tiers(make_tier_union()))
  expect_equal(m$cross_strategy_count,
               as.integer(m$strategy_a) + as.integer(m$strategy_b) + as.integer(m$strategy_c))
})

test_that("stratify_tier0 splits Tier 0 into 0a/0b/0c by support 3/2/1 and partitions it", {
  skip_scaffold()
  out <- assign_tiers(make_tier_union())
  sub <- stratify_tier0(out)
  expect_true(all(mcols(sub$Tier_0a_3strategy_HIGH)$cross_strategy_count == 3))
  expect_true(all(mcols(sub$Tier_0b_2strategy_MEDIUM)$cross_strategy_count == 2))
  expect_true(all(mcols(sub$Tier_0c_1strategy_LOW)$cross_strategy_count == 1))
  n0 <- sum(mcols(out)$tier == "Tier0_BiologicalPriority")
  expect_equal(sum(lengths(sub)), n0)                 # the three buckets partition Tier 0
})

test_that("keep_tier0c gates on the PRIMARY correlation threshold (default 0.4)", {
  skip_scaffold()
  expect_true(keep_tier0c(0.45))
  expect_true(keep_tier0c(0.40))                      # boundary inclusive
  expect_false(keep_tier0c(0.39))
  expect_false(keep_tier0c(0.60, threshold = 0.70))   # custom threshold
})
