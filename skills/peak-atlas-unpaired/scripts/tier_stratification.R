# tier_stratification.R - Tier assignment + Tier 0 sub-stratification for the
#                         unpaired A/B/C consensus.
#
# Provenance (reproduced faithfully):
#   Tier assignment (Tier 0..4, NotReproducible)
#     -> 02_analysis/00_main/.archive/peak_atlas_superseded_20250106/
#        S6.1a_multitier_reproducibility.R  (lines ~398-451)
#   Tier 0 sub-stratification (0a/0b/0c by cross_strategy_count)
#     -> foundation/S6.2a_tier0_stratification.R
#
# Tiers are assigned in PRIORITY ORDER (first match wins). Tier 0 protects rare
# condition-specific biology; Tier 0c (single-strategy) is kept ONLY when the
# external-correlation PRIMARY test passes (see references/external-validation.md).
#
# A newer pipeline drops explicit tiers in favor of the framework's summit-
# faithful iterative winner-selection (peak-atlas-framework/scripts/
# iterative_overlap.R) with strategy/source as metadata; this tier scheme
# remains the auditable fallback and the way to read older atlases.

suppressPackageStartupMessages({
  library(GenomicRanges)
})

#' Assign reproducibility tiers to a union of A/B/C peaks
#'
#' Expects per-peak metadata columns (compute these by overlapping the union
#' against each strategy's merged set; see assemble notes at bottom):
#'   strategy_a, strategy_b, strategy_c : logical, found by that strategy
#'   has_rare_celltype                  : logical, peak belongs to a rare type
#'   has_celltype_of_interest           : logical, peak from a target celltype
#'   a_max_n_cells, b_max_n_cells       : numeric, max cells in any A/B group
#'   b_max_score, c_max_score           : numeric, Strategy B/C support count
#'   a_n_groups, b_n_groups             : numeric, n groups per strategy
#'   a_n_clusters, b_n_clusters         : numeric, n distinct clusters
#'   b_n_conditions, c_n_conditions     : numeric, n distinct conditions
#'
#' @param peaks GRanges union with the columns above.
#' @param high_cell_fraction Cell count for Tier 2 "high fraction" (default
#'   4600 ~= 20% of ~23k cells in the source; set to 0.20 * n_cells for yours).
#' @return GRanges with `tier` and `cross_strategy_count` columns added.
assign_tiers <- function(peaks, high_cell_fraction = 4600) {
  m <- mcols(peaks)
  n <- length(peaks)

  cross_strategy_count <- as.integer(m$strategy_a) +
                          as.integer(m$strategy_b) +
                          as.integer(m$strategy_c)

  tier <- character(n)

  # Tier 0: Biological priority - rare celltype with stringent support
  tier0 <- m$has_rare_celltype & (
    (m$b_max_score > 15 & m$b_n_groups > 0) |   # condition-specific (B), score >15
    (m$a_n_groups > 0 & m$b_n_groups > 0)        # present in BOTH A and B
  )
  tier[tier0] <- "Tier0_BiologicalPriority"

  # Tier 1: Cross-strategy (>= 2 of A/B/C)
  t1 <- tier == "" & cross_strategy_count >= 2
  tier[t1] <- "Tier1_CrossStrategy"

  # Tier 2: High cell fraction OR high Strategy-B score
  t2 <- tier == "" & ((m$a_max_n_cells >= high_cell_fraction) | (m$b_max_score > 25))
  tier[t2] <- "Tier2_HighCellFraction"

  # Tier 3a: Cross-condition (>= 2 conditions & score >15)
  t3a <- tier == "" & (
    (m$b_n_conditions >= 2 & m$b_max_score > 15) |
    (m$c_n_conditions >= 2 & m$c_max_score > 15)
  )
  tier[t3a] <- "Tier3a_CrossCondition"

  # Tier 3b: Condition-specific biology (1 condition, target celltype, score>15, >=20 cells)
  t3b <- tier == "" & m$has_celltype_of_interest & (
    (m$b_n_conditions == 1 & m$b_max_score > 15 & m$b_max_n_cells >= 20) |
    (m$a_n_groups == 1 & m$a_max_n_cells >= 2000)
  )
  tier[t3b] <- "Tier3b_ConditionSpecificBiology"

  # Tier 4: Cross-cluster (>= 2 clusters; score gate for B)
  t4 <- tier == "" & ((m$a_n_clusters >= 2) | (m$b_n_clusters >= 2 & m$b_max_score > 20))
  tier[t4] <- "Tier4_CrossCluster"

  # Remainder: NotReproducible (excluded from the atlas)
  tier[tier == ""] <- "NotReproducible"

  mcols(peaks)$tier <- tier
  mcols(peaks)$cross_strategy_count <- cross_strategy_count
  peaks
}

#' Sub-stratify Tier 0 by cross-strategy support into 0a / 0b / 0c
#'
#' 0a = 3-strategy (HIGH, keep), 0b = 2-strategy (MEDIUM, keep+validate),
#' 0c = 1-strategy (LOW, keep ONLY if external-validation correlation r >= 0.4).
#'
#' @param peaks GRanges already passed through assign_tiers().
#' @return list(Tier_0a_HIGH, Tier_0b_MEDIUM, Tier_0c_LOW) of GRanges subsets.
stratify_tier0 <- function(peaks) {
  t0 <- peaks[mcols(peaks)$tier == "Tier0_BiologicalPriority"]
  csc <- mcols(t0)$cross_strategy_count
  list(
    Tier_0a_3strategy_HIGH   = t0[csc == 3],
    Tier_0b_2strategy_MEDIUM = t0[csc == 2],
    Tier_0c_1strategy_LOW    = t0[csc == 1]
  )
}

#' Decide whether to keep Tier 0c given the PRIMARY correlation result
#'
#' @param primary_r Pearson r from the named external-correlation hypothesis
#'   test (references/external-validation.md).
#' @param threshold Retention threshold (default 0.4).
#' @return TRUE to keep Tier 0c, FALSE to exclude it.
keep_tier0c <- function(primary_r, threshold = 0.4) {
  isTRUE(primary_r >= threshold)
}

# ---------------------------------------------------------------------------
# Assembling the per-peak metadata (provenance: archived S6.1a):
#   For each strategy's merged GRanges, findOverlaps(union, strategy_merged),
#   then aggregate score_max (max), n_groups / n_conditions / n_clusters, and
#   max n_cells onto the union peak. strategy_a/b/c are overlapsAny() flags.
#   has_rare_celltype / has_celltype_of_interest come from the per-group celltype
#   labels carried on Strategy A/B peaks (the rare-type list is project-specific).
# ---------------------------------------------------------------------------
