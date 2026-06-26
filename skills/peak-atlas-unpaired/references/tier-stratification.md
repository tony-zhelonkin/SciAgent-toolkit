# Tier stratification: selecting the best peak set

After A/B/C are called and merged into a union, each peak must be judged on
*how trustworthy* it is. The source uses a biologically-informed tier scheme
that assigns each union peak to exactly one tier in **priority order** (first
match wins), then sub-stratifies the protected Tier 0. Provenance: archived
`S6.1a_multitier_reproducibility.R` (tiers), `S6.2a_tier0_stratification.R`
(0a/0b/0c). Code: `scripts/tier_stratification.R`.

## Per-peak metadata the tiers key off

Overlap the A/B/C union against each strategy's merged set and aggregate onto
each union peak:
- `strategy_a / strategy_b / strategy_c` — logical, found by that strategy;
  their sum is `cross_strategy_count` (1..3).
- `a_max_n_cells`, `b_max_n_cells` — max cells in any contributing group.
- `b_max_score`, `c_max_score` — Strategy B/C support count (the peak's score).
- `a_n_groups`, `b_n_groups`; `b_n_conditions`, `c_n_conditions`;
  `a_n_clusters`, `b_n_clusters` — breadth counts.
- `has_rare_celltype` — peak belongs to a curated rare type (project-specific
  list, e.g. the BATF3-independent subtype, migratory DCs, a CD209a+ subtype).
- `has_celltype_of_interest` — broader target-celltype membership.

## The tiers (priority order)

```
Tier 0  Biological priority   has_rare_celltype AND
                              ( (b_max_score > 15 & b_n_groups > 0)        # B condition-specific, score>15
                                OR (a_n_groups > 0 & b_n_groups > 0) )      # present in BOTH A and B
Tier 1  Cross-strategy        cross_strategy_count >= 2
Tier 2  High cell fraction    a_max_n_cells >= ~20% of cells   OR  b_max_score > 25
Tier 3a Cross-condition       (b_n_conditions >= 2 & b_max_score > 15)
                              OR (c_n_conditions >= 2 & c_max_score > 15)
Tier 3b Condition-specific    has_celltype_of_interest AND
        target biology        ( (b_n_conditions == 1 & b_max_score > 15 & b_max_n_cells >= 20)
                                OR (a_n_groups == 1 & a_max_n_cells >= 2000) )
Tier 4  Cross-cluster         a_n_clusters >= 2  OR  (b_n_clusters >= 2 & b_max_score > 20)
NotReproducible               everything else  -> EXCLUDED from the atlas
```

The ordering encodes the priority: **rare biology is protected first** (Tier 0
fires before the generic cross-strategy rule), so a rare condition-specific peak
that only one strategy found is not silently demoted to NotReproducible. Tiers 0
through 3b go into the atlas; Tier 4 is a lower-confidence "additions" set;
NotReproducible is dropped.

The "~20% of cells" in Tier 2 is `0.20 * n_cells` (4600 for the ~23k-cell
source) — set it from your own cell count. The `>= 2000` cells in Tier 3b's
Strategy-A arm is the source's threshold for "this lineage is abundant enough
that a single-condition Strategy-A peak is trustworthy."

## Tier 0 sub-stratification (0a / 0b / 0c)

Tier 0 dominates the consensus and mixes confidence levels, so split it by
`cross_strategy_count`:

| Sub-tier | Support | Confidence | Action |
|---|---|---|---|
| **0a** | 3-strategy | HIGH | keep unconditionally — core biology validated across all approaches |
| **0b** | 2-strategy | MEDIUM | keep, but validate separately in the PRIMARY test |
| **0c** | 1-strategy | LOW | keep ONLY if the PRIMARY external-correlation test passes (`r >= 0.4`) |

Tier 0c is the dangerous bucket: single-strategy (usually Strategy-B
condition-specific) rare-cell peaks that are *either* genuine rare biology
*or* noise. In the source these were ~2,464 peaks suspected to be a specific
rare subtype. The decision rule:

```
r >= 0.5   -> keep all of Tier 0 (0a + 0b + 0c); the rare biology is confirmed
r 0.4-0.5  -> keep 0a + 0b, monitor 0c in downstream analyses
r <  0.4   -> keep only 0a, exclude 0b and 0c (likely noise, not rare biology)
```

`keep_tier0c(primary_r, threshold = 0.4)` encodes the conservative gate. See
`references/external-validation.md` for the test itself.

## The modern alternative

A newer pipeline drops explicit tiers entirely. Instead it runs the framework's
**summit-faithful iterative winner-selection** over all A/B/C peaks (rank by
support-boosted score, keep the best per overlap cluster) and records strategy
and source as peak metadata rather than as a tier label
(`peak-atlas-framework/references/iterative-overlap-merge.md`). Prefer this for
new atlases — it is simpler and never produces the confusing empty tiers the
source noted. The tier scheme remains valuable as (a) an auditable, biologically
explicit fallback and (b) the way to read atlases that were built with it.

## Pitfalls

| Symptom | Cause | Fix |
|---|---|---|
| Rare-type peaks land in NotReproducible | tiers assigned without Tier 0 priority, or `has_rare_celltype` not propagated | assign in priority order with Tier 0 first; carry the rare-type flag from Strategy A/B peaks |
| Tier 0c kept and later looks like noise | retained without the correlation gate | gate on PRIMARY `r >= 0.4`; drop 0c if the named hypothesis fails |
| Most peaks fall into one tier | one criterion too permissive (e.g. cell-fraction threshold) | recompute the `0.20 * n_cells` and score thresholds for your dataset size |

## See also

- `references/abc-strategies.md` — where the per-strategy scores come from.
- `references/external-validation.md` — the PRIMARY test that gates 0c.
- `scripts/tier_stratification.R` — `assign_tiers`, `stratify_tier0`,
  `keep_tier0c`.
