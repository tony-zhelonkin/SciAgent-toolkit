# peak-atlas-unpaired — test scaffolds

Non-brittle test scaffolds for the unpaired (separate RNA + ATAC pools) peak
logic. They assert **invariants / properties** on synthetic fixtures and are
guarded by `skip_scaffold()` until wired up in an R env. See
`peak-atlas-framework/tests/README.md` for the shared testing philosophy.

## What is covered

| Script / check | Key invariants tested |
|---|---|
| `tier_stratification.R` | tiers **partition** the input; **priority order** (a peak meeting Tier 0 *and* Tier 1 is labelled Tier 0); `cross_strategy_count = a+b+c`; 0a/0b/0c split by support 3/2/1; `keep_tier0c` gate at `r >= 0.4` |
| `build_consensus.R` | `merge_replicates` keeps peaks in `>= ceil(n*0.5)` reps; consensus keeps `>= ceil(n*0.25)`; **two-gap-width rule** (70bp-apart peaks merge at gap=100 but NOT at gap=50); bias gate errors `> 50%`; monotone in threshold; empties dropped |
| `label_transfer_strategy.R` | 501bp/standard-chrom post-process; `group.by` set + min-cells gate + support=#groups (mock CallPeaks); LowConf gate on score/target set (mock TransferData) |
| `checks/check_consensus_bias.R` | `<= warn` → PASS exit 0; `(warn, fail]` → WARN exit 0; `> fail` → exit 1; RDS and CSV inputs |

`merge_replicates`, `build_external_consensus`, and the whole tier scheme are
**pure GRanges logic** and fully testable with the synthetic fixtures here. The
Seurat/MACS3 paths (`transfer_labels`, `call_strategy_peaks`,
`add_gene_activity_bridge`) and the BSgenome N-content filter are mocked /
documented as TODO fixtures.

## The two-gap-width regression

`test-build_consensus.R` encodes the distinction promoted in
`references/external-consensus.md`: `gapwidth = 100` for the cross-dataset union
vs `gapwidth = 50` for within-study (replicate / cross-group) merges. Two peaks
70bp apart must merge under 100 but stay separate under 50 — a direct guard
against silently swapping the two.

## Running

```bash
bash tests/run_skill_tests.sh
```
