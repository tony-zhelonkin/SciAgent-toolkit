# checks/

Runnable validators for the unpaired atlas. Each exits non-zero on failure so it
can gate a pipeline. The structural / FRiP / marker gates are shared and live in
`peak-atlas-framework/checks/`.

| Check | Asserts |
|---|---|
| `check_consensus_bias.R` | No single external dataset contributes more than the fail threshold (default >50%) of consensus peaks; warns above 40%. Reads `build_external_consensus()$dataset_contribution` (RDS) or a `dataset,pct_contribution` CSV. |

Run: `Rscript checks/check_consensus_bias.R <contribution.rds|csv> [warn_pct] [fail_pct]`
