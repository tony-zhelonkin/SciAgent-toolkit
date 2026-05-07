# Per-celltype ANOVA against a user-named factor

Stage 7 asks: for each celltype × program, does the program score differ across factor levels? The reference 13403-YD uses `temp` (33/37/39 °C) as the factor; your skill must take an arbitrary `obs` column.

## The factor + celltype contract

The user picks (Decision Pause 5):

- `factor_column` — an `obs` column with a discrete factor (string or categorical)
- `celltype_column` — typically `celltype` from upstream annotation
- `factor_levels` — optional subset of levels to include (skip levels with low cell counts)

The ANOVA iterates `(celltype × program)`:

```python
import scipy.stats as stats
import numpy as np

results = []
program_cols = [c for c in adata.obs.columns if c.startswith("cNMF_")]

for program in program_cols:
    if adata.obs[program].isna().all():
        continue                                              # variant produced no scores; skip

    for celltype in sorted(adata.obs[celltype_column].dropna().unique()):
        subset = adata.obs[adata.obs[celltype_column] == celltype]
        if subset.empty:
            continue

        groups, levels_present = [], []
        for lvl in (factor_levels or sorted(subset[factor_column].dropna().unique())):
            g = subset[subset[factor_column] == lvl][program].dropna()
            if len(g) > 1:
                groups.append(g.values)
                levels_present.append(lvl)

        if len(groups) < 2:
            continue                                          # need ≥ 2 levels

        f_stat, p_val = stats.f_oneway(*groups)

        # eta-squared effect size
        all_vals = subset[program].dropna()
        overall_mean = float(all_vals.mean())
        ss_between = sum(len(g) * (g.mean() - overall_mean) ** 2 for g in groups)
        ss_total   = float(((all_vals - overall_mean) ** 2).sum())
        eta_sq = ss_between / ss_total if ss_total > 0 else 0.0

        means = {lvl: float(subset[subset[factor_column] == lvl][program].mean())
                 for lvl in levels_present}
        highest_at = max(means, key=means.get)

        results.append({
            "program": program,
            "celltype": celltype,
            "factor": factor_column,
            "levels_present": ",".join(map(str, levels_present)),
            "f_stat": f_stat,
            "p_value": p_val,
            "eta_squared": eta_sq,
            "highest_at": highest_at,
        })
```

## FDR correction

```python
from statsmodels.stats.multitest import multipletests
results_df = pd.DataFrame(results)
results_df["p_adj"] = multipletests(results_df["p_value"].fillna(1), method="fdr_bh")[1]
```

BH-FDR over the *entire* table (all programs × celltypes × factor combinations). Per-celltype FDR is also a valid choice; document the chosen one in the output schema.

## Effect size — η²

`eta_squared = SS_between / SS_total`. Reportable alongside the p-value; interpretable as "fraction of program-score variance explained by factor". Cohen's small/medium/large thresholds: 0.01 / 0.06 / 0.14.

The reference includes η² in the output table; downstream analysis uses it to rank programs by *biological importance*, not just statistical significance. A program with `p_adj < 0.05` and `η² = 0.001` is highly significant but biologically tiny; a program with `p_adj = 0.01` and `η² = 0.30` matters.

## What to report

Output CSV schema (`program_factor_anova.csv`):

| Column | Type |
|--------|------|
| `program` | str |
| `celltype` | str |
| `factor` | str (the column name) |
| `levels_present` | str (comma-separated) |
| `f_stat` | float |
| `p_value` | float (raw) |
| `eta_squared` | float |
| `highest_at` | str (level name with highest mean) |
| `p_adj` | float (BH-FDR) |

Top responsive programs:

```python
sig = results_df[results_df["p_adj"] < 0.05].sort_values("eta_squared", ascending=False)
print(sig[["program", "celltype", "eta_squared", "highest_at"]].head(10))
```

## Beyond plain ANOVA — when to switch

`scipy.stats.f_oneway` assumes homoscedasticity and independent observations. For unbalanced designs or batch confounds, switch to a linear-model framing:

```python
import statsmodels.api as sm
import statsmodels.formula.api as smf

model = smf.ols(f"{program} ~ C({factor_column}) + C(pool)", data=subset).fit()
print(model.f_pvalue)                                  # global F-test
```

This handles `pool` as a covariate, which `f_oneway` cannot. The shipped script exposes `--with-covariate <col>` for this.

For repeated measures (the same MouseID contributes cells to multiple factor levels — common in within-subject designs), use mixed models (`statsmodels.regression.mixed_linear_model.MixedLM`). Out of scope for the default skill but documented as the upgrade path.

## What "significant" means here

`p_adj < 0.05` in this ANOVA *does not* mean "the program is real" — it means "the program score depends on the factor". A truly real program might be uniformly expressed across all levels (no association); a wholly artefactual program might be batch-confounded with the factor and look highly significant. Cross-reference with the program's classification (`Biological` vs `Technical`) and confidence tier (High vs Low) before drawing conclusions.

## Pitfalls

- **`predicted_doublet` filter applied before this stage.** Already-applied if the QC variant ran; not applied if the raw variant was the one consolidated. Document which.
- **NaN cells.** Programs from a subset variant will be NaN for cells outside the subset; the loop's `groups = [g.dropna() ...]` handles this. Do not impute.
- **Single level present in a celltype.** A small celltype × low-frequency level combo may have only one factor level represented. The loop's `len(groups) < 2: continue` handles this; surface skipped pairs to the user.
