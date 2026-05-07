"""condition_anova.py — Stage 7, generic per-celltype ANOVA against any factor.

The reference 13403-YD pipeline hardcoded `temp` (33/37/39) as the factor and
Th1/Th17 as the per-celltype iterator. This generalises: factor is named by
the user (Decision Pause 5), celltypes are enumerated from a configurable
column.

Output: program_factor_anova.csv with f_stat, p_value, eta_squared,
highest_at, p_adj (BH-FDR over the entire table).

Optional: pass --covariate <col> to fit `program ~ C(factor) + C(covariate)`
via statsmodels.OLS instead of plain f_oneway. Useful when factor is
confounded with batch.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.stats as stats


def anova(
    in_h5ad: str | Path,
    factor_column: str,
    celltype_column: str = "celltype",
    factor_levels: list[str] | None = None,
    out_csv: str | Path = "program_factor_anova.csv",
    program_prefix: str = "cNMF_",
    covariate: str | None = None,
) -> pd.DataFrame:
    in_h5ad = Path(in_h5ad)
    out_csv = Path(out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(in_h5ad)
    if factor_column not in adata.obs:
        raise KeyError(f"factor_column {factor_column!r} not in obs; "
                       f"available: {list(adata.obs.columns)[:30]}...")
    if celltype_column not in adata.obs:
        raise KeyError(f"celltype_column {celltype_column!r} not in obs")

    program_cols = [c for c in adata.obs.columns if c.startswith(program_prefix)]
    if not program_cols:
        raise RuntimeError(f"no obs columns start with {program_prefix!r}; "
                           "did you run transfer_programs.py first?")

    if covariate is not None:
        try:
            import statsmodels.formula.api as smf  # noqa: F401
        except ImportError as e:
            raise ImportError("--covariate requires statsmodels. `pip install statsmodels`.") from e

    rows: list[dict] = []
    for program in program_cols:
        if adata.obs[program].isna().all():
            continue                          # variant produced no scores; skip

        for ct in sorted(adata.obs[celltype_column].dropna().unique()):
            ct_mask = adata.obs[celltype_column] == ct
            subset = adata.obs[ct_mask]
            if subset.empty:
                continue

            levels = factor_levels or sorted(subset[factor_column].dropna().astype(str).unique())
            groups: list[np.ndarray] = []
            present_levels: list[str] = []
            for lvl in levels:
                vals = subset[subset[factor_column].astype(str) == str(lvl)][program].dropna()
                if len(vals) > 1:
                    groups.append(vals.values)
                    present_levels.append(str(lvl))

            if len(groups) < 2:
                continue

            # plain ANOVA
            f_stat, p_val = stats.f_oneway(*groups)
            all_vals = subset[program].dropna()
            overall_mean = float(all_vals.mean())
            ss_between = float(sum(len(g) * (g.mean() - overall_mean) ** 2 for g in groups))
            ss_total = float(((all_vals - overall_mean) ** 2).sum())
            eta_sq = ss_between / ss_total if ss_total > 0 else 0.0
            means = {lvl: float(subset[subset[factor_column].astype(str) == lvl][program].mean())
                     for lvl in present_levels}
            highest_at = max(means, key=means.get)

            row = {
                "program": program, "celltype": str(ct), "factor": factor_column,
                "levels_present": ",".join(present_levels),
                "f_stat": float(f_stat), "p_value": float(p_val),
                "eta_squared": float(eta_sq), "highest_at": highest_at,
            }

            if covariate is not None:
                # Fit linear model with covariate; compare to factor-only via F-test
                import statsmodels.formula.api as smf
                df = subset[[program, factor_column, covariate]].dropna().rename(
                    columns={program: "y", factor_column: "factor_col", covariate: "cov_col"}
                )
                df["factor_col"] = df["factor_col"].astype("category")
                df["cov_col"] = df["cov_col"].astype("category")
                if df["factor_col"].nunique() < 2 or df["cov_col"].nunique() < 2:
                    row["covariate_p_value"] = np.nan
                else:
                    try:
                        full = smf.ols("y ~ C(factor_col) + C(cov_col)", data=df).fit()
                        reduced = smf.ols("y ~ C(cov_col)", data=df).fit()
                        from statsmodels.stats.anova import anova_lm
                        cmp = anova_lm(reduced, full)
                        row["covariate_p_value"] = float(cmp["Pr(>F)"].iloc[1])
                    except Exception as e:
                        row["covariate_p_value"] = np.nan
                        row["covariate_error"] = str(e)

            rows.append(row)

    if not rows:
        print("WARNING: no (program, celltype) pairs had >=2 levels with cells; nothing to test")
        empty = pd.DataFrame(rows)
        empty.to_csv(out_csv, index=False)
        return empty

    df = pd.DataFrame(rows)
    from statsmodels.stats.multitest import multipletests
    df["p_adj"] = multipletests(df["p_value"].fillna(1).values, method="fdr_bh")[1]
    df.to_csv(out_csv, index=False)
    print(f"  wrote {out_csv} ({len(df)} rows)")

    sig = df[df["p_adj"] < 0.05].sort_values("eta_squared", ascending=False)
    if not sig.empty:
        print(f"  top {min(10, len(sig))} factor-responsive programs:")
        print(sig[["program", "celltype", "eta_squared", "highest_at", "p_adj"]].head(10).to_string(index=False))
    else:
        print("  no programs with p_adj < 0.05; consider whether the chosen factor is the right axis.")
    return df


def _parse() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Per-celltype ANOVA of cNMF programs against a factor.")
    p.add_argument("--in", dest="in_h5ad", type=Path, required=True)
    p.add_argument("--factor", dest="factor_column", type=str, required=True)
    p.add_argument("--celltype-column", type=str, default="celltype")
    p.add_argument("--factor-levels", nargs="+", default=None,
                   help="Subset of levels to include; default uses all observed.")
    p.add_argument("--out-csv", type=Path, required=True)
    p.add_argument("--program-prefix", type=str, default="cNMF_")
    p.add_argument("--covariate", type=str, default=None,
                   help="obs column to include as a covariate (e.g., 'pool')")
    return p.parse_args()


def main() -> int:
    a = _parse()
    anova(
        in_h5ad=a.in_h5ad,
        factor_column=a.factor_column,
        celltype_column=a.celltype_column,
        factor_levels=a.factor_levels,
        out_csv=a.out_csv,
        program_prefix=a.program_prefix,
        covariate=a.covariate,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
