# Supplementary projection on Z — the L1 patch

The mitigation when MOFA is the upstream and the consumer wants a Benzécri-style biplot. Source: `docs/vision/latent/mofa.md` §5–§6 and `docs/vision/latent/synthesis.md` §5.

**Headline:** MOFA's Z is approximately PCA-shaped (Gaussian prior, linear bilinear coupling). The Level-1 supplementary-projection patch — continuous variables as arrows by Pearson correlation with each factor, categorical levels as labelled points at the barycenter of their cells — works on Z as-is. The numeric machinery already exists in MOFA2 as `correlate_factors_with_covariates` and `summarise_factors`; what is missing is the *display* layer.

This file gives:
- The math (whitepaper §1.4 → operationalised here).
- An R sketch (~50 lines) on top of `MOFA2::plot_factors`.
- A Python sketch (~50 lines) using `mofax.get_factors`.
- The caveats that determine whether the patch is faithful.

---

## 1. The math

Let `Z ∈ R^{N × K}` be the MOFA factor matrix (samples × factors). Given α, β ∈ {1, ..., K} you want to render `F_α × F_β`.

### 1.1 Continuous variable v → arrow

For each continuous metadata column `v ∈ R^N`:

```
g_α = cor(Z[:, α], v),  g_β = cor(Z[:, β], v)
arrow_tip = (g_α, g_β)
```

This is Pearson correlation of the variable with each factor. `MOFA2::correlate_factors_with_covariates` (`MOFA2/R/correlate_covariates.R:23-95`) already computes the full `K × p` matrix; the patch is `geom_segment` from origin to `(g_α, g_β)` overlaid on `plot_factors`'s scatter.

### 1.2 Categorical level ℓ → labelled point

For each level ℓ of a discrete grouping `labels ∈ {1, ..., L}^N`:

```
c_α(ℓ) = mean(Z[labels == ℓ, α])
c_β(ℓ) = mean(Z[labels == ℓ, β])
label_position = (c_α(ℓ), c_β(ℓ))
```

This is the barycenter of cells with label ℓ. `MOFA2::summarise_factors` (`MOFA2/R/correlate_covariates.R:114-171`) already computes the per-level **median** (line 132 — `dplyr::summarise(value = median(value, ...))`); the patch is `geom_text` at the barycenter overlaid on `plot_factors`'s scatter. Mean vs median is a minor choice — use mean for Benzécri orthodoxy, median for outlier robustness.

The full Benzécri form includes a `1/√λ_α` rescaling so the level point sits at the same scale as the factor's variance. For MOFA there is no `λ_α`; the closest substitute is `sqrt(var(Z[:, α]))`. Whether to apply the rescaling is the first caveat below.

---

## 2. R sketch — `mofa_biplot()` on top of `plot_factors`

```r
library(MOFA2)
library(ggplot2)
library(dplyr)

mofa_biplot <- function(
  model,
  factors = c(1, 2),
  color_by = NULL,
  arrows_continuous = NULL,    # character vector of obs columns
  points_categorical = NULL,   # character vector of obs columns
  arrow_scale = 1.0,
  arrow_color = "#444444",
  label_color = "#222222"
) {
  # Base scatter
  p <- plot_factors(model, factors = factors, color_by = color_by, dot_size = 1)

  Z <- get_factors(model, factors = factors)[[1]]   # N × 2
  meta <- samples_metadata(model)

  # Continuous arrows
  if (!is.null(arrows_continuous)) {
    for (v_name in arrows_continuous) {
      v <- as.numeric(meta[[v_name]])
      keep <- complete.cases(v, Z)
      g_alpha <- cor(Z[keep, 1], v[keep])
      g_beta  <- cor(Z[keep, 2], v[keep])
      p <- p +
        geom_segment(
          aes(x = 0, y = 0, xend = arrow_scale * g_alpha,
              yend = arrow_scale * g_beta),
          arrow = arrow(length = unit(0.18, "cm")),
          color = arrow_color, inherit.aes = FALSE
        ) +
        annotate("text",
          x = arrow_scale * g_alpha * 1.1,
          y = arrow_scale * g_beta  * 1.1,
          label = v_name, color = arrow_color, size = 3
        )
    }
  }

  # Categorical level points
  if (!is.null(points_categorical)) {
    for (cat_name in points_categorical) {
      labels <- meta[[cat_name]]
      cents <- data.frame(Z) %>%
        mutate(level = labels) %>%
        filter(!is.na(level)) %>%
        group_by(level) %>%
        summarise(across(everything(), mean), .groups = "drop")
      colnames(cents)[2:3] <- c("x", "y")
      p <- p +
        geom_text(data = cents,
                  aes(x = x, y = y, label = level),
                  color = label_color, size = 3.5,
                  fontface = "bold", inherit.aes = FALSE)
    }
  }

  p
}

# Usage:
# mofa_biplot(model, factors = c(1, 2),
#             color_by = "condition",
#             arrows_continuous = c("age", "qc_pct_mito"),
#             points_categorical = c("cell_type", "donor"))
```

---

## 3. Python sketch — using `mofax`

```python
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import mofax as mfx

def mofa_biplot(
    model_path,
    metadata,                    # pd.DataFrame, indexed by sample
    factors=(0, 1),              # zero-indexed for Python
    color_by=None,
    arrows_continuous=None,      # list of metadata column names
    points_categorical=None,     # list of metadata column names
    arrow_scale=1.0,
    figsize=(6, 6),
):
    model = mfx.mofa_model(model_path)
    Z = model.get_factors(factors=list(factors))      # N × 2 numpy array
    samples = model.samples                            # alignment with metadata

    df = pd.DataFrame(Z, index=samples, columns=[f"F{f+1}" for f in factors])
    df = df.join(metadata, how="inner")

    fig, ax = plt.subplots(figsize=figsize)
    c = df[color_by] if color_by else None
    ax.scatter(df.iloc[:, 0], df.iloc[:, 1], c=c, s=6, alpha=0.6)
    ax.axhline(0, color="lightgray", lw=0.5)
    ax.axvline(0, color="lightgray", lw=0.5)

    # Continuous arrows
    if arrows_continuous:
        for v_name in arrows_continuous:
            v = pd.to_numeric(df[v_name], errors="coerce")
            keep = v.notna()
            g_a = np.corrcoef(df.iloc[keep.values, 0], v[keep])[0, 1]
            g_b = np.corrcoef(df.iloc[keep.values, 1], v[keep])[0, 1]
            ax.annotate(
                "", xy=(arrow_scale * g_a, arrow_scale * g_b),
                xytext=(0, 0),
                arrowprops=dict(arrowstyle="->", color="#444"),
            )
            ax.text(arrow_scale * g_a * 1.08, arrow_scale * g_b * 1.08,
                    v_name, color="#444", fontsize=9)

    # Categorical level points
    if points_categorical:
        for cat_name in points_categorical:
            cents = df.groupby(cat_name).agg({df.columns[0]: "mean",
                                              df.columns[1]: "mean"})
            for lvl, row in cents.iterrows():
                ax.text(row.iloc[0], row.iloc[1], str(lvl),
                        color="#222", fontsize=10, fontweight="bold")

    ax.set_xlabel(f"Factor {factors[0]+1}")
    ax.set_ylabel(f"Factor {factors[1]+1}")
    return fig, ax

# Usage:
# fig, ax = mofa_biplot(
#     "model.hdf5", metadata=adata.obs,
#     factors=(0, 1), color_by="condition",
#     arrows_continuous=["age", "qc_pct_mito"],
#     points_categorical=["cell_type"],
# )
```

---

## 4. Caveats — what makes the patch (un)faithful

### 4.1 Sign convention

Re-fitting MOFA can flip the sign of any factor (`references/geometric-caveats.md` §5.1). The patch is sign-naive: an arrow pointing right in run 1 may point left in run 2. **Impose a canonical sign** before rendering — the simplest is largest-absolute-loading-positive across views. Document the convention in the figure caption.

### 4.2 The `(N - n) / (N - 1)` v.test correction is ambiguous

FactoMineR's classical supplementary-projection rendering normalises a level-point coordinate by a v.test:

```
v.test_α(ℓ) = (mean(Z[labels == ℓ, α]) - mean(Z[:, α])) /
              sqrt(var(Z[:, α]) * (N - n_ℓ) / (n_ℓ * (N - 1)))
```

For a Benzécri SVD, `var(Z[:, α])` is `λ_α`. For MOFA there is no `λ_α`; substituting `var(Z[:, α])` from the variational posterior is a *defensible approximation* but not the canonical Benzécri object. **Recommendation: drop the v.test correction in the L1 patch.** Render the raw barycenter `mean(Z[labels == ℓ, α])` and document that the rendering is a Z-space barycenter, not a v.test. If you want significance, do `MOFAcellulaR::get_associations` separately.

### 4.3 Stability gate attaches to the upstream fit, not the projection

The L1 patch is mechanically valid at any N. It is *interpretively* meaningless when the upstream MOFA fit is unstable (Tucker congruence on Z `< 0.85` across re-fits — see `synthesis.md:27`). MOFAcellulaR's typical N (10–100 donors) is in the danger zone (`mofa.md` §4). **Recommendation: gate biplot rendering on a stability check at the MOFA layer.** Re-fit MOFA with `B = 20` bootstrap re-seedings; abort the biplot if median Tucker congruence < 0.85.

### 4.4 High-cardinality categoricals saturate the plot

A donor metadata column with 50 levels produces 50 labelled points clouding the scatter. **Cap before projection** — render only levels with `n_ℓ >= 10` cells, or only the top-K levels by sample count, or aggregate rare levels into "other".

### 4.5 UMAP-on-Z compounds the problem

`MOFAcellulaR::plot_sample_2D` runs UMAP on Z (`MOFAcellulaR/R/plot_sample_2D.R:89` — `uwot::umap(factors)`). This compounds the duality loss: even the linear Z is now displayed non-linearly, breaking any chance of arrows working on the rendered axes. **The L1 patch must be applied to the linear F_α × F_β scatter, never to a UMAP of Z.**

### 4.6 What MOFA's `correlate_factors_with_covariates` already does (and doesn't)

`MOFA2/R/correlate_covariates.R:23-95` computes `cor(Z[:, α], v)` for every continuous metadata column and renders it as a `corrplot` heatmap (line 78) or BH-adjusted p-value heatmap (line 89). The **numeric output is the same** as the arrow tip in step 1.1 — only the display differs.

`MOFA2/R/correlate_covariates.R:114-171` (`summarise_factors`) does the median-per-level computation and renders it as a `geom_tile` heatmap. Same numeric, different display.

So the L1 patch is **purely a display refactor** on existing numerics. No new statistics are introduced.

---

## 5. When NOT to use this patch

- **You need a true Benzécri biplot with `cos²` / `ctr` budgets.** The patch gives a *map overlay* on a *modelled* Z; it does not give you inertia decomposition. Route to `factor-analysis-framework` for FAMD/MFA, which gives the budgets natively.
- **The upstream MOFA fit is unstable.** See caveat 4.3.
- **You are rendering a UMAP of Z.** See caveat 4.5. Apply L1 only to the linear factor scatter.

---

## Cross-reference

- `references/geometric-caveats.md` — why MOFA does not give you a biplot natively, sign / order indeterminacy.
- `references/variance-explained.md` — why there is no `λ_α` to plug into v.test.
- The parent `factor-analysis-framework` — for the trade-off between modelling fitness (MOFA) and map fitness (FAMD/MFA).
