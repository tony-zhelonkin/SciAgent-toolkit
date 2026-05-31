# Level-1 supplementary projection — the cheapest, highest-leverage patch

This reference page operationalises the Level-1 patch named in synthesis §5 and §6: project arbitrary metadata (continuous variables and categorical levels) onto an existing linear embedding, on the same axes as the sample scatter. Formulas, code receipts from the audited codebases, effort estimates, and the open questions that need adjudication before a host-neutral helper ships.

For the *conceptual* motivation (why this restores three of the four Benzécri properties on top of any honest fitter) see `benzecri-properties.md` §"Property 3" and `modelling-vs-map.md` §4. For per-family detail (which codebase already computes which ingredient) see `families-roster.md`.

---

## 1. The two formulas

### 1.1 Continuous arrow

For a continuous metadata variable `v` (length n) and an embedding `F: n × K`, place an arrow on the F1×F2 scatter with tip at

```
g_α = cor(v, F_α)         for α = 1, 2
tip = (g_1, g_2)
```

i.e. the Pearson correlation of `v` with each factor. The arrow lives inside the unit circle by construction; the closer the tip is to the circumference, the more of `v`'s variance is captured by the two factors plotted.

This is the "correlation circle" reading from MCA / FAMD / PCA. It is mathematically the supplementary-continuous projection from `R/PCA.R:139-140`.

### 1.2 Categorical level point

For a categorical metadata variable with levels `ℓ ∈ {ℓ_1, …, ℓ_L}` and an embedding `F: n × K`, place each level at the (optionally weighted) barycenter of the rows in that level:

```
bary_{ℓ,α} = Σ_{i : label_i = ℓ} (r_i · F_{i,α}) / Σ_{i : label_i = ℓ} r_i
```

where `r_i` is the row weight (defaults to `1/n`). For an unweighted embedding this is just

```
bary_{ℓ,α} = mean(F[labels == ℓ, α])
```

Optional rescaling by `1/√λ_α` matches FactoMineR's MCA convention; for non-FactoMineR embeddings the rescaling is conventional and should be documented per host.

### 1.3 `v.test`

The standardisation statistic for a categorical level on a single axis:

```
v.test_{ℓ,α} = (bary_{ℓ,α} - 0) / sqrt( (Var(F_α) / n_ℓ) · ((N - n_ℓ) / (N - 1)) )
```

where `n_ℓ` is the count in level `ℓ` and `N = Σ_i 1` is the total active sample count. `|v.test_{ℓ,α}| > 1.96` ⇒ level `ℓ` is significantly separated from the cloud centroid on axis `α` at α = 0.05.

The `(N - n_ℓ) / (N - 1)` factor is the finite-population correction. See `R/PCA.R:214` for the FactoMineR implementation; see Open Question §6.1 below for the question of whether to keep it for non-FactoMineR embeddings.

### 1.4 `eta²`

The between-group / total variance ratio for a categorical on a single axis:

```
eta²_α = Σ_ℓ n_ℓ · bary_{ℓ,α}² / (N · Var(F_α))
```

i.e. fraction of factor variance explained by the categorical. Reportable in a small table alongside the biplot.

### 1.5 `cos²` for samples and variables

Per-sample per-axis `cos²`:

```
cos²_{i,α} = F_{i,α}² / dist²_i,        where      dist²_i = Σ_β F_{i,β}²
```

i.e. the fraction of sample `i`'s squared distance from origin that lives on axis `α`. Equivalent for variables, with `G` in place of `F`.

`cos²` summing across axes is 1 per sample / per variable. Useful for "is this dot near the plane, or projected onto it from far away?" diagnostics — a high-`cos²` point can be read; a low-`cos²` point cannot.

---

## 2. Where the machinery already exists (audit receipts)

These are the per-host code paths that already compute the right numbers. The renderer that overlays them on the sample scatter is the part nobody wrote.

| Host | Continuous arrow | Categorical barycenter | Notes |
|---|---|---|---|
| **MOFA2** | `correlate_factors_with_covariates` — computes `cor(v, F_α)` per factor, renders as `corrplot` | `summarise_factors` (per-level Z median) — `correlate_covariates.R:114-171`, renders as `geom_tile` heatmap | Numbers exist; only `geom_segment` + `geom_text` overlay missing. Effort: half a day. |
| **mixOmics** | `cor(X, object$variates$X)` — one line; analogous to `plotVar` correlation-circle math | per-level mean of `object$variates` — not computed natively | Ingredients on every fit object; ~50-line helper adds `supvar.pca(object, factor)`. |
| **DIABLO** | per-block `cor(X[[b]], variates[[b]])` — analogous | per-level mean of `object$variates[[block]]` | Caveat: when projected metadata correlates with supervised Y, projection is mechanically valid but interpretively degenerate. Document. |
| **scITD** | `get_meta_associations` (`get_meta_associations.R:16-87`) — already runs the linear model; group means computed and discarded | same — already there, exported as R²/p-values | 20-line patch: export per-level means as coordinates instead of regression statistics. |
| **sciRED** | `get_factor_libsize_correlation` exists — generalise to any continuous metadata, ~10 lines | `pca_scores_varimax[labels == ℓ].mean(axis=0)` — one line of pandas | Caveats: varimax-rotated axes; covariates included in the GLM design matrix project to origin. |
| **FactoMineR** | first-class via `quanti.sup` (`R/PCA.R:187-205`) | first-class via `quali.sup` (`R/PCA.R:187-205`) | This is the ground truth. The math here is `R/PCA.R:139-205`. |
| **prince** | `pca.column_correlations` exists for PCA | partial — eta² for active categoricals in FAMD only | No `quali.sup` / `quanti.sup` analog in prince. |

---

## 3. Effort estimates

| Target | Effort | Why |
|---|---|---|
| Per-host renderer (one host) | ~50 lines | `geom_segment` + `geom_text` overlay on existing scatter; numbers already computed. |
| Host-neutral helper | ~30 lines of math + 100 lines typed wrapper | `project_metadata_onto_embedding(Z, metadata_df, row_w=None)` returning `{coord, cor, cos2}` per continuous and `{coord, cos2, v.test, eta2, dist}` per categorical level. |
| MOFA-specific overlay on `plot_factors` | half a day | `MOFA2::correlate_factors_with_covariates` + `MOFA2::summarise_factors` ⇒ one ggplot layer. |
| Full Level-1 patch across all four implemented MOFA repos | ~1 week | Math + typed wrapper + one renderer + provenance tagging + docs. |

The math is 30 lines transcribed from `R/PCA.R:139-205`. The typed-output wrapper, sign-canonicalisation, weight handling, and per-host adapters are the bulk of the work.

---

## 4. The host-neutral helper — recommended API

```python
project_metadata_onto_embedding(
    Z: ndarray,                   # n × K, the linear embedding
    metadata_df: pd.DataFrame,    # n rows; columns are either numeric (continuous) or categorical
    row_w: ndarray | None = None, # n; defaults to uniform 1/n
    rescale_categorical: Literal["none", "inv_sqrt_lambda"] = "none",
    lambda_per_axis: ndarray | None = None,
) -> LatentSpaceOverlay
```

`LatentSpaceOverlay` is a typed artifact carrying:

- `continuous: pd.DataFrame` indexed by variable name, columns `coord_axis_1, coord_axis_2, …, cor_axis_1, …, cos2_axis_1, …`
- `categorical: pd.DataFrame` MultiIndex on (variable, level), columns `coord_axis_1, …, cos2_axis_1, …, v_test_axis_1, …, n_level, dist`
- `eta2: pd.DataFrame` indexed by variable, columns `axis_1, axis_2, …`
- `provenance: dict` carrying the upstream embedding identifier, the row-weight scheme, the rescaling choice, the seed if randomized, and any host-specific caveats (varimax basis, supervised Y, mean-encoder linearity assumption).

The renderer is a separate function consuming the overlay + the sample scatter axes.

---

## 5. Caveats per host

- **PCA / FAMD / MFA from FactoMineR or scanpy:** straightforward; the linearity assumption holds.
- **MOFA's Z:** linearity holds (Gaussian prior, linear `Y ≈ ZW^T`). Sign / order of factors is non-canonical across runs; the overlay must enforce a sign convention at the L1 boundary.
- **scVI mean encoder:** approximately linear in a neighborhood (whitepaper §6.1). Use the mean encoder output, not the full encoder. Document the linearity-in-neighborhood caveat.
- **UMAP:** do NOT project. UMAP coordinates are nonlinear; barycenter / correlation interpretation is invalid.
- **sciRED (varimax PCA):** valid because varimax is an orthogonal rotation — barycenters and correlations transform consistently. But the axis interpretation is the *rotated* basis, not the principal-variance basis; document this.
- **DIABLO (supervised hosts):** projection of metadata correlated with Y is mechanically valid but interpretively degenerate. Emit a per-host warning when the projected metadata correlates above some threshold with the supervised driver.
- **scITD's donor scores:** valid; Tucker-2 produces a proper linear donor embedding. ICA rotation is orthogonal, so the math transfers; the per-factor `exp_var` non-additivity is a separate concern (see `benzecri-properties.md` Property 2).

---

## 6. Open questions (synthesis §9)

These are unresolved and need adjudication before a host-neutral helper ships.

### 6.1 The `(N − n_ℓ) / (N − 1)` correction for non-FactoMineR embeddings

FactoMineR's `v.test` standardisation includes a finite-population correction that applies cleanly to *active* analyses (where the row weights sum to 1 and the embedding was fit on the same N). When projecting a supplementary categorical onto an *externally-fitted* embedding (MOFA Z, scVI mean), the right finite-population correction is ambiguous. The synthesis recommendation is "drop the correction and document"; this needs methodological resolution before publication.

### 6.2 Sign / rotation indeterminacy on cross-run comparison

MOFA's factor sign is non-canonical across runs. sciRED's varimax basis is non-unique. scITD's ICA rotation breaks variance ordering. Any cross-dataset L1 comparison requires Procrustes alignment, which no audited codebase ships. The overlay must enforce a sign convention at the L1 boundary so re-fits do not produce sign-flipped biplots; the convention itself (sign of the loading of a chosen high-`cos²` variable? sign of the first row?) needs to be picked.

### 6.3 High-cardinality categoricals

Single-cell metadata routinely has 10² donors and 10⁴ cell-barcode-derived levels. The chi-square `1/√p_k` rescaling penalises rare levels but does not prevent the dense indicator block from being p × K wide. The L1 helper must enforce a cardinality limit (~50 levels?) for supplementary categoricals and declare the limit in the API.

### 6.4 Cross-host comparability

If pathway-explorer ships L1 as a host-neutral helper, the same metadata column projected onto MOFA Z, scVI mean, and scanpy PCA will produce three different biplots. Which is the "true" one? The synthesis claim is that all three are valid in their own basis but not directly comparable — and that is an interpretive principle the platform must teach, not engineer around. The overlay's `provenance` block should carry the upstream identifier so a downstream consumer cannot accidentally compare biplots across hosts.

### 6.5 Indicator rescaling drift between FactoMineR and prince

`factoMineR-prince.md §3.3` identifies that FactoMineR uses `prop` and prince uses `prop * 2` as the categorical-block rescaling constant. The L1 helper should standardise on FactoMineR's `prop` and document the choice.

### 6.6 Stability metadata inheritance

L1 inherits whatever stability the upstream embedding provides. The Tucker congruence < 0.85 stability gate from `latent-lens/synthesis.md` applies to F_1, F_2 *as fitted objects*; once accepted, the projection of `v` and `ℓ` onto them is a 2-line computation. The platform must carry an embedding-level stability score on its `LatentSpaceInput` so L1 overlays inherit it without re-computing — but the *computation* of that score is method-specific and not addressed by any audited codebase.

---

## 7. Why this is the highest-leverage move in the sequence

Nobody in single-cell ships this. The numeric machinery is hiding inside `summarise_factors`, `correlate_factors_with_covariates`, `get_meta_associations`, `get_factor_libsize_correlation` — the renderer that puts them on the same axes as the cell / sample scatter is the part nobody wrote. The intellectual leverage per line of code is the highest in the whole synthesis sequence (L1 < L2 < L3 < L4 < L5):

- **L1 effort:** half-day to one week.
- **L1 unlocks:** the entire biplot reading vocabulary the field forgot — distance, angle, arrow length, level separation, factor-variable association on the same axes.
- **L1 risk:** linearity-of-scaffold assumption + sign convention + cardinality cap. All manageable per §5–§6 above.

It is the move synthesis §6 recommends as the H1 deliverable for the latent-lens layer of pathway-explorer, pivoting the H1 scope from "build a new linear factor model + biplot" to "add supplementary-projection arrows + level points to whatever scatter the team is already rendering."
