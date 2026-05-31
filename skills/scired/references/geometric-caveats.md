# sciRED — geometric caveats against the Benzécri properties

This page summarises sciRED's status against each of the four Benzécri properties that the super-family router (`factor-analysis-framework`) names. The audit is one-line summary + evidence in the code + a "what would have to change" note per property. The long-form audit, with full reasoning, lives at `docs/vision/latent/sciRED.md`.

| Property | sciRED status |
|---|---|
| 1. Duality (rows + columns on the same axes via the transition formula) | **ABSENT** |
| 2. Inertia decomposition (`λ_α`, `ctr`, `cos²` partition) | **PARTIAL / UNUSED** |
| 3. Categorical-as-barycenter | **ABSENT** (closest analogue is variance per level, not mean per level) |
| 4. Distribution-free framing | **MIXED** (PCA is, GLM is not, varimax breaks SVD uniqueness, FCAT is supervised) |

For the operational definitions, defer up to `factor-analysis-framework/references/benzecri-properties.md`. This page only states sciRED's verdict and shows the one-afternoon patches that fix three of the four.

---

## 1. Duality — ABSENT

**One-line verdict.** No co-projection helper exists. Gene scatters and cell scatters live on separate matplotlib figures with no shared rescaling.

**Evidence.**

- `plot_factor_loading` (`utils/visualize.py:70-117`) plots **genes only** on `F_x × F_y` with top-N / bottom-N labels. No arrows from the origin; no inertia-rescaled length; no cell points; no level points.
- `plot_factor_scatter` (`utils/visualize.py:47-66`) plots **cells only** on the same factor pair, coloured by metadata. Separate `plt.figure(...)` invocation — never the same axes.
- The varimax rotation matrix is applied to both loadings and scores (`rotations.py:83-89`, `factor_scores @ rotmat`), but this is a coordinate change for visual consistency, **not** Benzécri's transition formula `F_i = λ_α^{-1/2} Σ_j m_{ij} g_{j,α}` — and the rotation step deliberately breaks SVD uniqueness, so no transition formula can hold.

**What would close the gap.** A `plot_biplot(scores, loadings, x_i, y_i, scaling="symmetric")` helper that:
1. Rescales scores by `1 / √λ_α` (the symmetric biplot convention) — but note `λ_α` is the *unrotated* eigenvalue; after varimax there is no λ.
2. Plots both clouds on one `plt.gca()` with a secondary axis if the magnitudes differ.
3. Optionally draws `geom_segment(0, 0, loading[g, x_i], loading[g, y_i])` for top-N genes.

In sciRED specifically, the cleanest path is to render the **unrotated** scores and loadings as a biplot, then run varimax for the *loading-interpretation* readout in a separate figure. Mixing the two is geometrically incoherent.

---

## 2. Inertia decomposition — PARTIAL / UNUSED

**One-line verdict.** `sklearn.decomposition.PCA.explained_variance_ratio_` is computed but never read. `factor_variance` is reported only in FIST as a min-max-scaled row.

**Evidence.**

- `example_scMixology.py:75-76` does `pca.explained_variance_ratio_` with an immediately-commented-out plot — the value is never consumed by any FCAT or FIST helper.
- `factor_variance` (`metrics.py:96-105`) returns `[np.var(factor_scores[:, i]) for i in range(K)]`. This is the post-rotation column variance — equivalent to the eigenvalue `λ_α` *only before* varimax. After varimax, the per-factor variance is non-orthogonal and non-additive.
- `FIST` includes it as a row labelled (in upstream usage) "Effect size", then min-max-scales it via `get_scaled_metrics` (`metrics.py:355-367`). The min-max scale depends on the *other* factors in the analysis — changing `NUM_COMPONENTS` changes the colour without changing the biology.
- The Benzécri quantities `ctr_{j,α}` (gene contribution to axis α) and `cos²_{j,α}` (axis quality of representation per gene) are not computed anywhere in the package.

**What would close the gap.**

```python
# Unrotated inertia partition
lambda_alpha = pca.explained_variance_              # (K,)
total_var    = lambda_alpha.sum()
pct_inertia  = 100 * lambda_alpha / total_var

# Per-gene contribution to axis α (only valid on UNROTATED loadings)
loadings_unrot = pca.components_                    # (K, n_genes)
ctr_g_alpha    = 100 * loadings_unrot**2 / (loadings_unrot**2).sum(axis=1, keepdims=True)
```

This is a 4-line patch on the unrotated PCA. It cannot be done on the varimax basis because varimax breaks the eigenvalue partition by design.

---

## 3. Categorical-level-as-barycenter — ABSENT

**One-line verdict.** Categoricals enter sciRED as nuisance regressors, supervised classifier targets, or cell colourings — never as barycenters on the cell embedding. The mean factor score per level is *not computed anywhere* despite being one line of pandas.

**Evidence.** Categorical metadata has three roles in sciRED, none of which is "level as a point on the F_x × F_y axes":

1. **Nuisance regressor in `poissonGLM`** (`glm.py:5-42`). Any covariate in the design matrix has its main effect *stripped* from the residuals before PCA — by construction it cannot have an axis in the embedding.
2. **Supervised target in FCAT** (`ensembleFCA.py:198-218`). Each level becomes a binary indicator `y = a_binary_cov`; classifiers learn `X = factor_scores → y`. The output is per-level **classifier importance** vectors, not level coordinates.
3. **Cell colouring** in `plot_factor_scatter` and `plot_pca` (`utils/visualize.py:47-66`, `:8-43`). Levels appear as colours in the cell cloud — exactly the whitepaper's "categorical as colour on a UMAP" anti-pattern.

The closest analogue in the codebase is `get_scaled_variance_level` (`metrics.py:141-152`):

```python
scaled_variance = np.var(a_factor[covariate_vector == covariate_level]) / np.var(a_factor)
```

This is the **dispersion** of a factor *within* cells of a given level — useful for "is this factor uniform within this batch?" — but not the **mean** ("where does this batch sit on this factor?"). The mean factor score per level is the categorical-as-barycenter quantity and it is not computed anywhere.

**What would close the gap (the L1 patch).**

For continuous metadata (~10 lines, generalises `corr.py:4-13`):

```python
def supplementary_continuous(scores, metadata_df):
    out = {}
    for col in metadata_df.select_dtypes(include="number"):
        v = metadata_df[col].to_numpy()
        out[col] = np.array([np.corrcoef(scores[:, k], v)[0, 1]
                             for k in range(scores.shape[1])])
    return out
```

For categorical levels (~10 lines, the move sciRED never makes):

```python
def supplementary_categorical(scores, labels):
    return {lvl: scores[labels == lvl].mean(axis=0) for lvl in pd.unique(labels)}
```

Both can be rendered on the same `plot_factor_scatter` axes — continuous columns as arrows from the origin to `(corr_x, corr_y)`, categorical levels as labelled dots at the barycenter coordinates. The barycenter math is consistent in the post-varimax frame because varimax is orthogonal.

**Two caveats** (the audit at `docs/vision/latent/sciRED.md §6` is explicit):

1. The varimax rotation has chosen its own axis orientation by the time L1 runs. Barycenter coordinates are valid but interpretable only in the rotated basis. This is a labelling note, not a bug.
2. The GLM residualisation step strips out any covariate included in the design matrix. Barycenters of *those* levels will collapse to the origin by construction. If you want to *see* a covariate as a barycenter, leave it OUT of `poissonGLM`'s design matrix.

---

## 4. Distribution-free framing — MIXED

**One-line verdict.** PCA core is distribution-free; pre-PCA residualisation is likelihood-based; varimax breaks SVD uniqueness deliberately; FCAT is supervised and classifier-dependent.

**Evidence per layer.**

| Layer | Where | Distribution-free? |
|---|---|---|
| PCA core | `Pipeline([StandardScaler, PCA(n_components=K)])` in `example_scMixology.py:71` | **Yes** — pure SVD on standardised residuals |
| Pre-PCA Poisson GLM | `sm.GLM(..., family=sm.families.Poisson())` per gene in `glm.py:25` | **No** — likelihood-based; residuals depend on the Poisson assumption |
| Rotation | `varimax` / `promax` in `rotations.py:8-79` | **No** — deliberately non-unique; rotation chosen to maximise the squared-squared-loadings heuristic (`rotations.py:23-35`) |
| FCAT | Ensemble of supervised classifiers in `ensembleFCA.py:99-194` | **No** — model-dependent; importance scores vary with classifier choice and scaling rule |

**Diagnostics that exist.**

- `example_scMixology.py:182-214` measures the cross-correlation of varimax and promax factors as a tell that the rotation is over-determined.
- `permutation.py:176-183` (`shuffle_covariate`) provides an explicit shuffle-the-label null distribution check — the most distribution-aware artefact in the package, but not exercised in any example.

**What this means in practice.**

- If you want a **distribution-free reading** in the strict Benzécri sense (single SVD, sign-equivariant only), drop the GLM step (residualise upstream with scanpy / Seurat) and drop the varimax step (keep the unrotated PCA). What remains is just PCA on Pearson-residual-equivalent data — geometric, unique up to sign.
- If you want the sciRED interpretation engine (FCAT + FIST), accept that the output is supervised and classifier-dependent. Run `time_eff=True` and `time_eff=False` separately and average; cross-validate against `permutation.py`'s shuffle null.

---

## Summary

sciRED's headline interpretive contract is **"factors named by classifier importance against known covariates, with rotated loading scatters as the gene-side readout"** — a different artefact than what a Benzécri biplot delivers. Three of the four properties are either absent or only partially preserved; the fourth (distribution-free) holds only in the PCA core.

The cheapest move that restores three of the four properties on top of sciRED's existing outputs is the **Level-1 supplementary-projection patch**: ~20 lines total to compute per-level barycenters and per-continuous-column correlation arrows, rendered on the same `plot_factor_scatter` axes. See `factor-analysis-framework/references/supplementary-projection.md` for the canonical form and `docs/vision/latent/sciRED.md §6` for the sciRED-specific caveats.
