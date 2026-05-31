# sciRED — end-to-end pipeline walkthrough

This page walks the four steps of the canonical sciRED pipeline as documented in `README.md:16-20` and exercised in `example_scMixology.py:42-89`. Each step gets: what it does, the relevant file:line, the arguments that matter, the output shape, and the gotchas. All citations resolve against `/workspaces/DC_hum_verse/01_modules/.ref/sciRED/`.

For the conceptual framing (modelling-vs-map; the four Benzécri properties), defer up to `factor-analysis-framework/references/`. For the per-property audit of sciRED specifically, see `references/geometric-caveats.md`.

---

## Step 0 — Inputs

| Required | Type | Notes |
|---|---|---|
| Raw counts | sparse / dense numeric matrix in `data.X` | cells × genes, **integer** (Poisson assumption) |
| Cell metadata | `data.obs` | At minimum: library size column (e.g. `nCount_RNA`) + one categorical for FCAT |
| Gene names | `data.var_names` | Used by `plot_factor_loading` to label top loaders |

`utils/preprocess.py:46-71` (`get_sub_data`) is the recommended preconditioner: it drops cells with `cell_sums == 0`, drops zero-count genes, and HVG-subsets to `num_genes=2000` by raw-count variance. Note that this HVG rule is **not** Seurat VST / scanpy dispersion — it picks high-mean genes preferentially. For datasets with a wide expression dynamic range, run `scanpy.pp.highly_variable_genes(..., flavor="seurat_v3")` upstream and skip `get_sub_data`.

---

## Step 1 — Log-normalise (optional, upstream-tangled)

`example_scMixology.py:43` does:

```python
y = np.log(y + 1e-10)
```

before passing `y` to `poissonGLM`. This is **not standard** and is at odds with the GLM's own Poisson assumption — `poissonGLM` expects raw integer counts as the response. The `+ 1e-10` offset is numerically unstable and unmotivated.

**Recommendation:** skip this step. Pass raw counts directly to `poissonGLM`. If a downstream consumer needs log-scale, apply `np.log1p` to a copy after the GLM, not before.

If you must replicate the upstream example bit-for-bit, the `1e-10` offset is what `tutorial1_scMixology.ipynb` and the published figures use.

---

## Step 2 — Build a confounder design matrix

Two helpers in `utils/preprocess.py`:

```python
# Library size only
x = proc.get_library_design_mat(data, lib_size="nCount_RNA")  # preprocess.py:119-125
# Stacks: column of 1s (intercept) + column of library sizes
# Shape: (n_cells, 2)

# Categorical covariate (one-hot)
x_protocol = proc.get_design_mat(metadata_col="protocol", data=data)  # preprocess.py:102-115
# Returns one-hot, one column per level. Shape: (n_cells, n_levels)

# Combine (typical: library + categorical + intercept)
x = np.column_stack((data.obs.nCount_RNA, x_protocol))
x = sm.add_constant(x)                                       # adds the intercept
```

**Arguments that matter:**

- `lib_size` keyword must match the column name on `data.obs`. Upstream examples use `nCount_originalexp` (scMixology), `nCount_RNA` (Seurat default), or arbitrary other names — pass the right one.
- One-hot encoding via `get_design_mat` does **not** drop a reference level, so add the intercept explicitly with `sm.add_constant`.

**Output shape:** `(n_cells, n_design_cols)` numpy array.

**Gotcha:** any covariate in `x` is **regressed out before PCA**. If you include `condition` in the design matrix, no factor will track condition — by construction. Decide deliberately:
- Residualise away *technical* confounders (library size, protocol, batch).
- *Leave biological covariates in the metadata only* and let FCAT name the factors that track them.

---

## Step 3 — Per-gene Poisson GLM (the residualisation)

`sciRED/glm.py:5-42`:

```python
glm_fit = glm.poissonGLM(y, x)
# y: (n_cells, n_genes), raw counts as response
# x: (n_cells, n_design_cols), confounder design

resid_pearson = glm_fit["resid_pearson"]   # (n_genes, n_cells) — note the transpose
```

Internally (`glm.py:23-29`):

```python
for i in range(len(y[0])):
    y_a_gene = y[:, i]
    model = sm.GLM(y_a_gene, x, family=sm.families.Poisson())
    result = model.fit()
    resid_pearson.append([result.resid_pearson])
```

**One independent IRLS fit per gene**, sequential, no joint count likelihood. This is the residualise-then-PCA pattern (cf. limma `removeBatchEffect` + PCA), not GLM-PCA.

**Outputs (`glm.py:38-40`):** dict with three keys — `resid_pearson`, `resid_deviance`, `resid_response`. Each is shape `(n_genes, n_cells)`. Only the Pearson residuals are used downstream in the examples.

**Gotchas:**

- **Slow.** ~5–10 minutes for `(5000 cells × 2000 HVGs)` on one CPU; scales linearly with `n_genes`. HVG-subset aggressively.
- **No parallelism in the upstream loop.** A trivial `joblib.Parallel(n_jobs=-1)` wrapping the loop would help, but the package does not provide one.
- **Convergence warnings flood stdout.** statsmodels GLM emits per-gene `IRLS did not converge` for high-zero genes. Suppress with `warnings.catch_warnings()` if running in a notebook.
- **Don't pass log-normalised data here.** The Poisson likelihood assumes integer counts. Step 1 is at odds with this step; honest practice is to skip Step 1.

---

## Step 4 — PCA on residuals + varimax rotation

`example_scMixology.py:71-75` and `:118-121`:

```python
y = resid_pearson.T                          # (n_cells, n_genes)
pipeline = Pipeline([("scaling", StandardScaler()),
                     ("pca",     PCA(n_components=NUM_COMPONENTS))])
pca_scores = pipeline.fit_transform(y)        # (n_cells, K)
pca_loading = pipeline.named_steps["pca"].components_   # (K, n_genes)

# Varimax
vr = rot.varimax(pca_loading.T)               # rotations.py:8-42
varimax_loading = vr["rotloading"]            # (n_genes, K)
pca_scores_varimax = rot.get_rotated_scores(pca_scores, vr["rotmat"])  # rotations.py:83-89
```

**Default in upstream examples (`example_scMixology.py:24`):** `NUM_COMPONENTS = 10`. There is no automatic K selection; sweep manually if needed.

**`StandardScaler` before PCA** — column-wise standardisation. This makes the PCA equivalent to *correlation* PCA on the residual matrix, with each gene re-scaled to unit variance.

**`varimax(x, normalize=True, eps=1e-5)`** is the port of base-R `varimax()`. It iterates up to 1000 times and stops when `d < dpast * (1 + eps)` (`rotations.py:34`). The rotation matrix `vr["rotmat"]` is applied to the scores via `get_rotated_scores` (`rotations.py:83-89`), which is just `factor_scores @ rotmat`.

**Gotcha — varimax is non-unique.** Sign and order flip between runs even with fixed seeds in upstream code. Force a canonical sign by hand if you care about run-to-run comparability:

```python
for k in range(varimax_loading.shape[1]):
    if varimax_loading[np.argmax(np.abs(varimax_loading[:, k])), k] < 0:
        varimax_loading[:, k] *= -1
        pca_scores_varimax[:, k] *= -1
```

Reorder by post-rotation variance:

```python
order = np.argsort(-pca_scores_varimax.var(axis=0))
pca_scores_varimax = pca_scores_varimax[:, order]
varimax_loading = varimax_loading[:, order]
```

**`promax` is available** (`rotations.py:46-79`) as an oblique alternative. `example_scMixology.py:139-141` runs both and computes cross-correlations as a diagnostic — they end up well-correlated, which is mild evidence that the rotation policy is over-determined for the scMixology dataset.

---

## Step 5 — FCAT (factor-covariate association table)

`sciRED/ensembleFCA.py:198-218`:

```python
fcat = efca.FCAT(
    covariate_vec=data.obs["cell_line"],
    factor_scores=pca_scores_varimax,
    scale="standard",        # 'standard' | 'minmax' | 'rank'
    mean="arithmatic",       # 'arithmatic' | 'geometric'
    time_eff=True,           # True → drop RandomForest + KNN-perm
)
# Returns DataFrame: (n_levels, K), rows indexed by covariate level, columns F1..FK
```

**Internally per covariate level:**

1. Build 1-vs-rest binary indicator `a_binary_cov` (`ensembleFCA.py:25-35`).
2. Train an ensemble of classifiers with `X = factor_scores`, `y = a_binary_cov` (`ensembleFCA.py:113-114`):
   - **Default** (`time_eff=True`): LogisticRegression, DecisionTreeClassifier, XGBClassifier (`ensembleFCA.py:99-108`).
   - **Full**: + RandomForestClassifier + KNeighborsClassifier (permutation importance).
3. Extract feature importance (`ensembleFCA.py:117-129`):
   - LogReg → `|coef_|`
   - Tree models → `feature_importances_`
   - KNN → `permutation_importance` (mean)
4. Append Mann-Whitney AUC, rescaled as `2 · |AUC − 0.5|` (`ensembleFCA.py:38-86`).
5. Aggregate model rows into one (1 × K) row via `get_mean_importance_level` (`ensembleFCA.py:142-194`):
   - `scale` rescales each model's importances (`standard` zero-mean unit-var; `minmax` to [0, 1]; `rank` to rank/K).
   - `mean` then arithmetic- or geometric-averages across models.

**Threshold + matched-factor diagnostics (`ensembleFCA.py:221-247`):**

```python
thresh = efca.get_otsu_threshold(fcat.values.flatten())
matched_factors_per_level, pct_factors = efca.get_percent_matched_factors(fcat, thresh)
matched_levels_per_factor, pct_levels  = efca.get_percent_matched_covariates(fcat, thresh)
```

**Gotchas:**

- **FCAT is per-covariate.** For *N* metadata columns, call `efca.FCAT` *N* times and concatenate vertically with a multi-index.
- **Run with both `time_eff=True` and `False`** if classifier robustness matters. The three default classifiers span different inductive biases; their importance scores need not agree.
- **`permutation.py` exists** for shuffle-the-label null distributions (`permutation.py:176-183`, `shuffle_covariate`); not exercised in any example but available for sanity checks.

---

## Step 6 — FIST (factor interpretability scoring tool)

`sciRED/metrics.py:355-383`:

```python
all_metrics = {
    "bimodality_silhouette": met.kmeans_bimodal_score(pca_scores_varimax),
    "bimodality_index":      met.bimodality_index(pca_scores_varimax),     # GMM-based
    "factor_variance":       met.factor_variance(pca_scores_varimax),       # post-rotation var
    "simpson_diversity":     met.simpson_diversity_index(fcat),
    "asv_celltype":          met.average_scaled_var(pca_scores_varimax, data.obs["cell_line"]),
    "dip_statistic":         met.get_dip_test_all(pca_scores_varimax)[0],   # returns (dip, pval); take dip
}
fist = met.FIST(all_metrics)
```

`FIST` constructs a DataFrame from the dict and **min-max-scales each column independently** via `get_scaled_metrics` (`metrics.py:355-367`). The output is intended for `plot_FIST` (`utils/visualize.py:314-349`), a clustermap.

**Gotchas:**

- **Each metric is min-max scaled column-wise.** Adding or removing factors changes the colour scale of the *other* factors without changing the underlying biology.
- **`factor_variance` is post-rotation.** It is not `λ_α` from the original SVD. Varimax preserves total variance of the rotated subspace but redistributes it across columns; the per-factor number is informative as a *rank* but not as a fraction of total inertia.
- **Treat FIST as a visual ranking aid**, not as a calibrated scoreboard. "F7 has high bimodality and high specificity" is a useful triage signal; do not report the FIST cell values as numbers in a paper.

---

## Step 7 — Exporting

`example_scMixology.py:235-245` writes CSVs that are exactly what the L1 supplementary-projection patch (see `factor-analysis-framework/references/supplementary-projection.md`) consumes:

```python
pd.DataFrame(varimax_loading, index=genes,
             columns=[f"F{i}" for i in range(1, K+1)]).to_csv("loadings.csv")

scores_df = pd.DataFrame(pca_scores_varimax, index=data.obs.index,
                         columns=[f"F{i}" for i in range(1, K+1)])
pd.concat([data.obs.reset_index(drop=True),
           scores_df.reset_index(drop=True)], axis=1).to_csv("scores_with_meta.csv")
```

---

## End-to-end verification

After running the pipeline, confirm:

- `pca_scores_varimax.shape == (n_cells, K)` and contains no NaNs.
- `varimax_loading.shape == (n_genes, K)`; rows match `data.var_names`.
- `fcat.shape == (n_unique_levels_in_covariate, K)`.
- The cell scatter `vis.plot_factor_scatter(pca_scores_varimax, x_i=0, y_i=1, cell_color_vec=...)` shows known cell types separating along recognisable axes.
- `vis.plot_factor_loading(varimax_loading, genes, x_i=0, y_i=1)` shows recognisable marker genes at the tips.
- FCAT row for a known biological covariate (e.g. `cell_line`) has one high-importance factor per level, not a diffuse spread.

If FCAT scores are diffuse or factor loadings show no recognisable markers, the most likely cause is **a confounder accidentally regressed in** (Step 2) — re-check that the design matrix only includes technical covariates.
