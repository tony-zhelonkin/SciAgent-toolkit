# The four Benzécri properties — operational definitions

This is the reference page for the geometric axis. It defines, in code-citable terms, the four properties classical factor analysis (Benzécri 1973; FactoMineR as ground truth) satisfies and which the modern single-cell stack each ablates. Every property is stated as: what it is, the formula, the FactoMineR code path, why the property matters, what it lets you read off the biplot.

The companion table summarising which codebase preserves which property is in `factor-analysis-framework/SKILL.md` §"The Benzécri four-property frame" and in `families-roster.md`. The reconciliation of these properties with the orthogonal "modelling fitness" axis is in `modelling-vs-map.md`.

Code citations point at FactoMineR R sources; the same math is transcribed in prince (Python) for the PCA / CA / MCA paths.

---

## Property 1 — Duality

**What it is.** A *single* SVD of a metricized data table produces row coordinates and column coordinates in the *same* Euclidean k-dimensional space, sharing axis scaling `sqrt(λ_α)`. The transition formula — "row equals weighted barycenter of its columns, with the eigenvalue as the rescaling constant" — holds quantitatively, not just visually.

**Formula.** Let `X` be the metricized table (centred, scaled, with row weights `r` and column weights `c`). The triplet SVD writes

```
X = U · diag(σ) · V^T,    where     U^T diag(r) U = I,   V^T diag(c) V = I,
```

with eigenvalues `λ_α = σ_α²`. Row coordinates and column coordinates are

```
F = diag(1/sqrt(r)) · U · diag(σ)            # row × k
G = diag(1/sqrt(c)) · V · diag(σ)            # column × k
```

The transition formula:

```
F_{i,α} = (1/sqrt(λ_α)) · Σ_j X_{ij} · c_j · G_{j,α}
G_{j,α} = (1/sqrt(λ_α)) · Σ_i X_{ij} · r_i · F_{i,α}
```

i.e. row `i` lives at the `c_j`-weighted barycenter of its columns, rescaled by `1/sqrt(λ_α)`, and symmetrically for columns. This is duality.

**FactoMineR code path.**
- `R/svd.triplet.R:16-86` — the weighted triplet SVD, producing the single `(U, V, vs)` triplet from one metricized matrix.
- `R/PCA.R:114-115` — row and column coordinates derived from the same triplet sharing `sqrt(eig)`.
- `R/PCA.R:139-140, 162-163` — the supplementary-projection code path, which is the transition formula written explicitly for supplementary rows / columns.
- `R/predict.PCA.R:14` — the transition formula reused for out-of-sample prediction.

**Why it matters.** Without duality, "the gene arrow points at the cell cluster" is metaphor, not measurement. With duality, the angle between a row vector and a column vector encodes the value of that entry in `X` up to the inertia in the chosen axes. Reading a biplot — distance, angle, length — is *only* a valid operation when duality holds.

**Quote.** "`coord.ind` (n × k) and `coord.var` (p × k) live in the **same** k-dimensional Euclidean space, sharing axis scaling `sqrt(eig)`. That co-projection is duality." (`factoMineR-prince.md:§2.1`)

---

## Property 2 — Inertia decomposition

**What it is.** The total weighted variance of the metricized table ("total inertia") is partitioned exactly across the SVD axes. Each axis owns an eigenvalue fraction; each variable on each axis owns a contribution percentage; each sample on each axis owns a contribution percentage; and per-sample / per-variable `cos²` partition the row and column distances across axes.

**Formulas.**

- Total inertia: `TotalInertia = Σ_α λ_α = trace(X^T diag(r) X diag(c))`.
- Per-axis eigenvalue fraction: `λ_α / TotalInertia` (sums to 1 across axes).
- Per-variable contribution per axis: `ctr_{j,α} = c_j · G_{j,α}² / λ_α`; `Σ_j ctr_{j,α} = 1` (or 100%) per axis.
- Per-variable cos² per axis: `cos²_{j,α} = G_{j,α}² / Σ_α G_{j,α}²`; `Σ_α cos²_{j,α} = 1` per variable.
- Per-individual contribution per axis: `ctr_{i,α} = r_i · F_{i,α}² / λ_α`; `Σ_i ctr_{i,α} = 1` per axis.
- Per-individual cos² per axis: `cos²_{i,α} = F_{i,α}² / dist²_i`, where `dist²_i = Σ_α F_{i,α}²` is the squared distance from origin in the full k-space.

In CA the chi-square identity is exact: `sum(inertia.row) == sum(inertia.col) == sum(eig)` (`R/CA.R:102-103`).

**FactoMineR code path.**
- `R/PCA.R:103-110` — eigenvalue table (`$eig`: eigenvalue, percentage, cumulative percentage).
- `R/PCA.R:116-119` — per-variable `contrib`, `cos²`.
- `R/PCA.R:127-128` — per-individual `contrib`, `cos²`, `dist`.
- `R/CA.R:78` — single `base::svd` in CA.
- `R/CA.R:102-103` — the chi-square identity in CA.

**Why it matters.** The inertia budget is the *meaning* of "variance explained" in classical factor analysis. When the budget partitions, you can say: "this axis is 23% of the total variance, and gene X contributes 7% of that axis, and donor Y is 84% explained by these two axes." None of those statements is well-defined when the partition fails. The modern stack mostly retains the *name* "variance explained" and breaks the partition under the hood:
- mixOmics's `prop_expl_var` is *redundancy* `Rd(X, t_h) = (1/p) Σ_j cor²(X_j, t_h)` — a diagnostic, not a partition (`R/explained_variance.R:25`).
- MOFA's per-(view, group, factor) R² does not sum to `r2_total`; `calculate_variance_explained.R:73-76` applies `max(0, .)` precisely because the partition does not hold.
- scITD's per-factor `exp_var` is rank-1 reconstruction loss, made non-additive by ICA rotation.
- sciRED's `pca.explained_variance_ratio_` is computed and never read; FIST min-max-scales `factor_variance` out of recognition.

---

## Property 3 — Categorical-as-barycenter

**What it is.** A categorical variable's *levels* receive coordinates on the same axes as the rows. Each level `ℓ` is placed at the centred-and-scaled barycenter of the active rows that belong to that level. The level is a *point* on the biplot, not a colour, not a heatmap row, not a regression coefficient.

This applies to two cases:
- **Active categoricals** — the categorical participates in the SVD (FAMD, MCA path). Its levels' coordinates are first-class outputs of the fit.
- **Supplementary categoricals** — the categorical is *not* used to fit the axes; it is projected onto them after the fact. This is the diagnostic move: "where do the donor / batch / condition levels land on the F1×F2 plane that gene expression already built?"

**Formula.** For a supplementary categorical with levels `ℓ`, given row coordinates `F` (n × k) and row weights `r`:

```
bary_{ℓ,α} = Σ_{i : label_i = ℓ} (r_i · F_{i,α}) / Σ_{i : label_i = ℓ} r_i
```

i.e. the (row-weight-weighted) mean factor score over the rows in level `ℓ`.

Standardisation: the `v.test` statistic for level `ℓ` on axis `α` is

```
v.test_{ℓ,α} = bary_{ℓ,α} / sqrt( (Var(F_α) / n_ℓ) · ((N - n_ℓ) / (N - 1)) )
```

where `n_ℓ` is the count in level `ℓ` and `N` is the total active sample count. `|v.test| > 1.96` ⇒ level significantly separated from the cloud centroid on that axis.

Per-categorical eta² on axis `α`:

```
eta²_α = Σ_ℓ n_ℓ · bary_{ℓ,α}² / (N · Var(F_α))      # between-group / total variance
```

**FactoMineR code path.**
- `R/PCA.R:187-205` — supplementary categorical: coordinates, `v.test`, eta², dist.
- `R/PCA.R:221` — output object carries `coord, cos2, v.test, eta2, dist` per level.
- `R/FAMD.R:130-152` — active-categorical handling in FAMD (`bary_ℓ · diag(col.w) · V`).
- `R/MCA.R:235` — MCA's level coordinates are the column-points of the disjunctive-table CA by construction.

**Why it matters.** This is the property that fell off the cliff in the modern stack. Categorical metadata is the wet-lab handoff: condition, donor, batch, biopsy site, treatment, response. The biplot's headline question is "where does my categorical fall on the axes my data already drew?" — and that question is *only* answerable when categoricals can be points on the same axes as the rows. Three modern codebases compute the right number and refuse to render it as a point:

- **MOFA**: `summarise_factors` (per-level Z median; `correlate_covariates.R:114-171`) computes the barycenter and renders as a `geom_tile` heatmap. `correlate_factors_with_covariates` computes the supplementary-continuous formula and renders as `corrplot`. The numeric machinery exists; the geometric display layer is amputated.
- **scITD**: `get_meta_associations` (`get_meta_associations.R:16-87`) regresses metadata on donor scores and returns R²/p-values. Per-level group means are computed and discarded.
- **sciRED**: the mean factor score per level is not computed anywhere despite being one line of pandas; categoricals enter as nuisance regressors (`glm.poissonGLM`), supervised classifier targets (FCAT), or cell colourings.

This is what the Level-1 patch (see `supplementary-projection.md`) restores in one afternoon of code per host.

---

## Property 4 — Distribution-free framing

**What it is.** No prior, no likelihood, no ELBO. One call to `base::svd()`. Output is deterministic up to random seed for randomized SVD; otherwise fully reproducible. The metric is built into the matrix (chi-square for indicator blocks, standardised-Euclidean for continuous), *not* into a model of how the data was generated.

This is the operational form of Benzécri's epistemic stance: "Le modèle doit suivre les données, et non l'inverse." — the model must follow the data, not the other way around. The data table is metricized in a way appropriate to its type, the SVD is taken, and you read what is there.

**FactoMineR code path.**
- `R/PCA.R:102` — single `base::svd()` call.
- `R/CA.R:78` — single `base::svd()` in CA.
- `R/svd.triplet.R:16-86` — the entire SVD machinery, with no iterative outer loop.
- `R/FAMD.R, R/MCA.R, R/MFA.R` — all dispatch to `svd.triplet`. No priors, no likelihood, no autograd.

**Why it matters.** The property is what makes the geometry inspectable. When the fit depends on a prior, an init, a seed, an iteration count, and a convergence tolerance, you cannot point at a feature of the biplot and say "this is in the data." You can only say "this is what my model produced given my priors on this random seed." Distribution-free framing is what makes the biplot a *measurement* of the table rather than a *summary of a model fit*.

The trade is that distribution-free framing forbids the moves that bought the modern stack scale and sparsity:
- **MOFA** trades distribution-free framing for ARD sparsity priors, per-view likelihoods, missing-data marginalisation, SVI at 10⁶ cells. This is its architectural axis; the four geometric properties were sacrificed for it (`mofa.md:§2.4`).
- **sciRED** is mixed: the PCA core is distribution-free, but Poisson GLM residualisation per gene is likelihood-based, varimax rotation breaks SVD uniqueness deliberately (in exchange for "interpretable" sparse loadings), and FCAT is a supervised classifier ensemble (`sciRED.md:§2.4`).

For mixOmics, DIABLO, and scITD this property is genuinely preserved — they are deterministic SVD / power-iteration / Tucker-ALS with no Bayesian machinery. This is the corner of the audit table where the modern stack inherits the French school cleanly.

---

## How the four properties interact

The properties are not independent. Distribution-free framing makes duality possible (no iterative outer loop, so the SVD is the whole story). Duality makes the inertia partition exact (`Σ_α λ_α = TotalInertia`). The inertia partition makes the `v.test` standardisation well-defined (it needs `Var(F_α)` per axis). Categorical-as-barycenter is the move that puts metadata onto the *same* axes the rows live on — which is the affordance the other three properties exist to support.

Strip one: the others may survive but their interpretive value degrades. Strip distribution-free (as MOFA does, deliberately): duality is replaced by two separate posterior expectations linked by a likelihood, the partition fails by construction, categoricals lose their natural projection. Strip duality (as sciRED does via varimax): the partition becomes axis-dependent, categorical barycenters live in a rotated basis the user must remember to undo. Strip the partition (as everyone does in practice): the "variance explained" number on the axis becomes a diagnostic rather than a budget.

This is why the patch (`supplementary-projection.md`) is so high-leverage: it restores property 3 on top of any linear scaffold, *and* property 3 carries enough of properties 1–2 with it to make the biplot reading rules work in practice. You do not need to fix MOFA's likelihood to read MOFA's Z geometrically; you need to put the categorical barycenters and the continuous arrows on the same axes the cells already live on.
