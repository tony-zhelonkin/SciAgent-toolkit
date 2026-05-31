# Variance explained — what R² is and is NOT in MOFA

The most-cited and most-misread MOFA output. Source: `docs/vision/latent/mofa.md` §2.2.

The headline: **MOFA's per-(view, group, factor) R² is not an inertia partition.** It looks like one (PCA conditioning), tile-plots like one (`plot_variance_explained` in MOFA2), but the budget properties of a Benzécri inertia decomposition do not hold. Read it as a per-view diagnostic, not as a percentage of total variance carved up among factors.

---

## How MOFA computes R²

The engine routine is `BayesNet.calculate_variance_explained` (`mofapy2/core/BayesNet.py:135-176`):

```python
for m in range(self.dim["M"]):
    for g in range(self.dim["G"]):
        SS = np.square(Y[m][gg, :]).sum()
        if total:
            Ypred = np.dot(Z[gg, :], W[m].T)
            Res = np.sum((Y[m][gg, :] - Ypred) ** 2.0)
            r2[g][m] = 1.0 - Res / SS
        else:
            for k in range(self.dim["K"]):
                Ypred = np.outer(Z[gg, k], W[m][:, k])
                Res = np.sum((Y[m][gg, :] - Ypred) ** 2.0)
                r2[g][m, k] = 1.0 - Res / SS
```

The R wrapper (`MOFA2/R/calculate_variance_explained.R:80-84`) mirrors this exactly:

```r
a <- sum((as.matrix(Y[[m]][[g]]) - tcrossprod(Z[[g]][,k], W[[m]][,k]))**2, na.rm = TRUE)
b <- sum(Y[[m]][[g]]**2, na.rm = TRUE)
return(1 - a/b)
```

Per-factor R² is the fraction of `SS(Y)` explained by a rank-1 approximation `Z[, k] · W[m, :, k]^T`, computed **independently against the full Y** — not against the residual after the previous k−1 factors.

---

## Why this is NOT a partition

Four reasons:

### 1. Per-factor R² shares a denominator

`Σ_k r2[g][m, k] ≠ r2_total[g][m]`. Two perfectly correlated factors would each get R² ≈ 0.5 in a true partition; in MOFA's formulation they would each get R² close to the total R², approximately doubling. The `max(0, .)` floor at `MOFA2/R/calculate_variance_explained.R:73-76, 94-97` exists precisely because the formulation can produce *negative* per-factor R² values when factors interfere — which would be impossible in a real partition.

### 2. Z columns are not orthogonal through training

The prior on Z is per-factor independent Gaussian (`Z_nodes.py:158-167`). The posterior is updated via coordinate-wise VI in a loop `for k in range(K)` at `Z_nodes.py:133`. Orthogonality is offered as an **initialisation** option only — `init_model.py:78, 106, 206, 356` show `elif qmean == "orthogonal":` branches — not enforced through training. The MOFA2 docs state at `MOFA2/R/plot_factors.R:471` that "the model encourages the factors to be uncorrelated, so this function usually yields a diagonal correlation matrix" — note the word *usually*. `plot_factor_cor` (`plot_factors.R:487`) exists so users can verify.

Without exact orthogonality, R² values are not additive. They can sum to less than the total (orthogonal factors), more than the total (correlated factors), or anywhere in between.

### 3. Factor ordering is by total R², not eigenvalue

`mofapy2/build_model/save_model.py:115-119`:

```python
order_factors = np.argsort(
    np.array(self.r2).sum(axis=(0, 1), where=~np.isnan(np.array(self.r2)))
)[::-1]
self.r2 = [x[:, order_factors] for x in self.r2]
```

Factors are sorted descending by total R² across groups and views. This *mimics* PCA's eigenvalue ordering but is not the same:

- **PCA:** `λ_1 ≥ λ_2 ≥ ...` is a property of the SVD. Reproducible to numerical precision.
- **MOFA:** the ordering depends on the converged variational posterior, which depends on seed, init, and prior. Two re-fits with different seeds can produce different orderings (and sign flips — see `references/geometric-caveats.md`).

**Implication:** "Factor 3" in one run is not necessarily "Factor 3" in another. Match factors across runs by Tucker congruence on Z, not by index.

### 4. No `cos²` / `ctr` vocabulary

A Benzécri inertia decomposition gives, for each factor α:

- `λ_α` — the inertia contribution of factor α (the budget per axis).
- `cos²_{j, α}` — the *quality of representation* of variable j on axis α; sums to 1 across α.
- `ctr_{j, α}` — the *contribution* of variable j to axis α; sums to 100% across j.

MOFA produces none of these. The variance-explained matrix is `(group × view × factor)` — three-dimensional and useful as a per-view diagnostic ("did factor k capture signal in view m?"), but not a budget. `MOFA2/R/calculate_variance_explained.R:111-...` also exposes `calculate_variance_explained_per_sample`, one step closer to a per-individual budget, but still suffers from non-orthogonality.

---

## Factor dropping (`removeInactiveFactors`)

`mofapy2/core/BayesNet.py:178-220` drops factors whose total R² is below `min_r2` (the `drop_factor_threshold` train option) across all views. This is the **good** use of R²: as a per-view regression diagnostic answering "does this factor capture any signal in any view?". When `drop_factor_threshold` is set (`entry_point.py:867, 938-948`), MOFA can shrink K dynamically.

This is a regression-style filter, not a budget allocation. The dropped factors are not "the bottom 20% of inertia" — they are "the factors that fit no view".

---

## How to read R² safely

**Do:**
- Read the `(view × factor)` heatmap as a diagnostic of factor *activity* per view. "Factor 3 explains 30% of view A and 2% of view B" → factor 3 is largely view-A-specific.
- Use it to identify shared (high R² in many views) vs view-specific (high in one) factors.
- Use it to gauge how much total signal each view contributes to the fit (`Σ_k r2[m, k]` is a rough total-explained-in-this-view, with the caveats above).

**Do not:**
- Present "Factor 1 = 35% of variance" as if it were PCA. Phrase as "Factor 1 explains 35% of view A".
- Sum R² values across factors and call the result "total variance explained".
- Compare R² values across views as if they were on the same scale — `SS(Y[m])` is view-specific.
- Compare R² values across runs without first aligning factors by Tucker congruence.

---

## Reference: equivalent fields in the output

When you read a fitted model:

- **mofapy2 / HDF5:** `model.hdf5/variance_explained/r2_per_factor/<group>` is `(view × factor)`; `model.hdf5/variance_explained/r2_total/<group>` is `(view,)`. They are *not* expected to satisfy `Σ_k r2_per_factor[m, k] == r2_total[m]`.
- **MOFA2 R:** `model@cache$variance_explained$r2_per_factor[[group]]` and `r2_total[[group]]`.
- **mofax Python:** `model.get_r2()` returns a long-format DataFrame.

All three flow from the same engine routine. The semantics are identical.

---

## Cross-reference

- `references/mofa-architecture.md` — why Z is not orthogonal through training.
- `references/geometric-caveats.md` — why no biplot can be built on this R² matrix.
- `references/troubleshooting.md` — "factor pruning unexpected" entry.
