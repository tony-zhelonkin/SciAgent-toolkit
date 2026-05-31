# Geometric caveats — duality, biplots, sign / order ambiguity

This file collects the four failure modes against the Benzécri-school properties. Source: `docs/vision/latent/mofa.md` §2.1–2.4 and the linked whitepaper `docs/vision/research/The Geometry We Left Behind.md`.

The audit verdict: MOFA fails **all four** Benzécri properties. Each failure is structural, not an implementation oversight, and follows from the Bayesian generative framing. The mitigation for the most common consumer need (a biplot) is the L1 supplementary-projection patch in `references/supplementary-projection-on-z.md`.

---

## 1. Duality — ABSENT

Duality in the Benzécri sense requires that rows and columns be co-projected into the *same* Euclidean coordinate system via a single SVD `M = UΣV^T`, with a transition formula `row_i = (1/√λ_α) Σ_j m_ij · col_j`.

MOFA has none of this. Z and W are two separate posterior expectations of two latent variables in a bilinear Gaussian likelihood. Five concrete consequences:

### 1.1 No shared singular values

A real SVD distributes scale via a single Σ; row coords `UΣ` and column coords `VΣ` share that Σ, so distances are commensurable. MOFA has no Σ. The R extractors apply independent rescalings:

- `get_factors(model, scale=TRUE)` divides Z by `max(abs(Z))` (`MOFA2/R/get_methods.R:209`).
- `get_weights(model, scale=TRUE)` divides each view's W by `max(abs(W))` (`MOFA2/R/get_methods.R:266`).

The two scalings are cosmetic normalisations, not the two halves of a unique decomposition. A sample at Z-coord `(1.0, 0.5)` and a feature at W-coord `(1.0, 0.5)` share factor *labels* F1 and F2 but not a metric. Their proximity is not interpretable.

### 1.2 No transition formula

The Benzécri identity "a row is the weighted barycenter of its columns" has no MOFA analogue. The closest object is `MOFAcellulaR::project_data` (`MOFAcellulaR/R/project_data.R:50-90`):

```r
inv_loading <- MASS::ginv(loading_matrix)
projected_factors <- test_data_mat %*% inv_loading
```

This projects new samples via the Moore–Penrose pseudoinverse `Z_new = X_new · W^+`. It is a least-squares regression onto W's column space — not a barycenter. It works on numeric features only (not categorical levels) and is post-hoc on already-trained W. It implicitly assumes the trained Z is approximately `X · W^+`, which holds only when noise precision Tau is high and prior precision Alpha is low.

### 1.3 Per-view W means no single feature cloud

With M views there are M separate loading matrices each living in `R^{D_m × K}` (`mofapy2/build_model/init_model.py:557-614`). Each W^(m) has its own ARD prior `AlphaW^(m, :)` (`Alpha_nodes.py:11-89`). Genes in view 1 and methylation features in view 2 share factor *indices* only, not a coordinate space. There is no joint feature cloud to put on a biplot.

### 1.4 Z and W are not jointly orthogonal

The Gaussian prior on Z is per-factor independent (`Z_nodes.py:158-167`); the posterior is updated via coordinate-wise VI. Orthogonality is *encouraged* by the prior but not enforced. `MOFA2/R/plot_factors.R:471` documents this directly: "the model encourages the factors to be uncorrelated, so this function usually yields a diagonal correlation matrix". `plot_factor_cor` (`plot_factors.R:487`) exists so users can verify empirically.

### 1.5 No biplot function exists

Verified by grep:

```bash
grep -rn "biplot" /workspaces/DC_hum_verse/01_modules/.ref/MOFA2/ \
                  /workspaces/DC_hum_verse/01_modules/.ref/MOFAcellulaR/ \
                  /workspaces/DC_hum_verse/01_modules/.ref/mofapy2/
# (returns empty)

ls /workspaces/DC_hum_verse/01_modules/.ref/MOFA2/man/ | grep -i biplot
# (returns empty)
```

The 25 `plot_*.Rd` man pages cover Z, W, factor correlations, variance explained, data heatmaps, enrichment, MEFISTO alignment — never a co-projection of Z and W on shared axes. **The absence is structural, not an oversight.**

If you need a biplot, you either build one yourself (see `references/supplementary-projection-on-z.md` for the L1 patch) or use a different family (FAMD/MFA via `factor-analysis-framework`).

---

## 2. Inertia decomposition — ABSENT

See `references/variance-explained.md` for the full treatment. Headline: per-(view, group, factor) R² is computed against the full Y (not residual), Z is not orthogonalised through training, factor ordering is by total R² (not eigenvalue), no `cos²` / `ctr` vocabulary exists. The R² matrix is a useful per-view diagnostic, not a budget.

---

## 3. Categorical-level-as-barycenter — ABSENT at fit-time

Categorical metadata never enters the MOFA fit. In `entry_point.set_data*` (`mofapy2/run/entry_point.py:201-388`), data is a list of numeric view matrices; sample metadata is held only as a row of strings outside the model. `set_covariates` (`entry_point.py:81-130`) is **not** a categorical hook — it is the MEFISTO smooth-covariate mechanism (continuous time/space coordinates fed to a Gaussian-process prior on Z; `mofapy2/core/nodes/Sigma_node.py`, `Z_nodes_GP.py`).

Categorical structure surfaces only **post-hoc**, in three ways, none of which are coordinates participating in the decomposition:

1. **As a colour aesthetic.** `plot_factors(model, color_by="condition")` (`MOFA2/R/plot_factors.R:290-379`) and `MOFAcellulaR::plot_sample_2D` (`MOFAcellulaR/R/plot_sample_2D.R:41-113`) consume the categorical metadata as a `ggplot` `color_by` / `shape_by` aesthetic. The level has no position; the cells coloured by it do.

2. **As an ANOVA outcome.** `MOFAcellulaR::get_associations` (`MOFAcellulaR/R/get_associations.R:49-115`) fits `aov(value ~ test_variable)` factor-by-factor against the categorical metadata and BH-adjusts p-values. This tests *whether* a factor separates conditions, not *where* the condition sits.

3. **As a median per level.** `MOFA2::summarise_factors` (`MOFA2/R/correlate_covariates.R:114-171`) computes the median of Z within each level of a discrete grouping. This is the closest object in the codebase to a "level barycenter" — and it is what a Benzécri supplementary projection would compute — but it is rendered as a tile heatmap (`correlate_covariates.R:149-156`), not as a labelled point on `plot_factors`.

The Level-1 patch turns operation (3) into a co-plot. See `references/supplementary-projection-on-z.md`.

---

## 4. Distribution-free framing — ABSENT

MOFA is fully Bayesian. The model declares priors on Z (Gaussian or spike-and-slab; `Z_nodes.py:16, 188`), on W (Gaussian or spike-and-slab; `W_nodes.py:20`), on ARD precisions (Gamma; `Alpha_nodes.py:11, 92`), and on noise precision Tau (Gamma; `Tau_nodes.py`). Inference is coordinate-ascent variational, monitoring an ELBO.

Every consequence the whitepaper enumerates lands:

| Consequence | Code location |
|---|---|
| Output depends on prior (spike-and-slab vs ARD-only) | `mofapy2/build_model/build_model.py:71-84`; `entry_point.py:987-1007` |
| Output depends on initialisation (`random` / `orthogonal` / `pca`) | `mofapy2/build_model/init_model.py:78` |
| Output depends on seed | `train_opts['seed']` at `entry_point.py:1017` |
| Convergence diagnostic required | ELBO trace at `freqELBO`, tolerance early-stop at `entry_point.py:961` |
| Sign-flip equivariance per factor | Inherent to bilinear `Y ≈ ZW^T` |
| Permutation equivariance across factors | Inherent; partially mitigated by sort-on-R² (`save_model.py:115-119`) |

The Bayesian framing is the **architecture**, not an accident.

---

## 5. Sign / order indeterminacy — what to do about it

### 5.1 Sign flips

The bilinear likelihood `Y ≈ Z · W^T` is sign-flip equivariant per factor: replacing `(Z[:, k], W[:, k])` with `(-Z[:, k], -W[:, k])` leaves the likelihood unchanged. Variational inference offers no preferred sign. Two consequences:

- **Cross-run biplots may be sign-flipped.** A factor that points "high → tumour" in run 1 may point "high → normal" in run 2.
- **The codebase does not enforce a canonical sign convention.** `save_model.py:115-119` reorders factors by total R² (fixing order, not sign).

**Recommendation:** before rendering, impose a canonical sign. The simplest convention is *largest-absolute-loading positive*:

```python
import numpy as np

def canonicalise_sign(Z, W_list):
    """Flip signs so the largest |W| entry per factor is positive."""
    for k in range(Z.shape[1]):
        max_w = max(np.max(np.abs(W[:, k])) for W in W_list)
        signs = [np.sign(W[np.argmax(np.abs(W[:, k])), k]) for W in W_list]
        # Pick the view that achieved the max
        view_with_max = np.argmax([np.max(np.abs(W[:, k])) for W in W_list])
        s = signs[view_with_max]
        if s < 0:
            Z[:, k] *= -1
            for W in W_list:
                W[:, k] *= -1
    return Z, W_list
```

Alternative conventions:
- Largest-loading **gene-symbol-first-alphabetical** positive (avoids dependence on which view wins the max).
- A fixed anchor: pick a marker gene known to be positive on the factor.

### 5.2 Order

`save_model.py:115-119` sorts factors descending by total R² across groups and views. This pins order *within* a run but not *across* runs — two re-fits can converge to permuted-then-resorted orderings if some factors swap R² rank.

**Recommendation:** match factors across runs by Tucker congruence on Z:

```python
from itertools import permutations
import numpy as np

def tucker_congruence(z_a, z_b):
    """Cosine similarity between matched factor pairs."""
    z_a = z_a / np.linalg.norm(z_a, axis=0, keepdims=True)
    z_b = z_b / np.linalg.norm(z_b, axis=0, keepdims=True)
    return np.diag(z_a.T @ z_b)
```

Tucker congruence `> 0.85` is the standard threshold for "the same factor"; `< 0.85` flags instability and is the stat reviewer's hard gate (`docs/vision/latent/synthesis.md:27`). This is **not** implemented in MOFA2 / MOFAcellulaR; the consumer must add it.

---

## 6. Open question: canonical sign at the platform layer

If a downstream platform (e.g. pathway-explorer's biplot layer) consumes MOFA output, it should either:

1. Impose a canonical sign convention before rendering and document the choice, **or**
2. Refuse to render Z below a stability gate (Tucker congruence `< 0.85` versus a reference fit).

The whitepaper's "biplot first" stance (`docs/vision/vision.md:33-50`) implies (1). The latent-lens stat reviewer's stability gate (`synthesis.md:27, 66`) implies (2). Both apply, and they are not mutually exclusive.

See `docs/vision/latent/mofa.md` §8 for the open-question list and `references/troubleshooting.md` for the "sign flips between runs" entry.
