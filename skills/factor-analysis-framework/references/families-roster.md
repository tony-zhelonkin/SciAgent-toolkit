# Factor-analysis families roster

This reference page lists each family that will eventually have a `*-framework` child skill in this super-family. For each family: name, primary codebases, languages, scope (single-cell vs pseudo-bulk), modality (single vs multi), strengths, weaknesses, and status in this skill toolkit.

Conceptual content (the four Benzécri properties, the MODEL/MAP axes, the Level-1 patch) is *not* repeated here — see `benzecri-properties.md`, `modelling-vs-map.md`, `supplementary-projection.md`.

For the cross-cutting audit table (which property each codebase preserves / ablates) see `factor-analysis-framework/SKILL.md` §"The Benzécri four-property frame".

---

## MOFA family — Bayesian multi-view, scale wins, geometry loses

**Codebases.** `mofapy2` (Python core), `MOFA2` (R bindings), `MOFAcellulaR` (R, pseudo-bulk wrapper for single-cell), `mofax` (Python loader).

**Languages.** Python (core, CLI), R (R user-facing API + pseudo-bulk wrapper).

**`.ref/` paths.**
- `/workspaces/DC_hum_verse/01_modules/.ref/mofapy2/`
- `/workspaces/DC_hum_verse/01_modules/.ref/MOFA2/`
- `/workspaces/DC_hum_verse/01_modules/.ref/MOFAcellulaR/`
- `/workspaces/DC_hum_verse/01_modules/.ref/mofax/`

**Scope.** Single-cell (MOFA+, via SVI + GPU; up to 10⁶ cells), pseudo-bulk (MOFAcellulaR; ~10²–10³ donors).

**Modality.** Multi-view by construction (one likelihood per view). Supports gaussian / bernoulli / poisson per view via `set_likelihoods` (`entry_point.py:371-382`).

**Strengths.**
- View-specific likelihoods — the only audited modern method that handles this honestly.
- ARD sparsity priors on per-view loadings (`Alpha_nodes.py:11-89`) produce readable factor loadings (~50 non-zero genes per factor instead of 30,000).
- Missing-data marginalisation (NaN mask in `Y_Node`; ELBO integrates over it; `Y_nodes.py:46-50`).
- Shared-vs-private factor decomposition — genuinely novel affordance not in FactoMineR.
- Scale: SVI + GPU reach 10⁶ cells. Largest scale ceiling in the audited stack.

**Weaknesses.**
- Strips all four Benzécri properties (see `benzecri-properties.md`).
- Z and W are posterior expectations of two latent variables linked by a likelihood `Y ≈ ZW^T`; no single SVD, no transition formula, no biplot function across the three repos.
- Per-(view, group, factor) R² does *not* sum to total R²; `max(0, .)` floor at `calculate_variance_explained.R:73-76` confirms it is not a partition.
- Fully Bayesian VI — output depends on prior, init, seed.
- No categorical-as-barycenter at fit; post-hoc numbers exist in `summarise_factors` / `correlate_factors_with_covariates` but are rendered as `geom_tile` / `corrplot`, not as biplot points.

**Status in toolkit.** **Implemented.** `mofa-framework` is the family-framework concept skill; `mofa-mofapy2`, `mofa-r`, `mofa-cellular`, `mofa-mofax` are the implementation children.

---

## FactoMineR / prince family — geometric ground truth

**Codebases.** `FactoMineR` (R, the canonical implementation), `prince` (Python port for PCA / CA / MCA / FAMD / MFA).

**Languages.** R (FactoMineR), Python (prince).

**`.ref/` paths.**
- `/workspaces/DC_hum_verse/01_modules/.ref/FactoMineR/`
- `/workspaces/DC_hum_verse/01_modules/.ref/prince/`

**Scope.** Pseudo-bulk / sample-level. FactoMineR caps at ~10⁴ × 10⁴ dense; prince's randomized SVD could reach 10⁶ rows but its FAMD path densifies on `pd.get_dummies` (`famd.py:57`) which blocks it.

**Modality.** Multi-block via MFA. Mixed continuous + categorical via FAMD. Categoricals-only via MCA. Counts / contingency tables via CA.

**Strengths.**
- All four Benzécri properties (see `benzecri-properties.md`). Ground truth for the geometric axis.
- `quali.sup`, `quanti.sup`, `ind.sup` are first-class arguments on every active function. `eta²`, `v.test` for every supplementary level.
- Single weighted triplet SVD (`svd.triplet.R:16-86`); transition formula explicit in `predict.PCA.R:14`.
- MFA's per-block `sqrt(λ_1^(b))` rescaling makes mixed-block analyses scale-invariant.

**Weaknesses.**
- Scale ceiling. ~10⁴ × 10⁴ dense in FactoMineR. prince's randomized SVD via `sklearn.utils.extmath.randomized_svd` is faster but the FAMD path densifies on one-hot expansion.
- No sparse path. No GPU. No SVI.
- Single Gaussian-on-everything contract — no per-view likelihoods.
- No ARD sparsity. Loadings on 30,000 genes are dense and unreadable without external feature ranking.
- prince's FAMD uses a different rescaling constant (`prop * 2` vs FactoMineR's `prop`) and exposes no `quali.sup` / `quanti.sup`; partial compatibility with FactoMineR.

**Status in toolkit.** **Planned** (`factomineR-framework`). The natural home for the L2 MFA-on-pseudo-bulk wedge and for full classical PCA / CA / MCA / FAMD recipes.

---

## mixOmics family — sparse multi-block via PLS

**Codebases.** `mixOmics` (R; includes `block.spls`, `block.splsda` = DIABLO).

**Language.** R.

**`.ref/` path.** `/workspaces/DC_hum_verse/01_modules/.ref/mixOmics/`

**Scope.** Sample-level (tens to low hundreds; caps in low-thousand range due to dense `base::svd`). Not single-cell.

**Modality.** Multi-block via `block.spls` and DIABLO. Sparse via soft-thresholding (`R/spca.R:262`).

**Strengths.**
- Sparse loadings via L1 / soft-thresholding — the chemometric / PLS tradition does this well.
- Per-block redundancy (`Rd(X, t_h)`) is a meaningful within-block diagnostic.
- Deterministic NIPALS / SVD — distribution-free framing retained.
- Correlation-circle plots (`plotVar`) and biplots (`R/biplot.R`) — visualisation surface is large.

**Weaknesses.**
- **Partial duality.** `R/biplot.R:77` uses a *cosmetic* rescaler `max(variates) / max(abs(loadings))` instead of `sqrt(λ_α)`; arrows fit inside the cloud but the "row equals weighted barycenter of columns" identity does not hold quantitatively.
- **Partial inertia.** `prop_expl_var` is redundancy, not the eigenvalue fraction. Docstring even warns this need not decrease with components (`R/pls.R:99-106`).
- **Fails categorical-as-barycenter.** `unmap()` produces 0/1 indicator without the `1/√p_k` chi-square rescaling; that indicator is the response matrix in PLS-DA, not a participating or supplementary variable. Centroids in `plotIndiv` are visual group means only.
- No sparse pipeline; n ≲ 10⁴.

**Status in toolkit.** **Planned** (`mixomics-framework`).

---

## DIABLO — multi-block supervised, lives inside mixOmics

**Codebase.** `mixOmics::block.splsda` (R). Audit at `docs/vision/latent/diablo.md`.

**Language.** R.

**`.ref/` path.** `/workspaces/DC_hum_verse/01_modules/.ref/diablo/` (and the parent `mixOmics/`).

**Scope.** Sample-level (≤ few thousand). Dense everywhere; `ginv()` on n×n when n<p_j (`internal_mint.block.R:740-760`).

**Modality.** Multi-block supervised — Y is the outcome, drives the axes.

**Strengths.**
- Generalised canonical correlation across blocks: `sum(design * cov2(variates.A))`.
- Per-block paired score + loading; useful for *within-block* duality (`internal_mint.block.R:596`).
- Sparse via soft-thresholding (inherited from mixOmics).
- Deterministic block-power-iteration; `init = "svd.single"`.
- Selection-commit-adjacent: `selectVar()` is closer to a commit contract than to a latent-lens scatter.

**Weaknesses.**
- **No across-block duality.** Per-block scores live in shared n-space but loadings live in block-specific p_q-spaces; no single SVD whose left and right singular vectors put all blocks' features and the samples on the same axes.
- **No λ_α partition.** AVE + per-block redundancy reported; loadings are L2-unit-normalised directions, not contributions.
- **Categorical-as-barycenter is structurally violated.** Y is the outcome — categorical drives geometry rather than receives coordinates. No `quali.sup` for unmodelled covariates (donor, batch, biopsy site); the design matrix is block-to-block coupling, not level projection.

**Status in toolkit.** **Planned** (either `diablo-framework` or absorbed into `mixomics-framework` as a sub-method). The audit suggests absorption — DIABLO's API is `block.splsda` inside the mixOmics package, and most of the contraindications (categorical-as-outcome inheritance) apply to the supervised side of mixOmics generally.

---

## scITD — Tucker + ICA on pseudo-bulk tensor

**Codebase.** `scITD` (R; tensor decomposition on donor × gene × cell-type pseudo-bulk).

**Language.** R.

**`.ref/` path.** `/workspaces/DC_hum_verse/01_modules/.ref/scITD/`

**Scope.** Pseudo-bulk × cell-type × donor tensor. O(10²–10³) donors. Sub-second on 100×500×8. n≪p in donor mode.

**Modality.** Single-omic (RNA) but multi-block-equivalent through the cell-type-as-view dimension.

**Strengths.**
- Tucker-2 + ICA on a donor × gene × cell-type tensor — genuinely novel decomposition for pseudo-bulk single-cell.
- Distribution-free: Tucker + HOSVD-ALS + ICA + varimax all SVD-based; no priors.
- `get_meta_associations` (`get_meta_associations.R:16-87`) is the post-hoc metadata regression; reusable as a barycenter source with a 20-line patch.
- A proper-biplot sibling exists in the codebase: `pca_unfolded` (`run_tucker_ica.R:318-361`), retained as a diagnostic only.

**Weaknesses.**
- **Partial duality.** Tucker-2 produces (donor scores, gene×celltype loadings), but the core tensor is fused into the loadings at `run_tucker_ica.R:137` before user-facing output; only `pca_unfolded` is a proper biplot.
- **Non-additive inertia.** `get_factor_exp_var` is rank-1 reconstruction loss, not a partition. ICA rotation makes per-factor exp_var non-orthogonal and non-additive; the current display at `plot_tucker.R:47-53` implies additivity that does not hold.
- **Fails categorical-as-barycenter.** `get_meta_associations` returns R²/p-values, not coordinates. Donors get score-vectors, cell-types get axis-labels; metadata appears as coloured columns in `plot_donor_matrix` — the whitepaper §4 "categorical as colour" pattern in literal form.

**Status in toolkit.** **Planned** (`scitd-framework`). The audit names a "20-line patch" (export per-level means as coordinates instead of regression R²) as the cheapest move bringing scITD halfway to L1.

---

## sciRED — PCA + varimax + Poisson-GLM nuisance regression

**Codebase.** `sciRED` (Python; cell-grain PCA with varimax rotation and Poisson-GLM residualisation).

**Language.** Python.

**`.ref/` path.** `/workspaces/DC_hum_verse/01_modules/.ref/sciRED/`

**Scope.** Cell-level. 10³–10⁵ cells. Dense `sklearn.decomposition.PCA`; no count-aware decomposition.

**Modality.** Single-omic. Categoricals enter as regressors / colourings / classifier targets.

**Strengths.**
- PCA core is distribution-free (sklearn SVD).
- Varimax rotation produces sparse interpretable loadings (a chemometric move).
- Poisson-GLM nuisance regression per gene (`glm.poissonGLM`) is a principled residualisation step.
- FCAT (`ensembleFCA.py:198-218`) — supervised factor-classifier scoreboard for marker-based factor interpretation.

**Weaknesses.**
- **Strips duality.** Gene scatter (`plot_factor_loading`) and cell scatter (`plot_factor_scatter`) are separate matplotlib figures with no shared rescaling. Varimax further breaks SVD uniqueness in exchange for sparse loadings.
- **Inertia partition unused.** `pca.explained_variance_ratio_` computed but never read; `factor_variance` in FIST is min-max-scaled out of recognition.
- **Fails categorical-as-barycenter.** Categoricals are nuisance regressors, supervised classifier targets, or cell colourings — never barycenters. The mean factor score per level is not computed anywhere despite being one line of pandas.
- **Distribution-free framing partly broken.** PCA core OK; Poisson GLM upstream is likelihood-based; varimax breaks SVD uniqueness deliberately; FCAT is supervised.

**Status in toolkit.** **Implemented** (skill: `scired`). Single-codebase implementation skill — sciRED has no sub-implementations, so it slots in as a direct sibling of `mofa-framework` under `factor-analysis-framework` rather than as its own family-router.

---

## Status summary

| Family | Skill name | Status |
|---|---|---|
| MOFA | `mofa-framework` (+ `mofa-mofapy2`, `mofa-r`, `mofa-cellular`, `mofa-mofax`) | **Implemented** |
| sciRED | `scired` (single-codebase implementation, no family-router) | **Implemented** |
| FactoMineR / prince | `factomineR-framework` | Planned |
| mixOmics | `mixomics-framework` | Planned |
| DIABLO | `diablo-framework` or absorbed into `mixomics-framework` | Planned |
| scITD | `scitd-framework` | Planned |

Today the super-family router points the multi-view branch to `mofa-framework` and the single-omic single-codebase branch to `scired`; the other family branches are stubs and emit "planned" notices in the decision tree (`factor-analysis-framework/SKILL.md` §"Family-selection decision tree"). The Level-1 supplementary-projection patch (see `supplementary-projection.md`) is the first concrete deliverable that works on top of *any* of these families' linear embeddings, including the implemented MOFA and sciRED ones.
