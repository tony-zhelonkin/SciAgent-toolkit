# Good MODEL vs good MAP — the two-axis frame

This reference page distinguishes the two orthogonal evaluation axes for factor-analysis methods in single-cell biology. They look like one axis from inside any one paper, but they are clearly two from outside, and the apparent contradiction between the ml-reviewer ranking (MOFA+ > NMF > sPLS > FAMD > raw CCA) and the geometric audit (FactoMineR ≫ MOFA) dissolves once the axes are named separately.

Content lifted with light editing from `docs/vision/latent/post-synthesis-discussion.md` §1–§2 and the synthesis reconciliation in `docs/vision/latent/synthesis.md` §4. Both files are the source of truth; quote them when communicating with reviewers.

---

## 1. What "modelling fitness" actually measures

"Modelling fitness" is an evaluation criterion the synthesis and ml.md used without defining. It asks: **how well does this method's generative model match the data's statistical structure?**

The criteria packed into that one phrase, applied to single-cell data:

- **Right likelihood per view.** Negative-binomial for UMI counts. Beta for methylation. Bernoulli for ATAC peaks. *Not* Gaussian-on-everything.
- **Missing-data handling.** Marginalises NaN entries out of the likelihood, vs requires complete-case input.
- **Sparsity.** Produces loadings with few nonzeros, so that a factor on 30,000 genes is readable instead of a flat vector of small numbers.
- **Scale.** Reaches 10⁶ cells without re-engineering the linear algebra (SVI, GPU, sparse path).
- **Shared-vs-private factor decomposition.** Per-view ARD priors that let some factors be active in one view and silent in another.
- **Posterior uncertainty.** Returns distributions, not point estimates; supports Bayesian DE / credible intervals.

The ml ranking **MOFA+ > NMF > sPLS > FAMD > raw CCA** is on this axis. It is correct on this axis. MOFA+ is a Bayesian group factor model with per-view likelihoods, ARD sparsity priors, stochastic variational inference, missing-data marginalisation, and shared/private factor decomposition. FAMD scores low on every one of those criteria: deterministic SVD on a Gaussian-metricized matrix, no likelihood, no sparsity mechanism, dense linear algebra, complete-case input only.

The ranking is asking: **how good a statistical model of the data-generating process is this method?**

That is a fundamentally different question from: how good a *map* of the data table is this method?

---

## 2. What "geometric interpretability" actually measures

The four Benzécri properties, defined operationally in `benzecri-properties.md`:

1. **Duality.** A single SVD produces row and column coordinates in the same Euclidean space, sharing axis scaling `sqrt(λ_α)`. The transition formula holds quantitatively.
2. **Inertia decomposition.** `Σ_α λ_α = TotalInertia`; per-variable `ctr` and `cos²` summing to 100% per axis.
3. **Categorical-as-barycenter.** Categorical levels receive coordinates on the same axes as the rows: `bary_{ℓ,α} = mean(F[labels == ℓ, α])`. `v.test`, `eta²`, supplementary projection of unfitted categoricals.
4. **Distribution-free framing.** No prior, no likelihood, no ELBO. One `base::svd()` call. Reproducible up to seed.

The whitepaper / cross-codebase audit ranks on this axis:

| Codebase | Geometric ranking |
|---|---|
| **FactoMineR / prince** | preserves all four |
| **DIABLO ≈ mixOmics** | preserves 1.5 of 4 (partial duality, partial inertia, fails categorical, distribution-free retained) |
| **scITD** | preserves 2 of 4 (partial duality, distribution-free; fails inertia partition under ICA rotation, fails categorical) |
| **MOFA ≈ sciRED** | strips three or four (MOFA fails all four; sciRED keeps PCA distribution-free but ablates duality, partition, categorical via Poisson-GLM + varimax) |

On *this* axis MOFA loses, by the same architectural choices that won it the modelling axis. The Bayesian bilinear-likelihood framing (Z and W as posterior expectations of two latent variables) makes a single SVD impossible; ARD on per-view AlphaW makes the inertia partition non-additive; full Bayesian VI breaks distribution-free framing by definition; the lack of a biplot function across MOFA2 / mofapy2 / MOFAcellulaR makes categorical-as-barycenter unreachable in standard usage.

The whitepaper claim is correct on *this* axis. It is also true that MOFA wins the modelling axis. Both statements are true. They are not in conflict.

---

## 3. The orthogonality, made explicit

Synthesis §4 puts it directly:

> "The ML ranking is about *modelling fitness*: which method handles n≪p, view-specific likelihoods, gene/pathway collinearity, and shared-vs-specific factor decomposition? On those axes MOFA+ wins — and the whitepaper agrees. The whitepaper audit is about *geometric interpretability*: does the method produce a single Euclidean coordinate system in which rows, columns, and categorical levels are all points with metric relationships? On that axis MOFA+ loses, by the same architectural choices that won it the modelling axes."

The two axes are orthogonal. A method can score high on one and low on the other:

|  | High modelling fitness | Low modelling fitness |
|---|---|---|
| **High geometric interpretability** | (the unicorn — does not exist in the audited stack) | FactoMineR/prince (PCA, FAMD, MFA, MCA) |
| **Low geometric interpretability** | MOFA+ | sciRED (varimax PCA), DIABLO (within-block) |

The Benzécri stance — "Le modèle doit suivre les données, et non l'inverse." — is precisely the rejection of "modelling fitness" as an evaluation criterion. Made operational:

- Do NOT posit a generative model. Don't choose a likelihood, don't fix a prior.
- Render the data as a cloud of points in a metric chosen for the data type — chi-squared for counts/categoricals, standardised-Euclidean for continuous, block-balanced for mixed blocks.
- Do SVD. Look. Read the geometry.
- If you infer anything, infer it from what the geometry shows, not from what your model assumed.

So the ml-ranking MOFA+ > FAMD is an Anglo-American ranking. From a strict Benzécri stance the ranking flips: FAMD wins because it commits to nothing about how the data was generated, and you read what is there; MOFA+ loses because it commits to a generative story you cannot inspect geometrically.

Two genuinely different epistemic stances. Neither is wrong; they answer different questions:

> "Is this method a good MODEL of my data? = does it correctly represent the data-generating process? — Anglo-American axis. The field's habit since 2010.
>
> Is this method a good MAP of my data? = does it preserve the data's structure in a readable geometry? — French school axis. The field's blind spot."
> — `post-synthesis-discussion.md §2`

---

## 4. You do not have to choose: the Level-1 patch

The synthesis's contribution — operationalised by the Level-1 supplementary-projection patch — is that you don't have to choose between the two axes for a given project.

Take a modelling-fit embedding (MOFA's Z, scVI's mean encoder output, GLM-PCA, PCA from scanpy) as the *scaffold*. Then put a French-school *overlay* on top of it: continuous variables as arrows at `(cor(v, F_1), cor(v, F_2))`, categorical levels as labelled points at `colMeans(F[labels == ℓ, :])`. The scaffold answers the modelling question; the overlay answers the map question. Both on one plot.

The math only needs the scaffold to be linear; it does not need the scaffold to be "right" in the French-school sense. The linearity assumption holds for PCA, MOFA's Z (Gaussian prior, linear `Y ≈ ZW^T`), GLM-PCA, prince's PCA / FAMD, FAMD/MFA scores. It fails for nonlinear encoders (scVI's full encoder; UMAP); for scVI use the mean encoder (linear in a neighborhood, per whitepaper §6.1); for UMAP do not project at all.

This is what makes L1 a one-afternoon patch with outsized leverage. Formulas, receipts, and open questions are in `supplementary-projection.md`.

---

## 5. How to use this in practice

Every time you read "modelling fitness", translate it as **good MODEL**. Every time you read "geometric interpretability", translate it as **good MAP**. Ask which question matters for the decision in front of you:

- **Predicting the next dataset, fitting under right likelihoods, scaling to 10⁶ cells:** the MODEL question. Route to `mofa-framework` or `scvi-framework`.
- **Understanding what is in this dataset, reading metadata onto the same axes as the samples, wet-lab handoff:** the MAP question. Route to `factomineR-framework` (planned), or apply the L1 patch on top of whichever MODEL fitter you have.
- **Both:** fit on the MODEL axis, read on the MAP axis via the L1 overlay. This is the default recommendation for a multi-omic single-cell project that wants both honesty and interpretability.

The wet-lab wedge in the pathway-explorer vision is unambiguously a MAP question. Donor / condition / treatment / response need to be points on the same axes the cells already live on, not regression coefficients in a separate table. The MODEL axis is necessary for honest fitting upstream; the MAP axis is necessary for the figure that lands on the bench.
