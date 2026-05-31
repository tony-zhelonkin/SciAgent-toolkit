# MOFA architecture — Markov blanket, bilinear likelihood, ARD, VI loop

This file describes the *model* and its inference loop. Use it when a child skill needs to explain why a knob does what it does (sparsity, factor pruning, group structure, view-weighting) and the answer involves the underlying Bayesian network. Source: `docs/vision/latent/mofa.md` §1.1, §2.1, §2.4 and the `mofapy2` source tree.

---

## 1. The Markov blanket

The variational Bayesian factor model is assembled from a small set of node classes under `mofapy2/core/nodes/`. The blanket is registered in `mofapy2/build_model/build_model.py:175-178`:

```python
nodes["Y"].addMarkovBlanket(Z=nodes["Z"], W=nodes["W"], Tau=nodes["Tau"])
nodes["Z"].addMarkovBlanket(Y=nodes["Y"], W=nodes["W"], Tau=nodes["Tau"])
nodes["W"].addMarkovBlanket(Y=nodes["Y"], Z=nodes["Z"], Tau=nodes["Tau"])
nodes["Tau"].addMarkovBlanket(Y=nodes["Y"], W=nodes["W"], Z=nodes["Z"])
```

This is the entire generative model. The remaining nodes (`AlphaW`, `AlphaZ`, `Theta`) attach as priors on `W`, `Z`, and the spike-and-slab mixture proportion.

| Node | File | Role | Shape |
|---|---|---|---|
| `Y_Node` | `mofapy2/core/nodes/Y_nodes.py:14-115` | Observed data per view; held as a *constant* node with NaN-mask for missingness | One per view, `D_m × N` |
| `Z_Node` | `mofapy2/core/nodes/Z_nodes.py:16-185` | Latent factor matrix; **single posterior shared across all views** | `N × K` |
| `SZ_Node` | `mofapy2/core/nodes/Z_nodes.py:188-468` | Spike-and-slab variant of Z (Bernoulli × Gaussian) | `N × K` |
| `W_Node` | `mofapy2/core/nodes/W_nodes.py:20-130` | Per-view loading matrix; **one independent W per view** | `D_m × K`, M copies |
| `SW_Node` | `mofapy2/core/nodes/W_nodes.py` (spike-and-slab) | Sparse W via Bernoulli × Gaussian | `D_m × K` |
| `AlphaW_Node` | `mofapy2/core/nodes/Alpha_nodes.py:11-89` | Gamma-distributed precision on W (ARD shrinkage) | One scalar per `(view, factor)` |
| `AlphaZ_Node` | `mofapy2/core/nodes/Alpha_nodes.py:92-` | Gamma precision on Z per sample-group | One scalar per `(group, factor)` |
| `Tau_Node` | `mofapy2/core/nodes/Tau_nodes.py` | Per-view noise precision (Gaussian likelihood only) | One scalar per `(view, feature)` |
| `Theta_Node` | `mofapy2/core/nodes/Theta_nodes.py` | Spike-and-slab mixture proportion | Per `(view, factor)` |

### 1.1 Why "one independent W per view" matters

`mofapy2/build_model/init_model.py:557-628` loops `for m in range(M): W_list[m] = ...`, instantiating M separate `W_Node` (or `SW_Node`) objects. Each has its own ARD precision vector `AlphaW[m, :]`. This is the structural reason MOFA can mix per-view likelihoods — each view has its own loading subspace, its own scale, and its own ARD prior. It is also the structural reason there is **no shared feature cloud**: genes in view 1 and pathway scores in view 2 do not occupy a common coordinate system; they share factor *indices* only.

### 1.2 Spike-and-slab variants

Two independent toggles control sparsity:

- `spikeslab_factors=True` (`build_model/build_model.py:92-100`) — sparse Z (e.g. each cell uses only a few factors).
- `spikeslab_weights=True` (`build_model/build_model.py:113-120`) — sparse W (e.g. each factor uses only a few genes).

The default for the multi-omics use case is `spikeslab_weights=True, spikeslab_factors=False`. The spike-and-slab adds a per-element Bernoulli mask `S` on top of the Gaussian slab, so the effective W is `S ⊙ W̃`. Implementation: `SZ_Node` at `Z_nodes.py:188-468` and `SW_Node` in `W_nodes.py`.

---

## 2. The bilinear likelihood

The single equation tying everything together is the per-view Gaussian-form ELBO term in `Y_Node.calculateELBO` (`mofapy2/core/nodes/Y_nodes.py:74-115`). Quoting the squared-error expansion at `Y_nodes.py:98-104`:

```python
tmp = Y**2 + ZZ.dot(WW) - Z2.dot(W2) + (ZW)**2 - 2*(ZW)*Y
elbo += 0.5 * (Tau_lnE * foo).sum() - (Tau_E * tmp).sum()
```

This encodes `E[(Y − ZW^T)²]` under the variational posterior `Q(Z)Q(W)`. For non-Gaussian views (`bernoulli`, `poisson`) the same bilinear coupling is wrapped in a Bohning / Jaakkola bound — see `mofapy2/core/distributions/` and `mofapy2/core/nodes/Tau_nodes.py`. Either way the *coupling* is bilinear: `Y_m ≈ Z · W_m^T`.

**Consequence: this is not an SVD.** An SVD `M = UΣV^T` is a single deterministic decomposition with a shared Σ that distributes scale between `U` and `V`. MOFA has no Σ. Z and W are two separate posterior expectations, computed by alternating coordinate-ascent updates. The closest thing to "singular values" is the per-view ARD precision `AlphaW[m, :]`, which controls how aggressively each factor is shrunk *in that view* — but it is not shared with Z.

---

## 3. ARD shrinkage drives factor pruning

`AlphaW_Node` (`mofapy2/core/nodes/Alpha_nodes.py:11-89`) is a Gamma-distributed precision per `(view, factor)`. Its variational update at `Alpha_nodes.py:60`:

```python
Qb = Pb + 0.5 * WW.sum(axis=0)
```

If a factor has small `W[m, :, k]` mass in view m, `Qb` stays close to the prior `Pb`, the posterior precision `Alpha[m, k]` stays large, and the prior on `W[m, :, k]` becomes tight at zero — pruning the factor out of that view. If `Alpha[m, k]` is large in *every* view, the factor is globally pruned.

This is the mechanism behind:

- **Sparse, view-decomposed factors.** A factor active in views 1 and 3 but pruned in view 2 reads as "shared by views 1 and 3, absent in 2".
- **`removeInactiveFactors`** (`mofapy2/core/BayesNet.py:178-220`) — at every iteration, factors with R² below `drop_factor_threshold` across all views are dropped. This dynamically shrinks K.
- **The factor-ordering choice** (`mofapy2/build_model/save_model.py:115-119`) — survivors are reordered by total R² descending, mimicking PCA's eigenvalue order. This is convenient but not the same as PCA: in PCA the ordering is a property of the SVD; here it is a property of the converged posterior, which depends on seed, init, and prior.

`AlphaZ_Node` (`mofapy2/core/nodes/Alpha_nodes.py:92-`) plays the analogous role for sample groups: if a factor has small Z mass in group g, `AlphaZ[g, k]` shrinks Z's prior tight at zero in that group. This is what makes sample-group structure (donor, batch) act as a *prior* rather than a hard partition.

---

## 4. The variational inference loop

The driver is `mofapy2/core/BayesNet.iterate` (called from `entry_point.py:run`). Each iteration:

1. Update `Z` (`Z_nodes.py:49-150`) — loop over factors `for k in range(K)` at line 133, compute Gaussian posterior conditional on the current W:
   ```python
   Qvar[:, k] = 1.0 / (Alpha[:, k] + foo[:, k])
   Qmean[:, k] = Qvar[:, k] * (bar + Alpha[:, k] * Mu[:, k])         # Z_nodes.py:146-147
   ```
   where `foo` sums `tau[m] · W[m]["E2"]` across views (`Z_nodes.py:126`) and `bar` aggregates `tau[m] · Y[m] · W[m]["E"]` minus cross-factor contributions (`Z_nodes.py:127-129, 140-142`).
2. Update `W` (`W_nodes.py:72-95`) — loop over `k`, compute Gaussian posterior conditional on the current Z. With stochastic mode, the natural-gradient interpolation `ro * new + (1−ro) * old` runs here (`W_nodes.py:87-95`).
3. Update `AlphaW`, `AlphaZ` (`Alpha_nodes.py`), `Tau` (`Tau_nodes.py`), and (if spike-and-slab) `Theta` (`Theta_nodes.py`).
4. Every `freqELBO` iterations (`entry_point.py:895`), compute the ELBO and check convergence against `tolerance` (`entry_point.py:961`). Early stop when `Δ ELBO < tolerance` for `convergence_mode`-dependent number of consecutive checks.

**Coordinate-ascent VI, not gradient descent.** Each node has a closed-form variational update given its Markov blanket. The ELBO is guaranteed non-decreasing **only when the variational family is fully conjugate and updates are exact** — for the Bohning-bounded non-Gaussian views, monotonicity is approximate. See `references/troubleshooting.md` for the "ELBO non-monotone" entry.

---

## 5. View weighting

The optional `weight_views=True` flag (`entry_point.py`, surfacing in `Z_nodes.py:114-118`):

```python
total_w = sum(Y[m].shape[1] for m in range(M))
weights[m] = total_w / (M * Y[m].shape[1])
```

This rescales each view's contribution to Z's update by the **inverse of its feature count**. Without it, a 30,000-gene RNA view dominates Z over a 200-feature methylation view simply by mass. With it, every view contributes the same total weight to Z's posterior. The weighting enters the ELBO terms in `Z_nodes.py:126, 129, 142`, not the input data — the data itself is not rescaled.

Spiritually analogous to MFA's `1/λ_1^(b)` block-balancing (Escofier–Pagès), but the rescaling factor is feature count rather than a top eigenvalue, and it enters the loss weighting rather than the metric.

---

## 6. What this architecture is, and is not

**Is:**
- A Bayesian generative model with per-view likelihoods.
- A coordinate-ascent VI fit with ELBO monitoring.
- A sparse, factor-pruned, view-decomposed latent space.
- Approximately linear and approximately PCA-shaped in Z (Gaussian prior, linear bilinear coupling).

**Is not:**
- An SVD. There is no shared Σ. Z and W are not the two halves of a unique decomposition.
- A Benzécri-school biplot generator. See `references/geometric-caveats.md`.
- A method with a unique answer. The posterior depends on prior, init (`init_model.py:78`: `random` / `orthogonal` / `pca`), and seed.

Read these alongside `references/variance-explained.md` (why R² is not a budget) and `references/geometric-caveats.md` (why no biplot).
