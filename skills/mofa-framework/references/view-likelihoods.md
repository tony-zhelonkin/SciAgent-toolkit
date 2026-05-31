# Per-view likelihood selection

The single most important per-view decision in MOFA. Source: `docs/vision/latent/mofa.md` §1.1 plus the `mofapy2/build_model/utils.py` validator.

MOFA accepts three likelihoods, set per view:

| Likelihood | Data shape | Typical inputs | Internal handling |
|---|---|---|---|
| `gaussian` | Continuous, real-valued | log-CPM, VST, log-normalised expression, beta-value methylation, pathway scores | Closed-form ELBO via `Tau_Node` (`mofapy2/core/nodes/Tau_nodes.py`) |
| `bernoulli` | 0 / 1 | Binarised ATAC peaks, methylation indicator (>50% beta) | Jaakkola lower bound |
| `poisson` | Non-negative integer | Raw UMI counts, raw fragment counts, untransformed read counts | Bohning lower bound |

The validator at `mofapy2/run/entry_point.py:371-382` enforces:

```python
assert set(likelihoods).issubset({"gaussian", "bernoulli", "poisson"})
```

There is no `nb` (negative binomial) and no `zinb`. If you need NB / ZINB, use `scvi-tools`; if you need beta, use the Gaussian likelihood on logit-transformed values and document the transform.

---

## When to use which

### `gaussian`

- **Use for:** log-normalised RNA (log-CPM, log1p-CPM), VST, log2(TPM+1), microarray, ATAC peak scores after log-normalisation, pathway scores, mass cytometry.
- **Per-group centring + per-view scaling** is applied in `process_data` (`mofapy2/build_model/utils.py:55-96`) for Gaussian views only. Centring is per-group along samples (`data[m][filt, :] -= np.nanmean(data[m][filt, :], axis=0)` at `utils.py:84`); scaling is **one scalar per view** (`data[m] /= np.nanstd(data[m])` at `utils.py:88`) — not per-feature. An optional per-group scale step lives at `utils.py:91-94`. This is the only "metric" operation in the pipeline.
- **Do not log-twice.** Calling `set_data_matrix(data=log_cpm, likelihoods=["gaussian"])` is correct. Calling `set_data_matrix(data=cpm, likelihoods=["gaussian"])` is silently wrong — the model will fit, but the linearity assumption is broken.
- **Negative values are fine** (log-fold-changes, methylation deviations).

### `bernoulli`

- **Use for:** binarised ATAC peaks (peak called / not called), methylation indicators, protein-presence flags, SNV genotype binarised.
- **Inputs must be 0 / 1 only.** A float column with values in `[0, 1]` is *not* Bernoulli — that is a continuous beta-distributed view; use Gaussian on its logit transform.
- The likelihood is handled via a Jaakkola variational lower bound. There is no `Tau_Node` for Bernoulli views (no noise precision needed — Bernoulli has fixed variance).
- **Per-feature variance varies by base rate.** Low-rate features (mostly zeros) contribute less to Z. Consider feature filtering by base rate to avoid hundreds of near-constant 0-columns.

### `poisson`

- **Use for:** raw integer counts where the dispersion is approximately Poisson (rare in scRNA — counts are typically overdispersed).
- **Inputs must be non-negative integers.** No `int` dtype check is enforced at the entry point, but the math assumes Poisson.
- For overdispersed counts (most scRNA), prefer Gaussian on log-normalised data. MOFA does **not** provide a negative binomial likelihood. This is one of the trade-offs versus scvi-tools — see `references/interop-matrix.md`.
- The likelihood is handled via a Bohning lower bound.

---

## The auto-detector

`mofapy2/build_model/utils.py:99-115` provides `guess_likelihoods`:

```python
def guess_likelihoods(data):
    lik = []
    for m in range(len(data)):
        if np.all((data[m] == 0) | (data[m] == 1) | np.isnan(data[m])):
            lik.append("bernoulli")
        elif np.all((data[m] == data[m].astype(int)) | np.isnan(data[m])):
            lik.append("poisson")
        else:
            lik.append("gaussian")
    return lik
```

Rules:
- All values in `{0, 1, NaN}` → bernoulli.
- All values integer (or NaN) → poisson.
- Otherwise → gaussian.

**Trap:** raw integer counts that you *intended* to log-normalise will silently be classified as `poisson`. The auto-detector cannot tell apart "I want Poisson on raw counts" and "I forgot to log-normalise". Always pass `likelihoods=` explicitly when the data is integer-valued.

---

## Per-view centring and scaling

`process_data` (`mofapy2/build_model/utils.py:55-96`) handles only the Gaussian views; Bernoulli and Poisson are passed through unchanged (the math assumes raw values). The transformations are driven by three `data_opts` flags:

```python
# build_model/utils.py:78-94 — Gaussian branch
if likelihoods[m] in ["gaussian"]:
    if data_opts["center_groups"]:
        for g in data_opts["groups_names"]:
            filt = [gp == g for gp in samples_groups]
            data[m][filt, :] -= np.nanmean(data[m][filt, :], axis=0)
    if data_opts["scale_views"]:
        data[m] /= np.nanstd(data[m])
    if data_opts["scale_groups"]:
        for g in data_opts["groups_names"]:
            filt = [gp == g for gp in samples_groups]
            data[m][filt, :] /= np.nanstd(data[m][filt, :])
```

Two points worth flagging:

- **Scaling is per-view, not per-feature.** `data[m] /= np.nanstd(data[m])` divides the whole view's matrix by a single scalar standard deviation. Features within a view retain their relative variance ratios.
- **Centring is per-group along samples.** When you have only one group, this is just feature-mean centring with `axis=0`.

**Recommendation:** for Gaussian views, do feature centring outside MOFA (so you control it explicitly) and pass `center_groups=False, scale_views=True` to MOFA. This makes the preprocessing auditable in your own pipeline.

---

## `weight_views` — when big views dominate

If RNA has 30,000 genes and methylation has 200 features, raw ELBO weighting gives RNA a 150× advantage in driving Z's posterior. The `weight_views=True` flag rescales each view's *loss contribution* by the inverse feature count:

```python
weights[m] = total_w / (M * Y[m].shape[1])      # Z_nodes.py:114-118
```

Implementation: the weight enters `Z_nodes.py:126, 129, 142` as a multiplier on the `tau[m] · Y[m] · W[m]` term in the Z update. The input data is **not** rescaled — the loss weighting is.

When to enable:
- Multi-view with very different feature counts (RNA + methylation, RNA + ATAC peaks).
- You care about a balanced contribution per view.

When to leave off:
- All views are comparably sized.
- One view is intentionally the "primary" view (e.g. RNA + a small panel of pathway scores).

---

## The big incompatibility this solves

MOFA's three-likelihood contract is the principal advantage over FAMD-style methods. A single FAMD / MFA SVD must Gaussian-metricise everything: continuous variables get standardized-Euclidean, categoricals get chi-squared via indicator coding. Real multi-omics breaks this:

- Continuous RNA log-counts: Gaussian-Euclidean works.
- Binary ATAC peaks: Bernoulli or chi-squared on indicator coding; FAMD lumps these as "categorical".
- Raw integer counts: Poisson; FAMD requires log-transformation first, which assumes Gaussian noise after the transform.

By dispatching one likelihood per view, MOFA avoids forcing every view through the same metric. This is the **modelling-fitness** gain that the geometric audit (`docs/vision/latent/mofa.md` §2.4 — "Distribution-free framing absent") simultaneously identifies as the **map-fitness** loss. See `references/geometric-caveats.md` for the other side of that trade.
