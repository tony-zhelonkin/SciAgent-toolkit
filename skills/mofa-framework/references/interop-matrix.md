# Implementation interoperability matrix

The four MOFA implementations share a single on-disk format — the HDF5 written by `mofapy2/build_model/save_model.py`. Every other tool reads from it. This file documents the supported transitions.

## Transitions

| From | To | How |
|---|---|---|
| `mofapy2` (Python) | `MOFA2` (R) | `MOFA2::load_model("model.hdf5")` |
| `MOFA2` (R) | `mofax` (Python) | Save HDF5 from R via `run_mofa(..., outfile = "model.hdf5")`, then `mfx.mofa_model("model.hdf5")` |
| `mofapy2` (Python) | `mofax` (Python) | Native: `mfx.mofa_model("model.hdf5")` |
| `MOFAcellulaR` | `MOFA2` | Identity. `MOFAcellulaR::run_MOFA_cell` returns the standard `MOFA` R object; downstream is identical to `MOFA2`. |
| `mofax` | `mofapy2` | **Not supported.** `mofax` is read-only on the HDF5 file; it cannot re-train or modify the model. |
| `mofax` | `MOFA2` | Indirect: save the HDF5 path, load it in R via `MOFA2::load_model`. |
| Trained `MOFA2` model | New samples (project_data) | `MOFAcellulaR::project_data(model, new_data)` — `Z_new = X_new · ginv(W)` (`MOFAcellulaR/R/project_data.R:50-90`). |
| `MOFA2` HDF5 | `MuData` / `AnnData` | Load via `mofax`, attach `Z` to `adata.obsm`, `W[m]` per-modality. |

## The HDF5 contract

The lingua franca written by `save_model.py:save_model`:

```
model.hdf5
├── data/                       (per group, per view: D_m × N)
├── expectations/
│   ├── Z/<group>               N × K
│   ├── W/<view>                D_m × K
│   ├── AlphaW/<view>           K-vector (ARD precisions)
│   ├── AlphaZ/<group>          K-vector
│   └── Tau/<view>              D_m-vector
├── variance_explained/
│   ├── r2_total/<group>        per-view scalar
│   └── r2_per_factor/<group>   per-view K-vector
├── samples_metadata/<group>    (optional)
├── features_metadata/<view>    (optional)
├── intercepts/<view>           per-feature scalar
├── training_opts                (seed, convergence_mode, etc.)
└── model_opts                   (likelihoods, num_factors, spikeslab flags)
```

Anything written by mofapy2 is readable by MOFA2 / mofax. The R wrapper writes the same shape via reticulate (MOFA2 R doesn't reimplement HDF5 I/O).

## When to chain

- **mofapy2 → mofax:** the canonical Python lifecycle. Fit with mofapy2, save HDF5, hand to mofax for downstream plotting and `get_factors` / `get_weights` / `get_r2`. mofax is a thin wrapper around `h5py` reads.
- **mofapy2 → MOFA2 R:** when the user team wants Python fitting (e.g. GPU on a Linux box) but R-side plotting for `plot_factors`, `plot_weights`, `correlate_factors_with_covariates`. Round-trip is lossless.
- **MOFAcellulaR fit → MOFA2 plotting → mofax for Python downstream:** entirely legitimate. MOFAcellulaR's model object *is* a `MOFA` object, and the HDF5 it writes is the standard format.
- **MOFA2 → MOFAcellulaR downstream wrappers:** if you fit with `MOFA2` but want MOFAcellulaR's `get_associations` / `plot_MOFA_hmap` / `plot_sample_2D`, call them on the `MOFA` object directly — they consume it as-is.

## Key rules that make interop work

1. **Versions must match.** `MOFA2` (R) calls `mofapy2` via `reticulate` / `basilisk`. The basilisk-pinned mofapy2 version (see `MOFA2/R/basilisk.R`) must match the Python `mofapy2` version if both are used directly. Mismatched versions can read older HDF5s but may emit warnings or fail on newer fields.
2. **Sample names round-trip as strings.** MOFA preserves `samples_metadata` columns; check that R factor columns are not coerced into numeric on save.
3. **Likelihoods are pinned at fit time.** They are stored in `model_opts/likelihoods`; downstream tools read them but cannot change them. Re-fit if you need a different likelihood.
4. **Sign / order indeterminacy is per-fit.** A re-fit (same HDF5 path or not) can produce a sign-flipped or reordered Z; see `references/geometric-caveats.md` §5. Do not assume Factor 1 in two HDF5s is the same factor.
5. **MOFAcellulaR's grain is donor-level.** When chaining MOFAcellulaR with downstream Python tools, the `Z` rows are donors, not cells. Don't `obs`-merge with a cell-level `AnnData`.

## See also

- `references/mofa-architecture.md` — what is in the model the HDF5 stores.
- `references/troubleshooting.md` — basilisk env mismatch, HDF5 read errors.
