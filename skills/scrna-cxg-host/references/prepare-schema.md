# Phase A — schema preparation reference

CellxGene's schema is strict but mostly undocumented in user-friendly form. Below is the operative version distilled from the reference helpers in `01_Scripts/Python_scripts/cxg_utils.py` of the 13403-YD reference, plus the failure modes each step prevents.

## Six steps, in order

### A.1 — `convert_obsm_to_arrays(adata)`

For each key in `adata.obsm`:

- If it is a `pd.DataFrame`: `adata.uns[f"{key}_columns"] = list(df.columns)`; replace with `df.values.astype("float32")`.
- If it is a `np.ndarray`: leave it.
- Otherwise (rare): warn and skip.

The float32 cast is intentional — CXG's frontend uses float32 internally; storing float64 doubles disk and RAM with no benefit.

### A.2 — `ensure_unique_varnames(adata)`

Preference for `var_names`: `gene_name` (symbols) → `gene_id` (Ensembl) → existing index. Duplicates get `__N` suffix where N starts at 1 (`Mboat2`, `Mboat2__1`, `Mboat2__2`). Original index is preserved in `adata.var["__orig_var_index"]`.

**Both `.var` AND `.raw.var` must be reindexed.** CellxGene queries dispatch through `.raw.X` whenever `.raw` is set, and the column lookup uses `.raw.var.index`. If only `.var` is reindexed, the gene autocomplete still returns Ensembl ids — the `.raw` shadow wins. Because `.raw` is read-only, the helper rebuilds it via `adata.raw.to_adata()` → set `var_names` → `adata.raw = raw`. The same `__N` dedup and `__orig_var_index` preservation applies on each side.

### A.3 — `ensure_unique_barcode_and_index(adata, joiner="_")`

This step has a sample-id-prefix detection step. Pseudocode:

```python
adata.obs["barcode_raw"] = adata.obs["barcode"].astype(str) if "barcode" in adata.obs else adata.obs.index.astype(str)
if "sample_id" in adata.obs:
    sample_ids = adata.obs["sample_id"].astype(str)
    # Sniff first 10 barcodes to detect already-prefixed cases (avoids 'sid-sid-AAA...')
    already_prefixed = all(bc.startswith(sid + delim)
                           for bc, sid in zip(adata.obs["barcode_raw"][:10], sample_ids[:10])
                           for delim in ("-", "_", "."))
    base = adata.obs["barcode_raw"] if already_prefixed else sample_ids + joiner + adata.obs["barcode_raw"]
else:
    base = adata.obs["barcode_raw"]
adata.obs["barcode"] = make_unique(base)            # __N suffix on duplicates
adata.obs["cell_id"] = adata.obs["barcode"].astype(str)
adata.obs.index = pd.Index(adata.obs["cell_id"].values, dtype=str, name="cell_id")
```

The string-typed index is load-bearing — CellxGene's index validator infers dtype via `pd.api.types.infer_dtype(idx)` and rejects anything that is not `'string'`.

### A.4 — `realign_aligned_mappings(adata)`

After A.3 changes the obs index, any DataFrame in `obsm`/`varm` keeps its old index — and CXG's final-checks alignment validator catches that. The helper iterates `adata._obsm._data` (or `adata._obsm` if `_data` is absent in the running anndata version), and for each DataFrame value: `df.index = pd.Index(target_obs_index.values, dtype=str, name="cell_id")`. NumPy arrays in obsm are left alone after a shape check.

### A.5 — `ensure_umap(adata, n_pcs=50)`

Triple-conditional:

```python
if "X_pca" not in adata.obsm or adata.obsm["X_pca"].shape[0] != adata.n_obs:
    sc.pp.pca(adata, n_comps=min(n_pcs, max(2, adata.n_vars - 1)))
if "neighbors" not in adata.uns or not all(k in adata.obsp for k in ("distances","connectivities")):
    sc.pp.neighbors(adata)
if "X_umap" not in adata.obsm or adata.obsm["X_umap"].shape != (adata.n_obs, 2):
    sc.tl.umap(adata)
```

Skip the recompute when valid; recompute when missing or mis-shaped.

### A.6 — `final_checks(adata)`

```python
assert adata.obs_names.is_unique
assert adata.var_names.is_unique
assert "barcode" in adata.obs and adata.obs["barcode"].is_unique
assert "X_umap" in adata.obsm and adata.obsm["X_umap"].shape[1] == 2
assert pd.api.types.infer_dtype(adata.obs.index) == "string"
assert pd.api.types.infer_dtype(adata.var.index) == "string"
for k, v in adata.obsm.items():
    if hasattr(v, "index"):
        assert v.index.equals(adata.obs.index), f"obsm[{k!r}] index drift"
        assert pd.api.types.infer_dtype(v.index) == "string"
    else:
        assert v.shape[0] == adata.n_obs
```

## After Phase A

```python
adata.uns["title"] = "<project> • <stage>"
adata.write_h5ad("03_results/objects/<NN>_explore_<name>.h5ad", compression="gzip")
```

Naming convention follows `scrna-pipeline-conventions`: `08_explore_full.h5ad`, `09_explore_<celltype>.h5ad`. The two-digit prefix lets multiple datasets sort meaningfully in `ls`.

## Idempotency

All six steps are idempotent — running them twice produces the same output. This matters because `prepare()` may be called once for the full dataset and again per subset (Phase A-bis). The skill never assumes "this AnnData is fresh".

## Where the helpers live

The skill ships these as Python functions in `scripts/prepare_for_cxg.py`. Other code in the project should import from there rather than re-implementing — the `__N` suffix logic, the sample-id prefix detection, and the `realign_aligned_mappings` wrapper-aware code are all subtle and tested.
