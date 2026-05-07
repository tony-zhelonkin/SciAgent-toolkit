# Consolidate programs onto the full barcode space

The transfer step (Stage 3) takes per-variant `usage` matrices and joins them into a single AnnData's `obs`. Cells absent from a variant get NaN. The output is `13_all_programs.h5ad`, ready for the cross-source correlation matrix in Stage 4.

## Per-variant usage matrices

Each variant's `<variant>_results.pkl` contains:

```python
{
    "usage":          pd.DataFrame,    # cells × K, values in [0, 1] simplex per row
    "spectra_scores": pd.DataFrame,    # genes × K, signed scores
    "spectra_tpm":    pd.DataFrame,    # genes × K, TPM-normalised
    "top_genes":      pd.DataFrame,    # rank × K, gene names sorted by spectra_score
    "k":              int,
}
```

`usage.index` is cell barcodes (matching the AnnData's `obs_names` after Phase A in `scrna-cxg-host` if that ran, or the raw barcode otherwise). Stage 3 joins on `usage.index ∩ adata.obs_names`.

## The transfer call

```python
def transfer(base_h5ad, sources: dict[str, Path], out_h5ad, correlation_csv):
    adata = sc.read_h5ad(base_h5ad)
    for variant_name, pkl_path in sources.items():
        with open(pkl_path, "rb") as f:
            results = pickle.load(f)
        usage = results["usage"].copy()
        # Rename columns: cNMF_<variant>_P<i>
        usage.columns = [f"cNMF_{variant_name}_P{i+1}" for i in range(usage.shape[1])]
        # Initialise NaN columns
        for col in usage.columns:
            adata.obs[col] = np.nan
        # Join on intersection of cells
        common = adata.obs_names.intersection(usage.index)
        adata.obs.loc[common, list(usage.columns)] = usage.loc[common].values
    # Cross-source correlation
    program_cols = [c for c in adata.obs.columns if c.startswith("cNMF_")]
    prog_data = adata.obs[program_cols].dropna(how="all")
    corr_matrix = prog_data.corr()
    corr_matrix.to_csv(correlation_csv)
    adata.write_h5ad(out_h5ad)
    return adata, corr_matrix
```

## Why per-variant prefixes matter for downstream merge

In Stage 5, programs from different variants are clustered by correlation. The merge logic uses prefix-matching to record `(program → source variant)`. The variant prefixes must be:

- Unique per variant
- Length-distinguishable (`cNMF_full` ≠ `cNMF_fullQC`)
- Match the merge step's `SOURCE_PRIORITY` list (sorted by length descending; the prefix-matching bug in 13403-YD's reference)

The shipped `transfer_programs.py` and `merge_programs.py` agree on the prefix scheme; do not edit one without the other.

## Coverage diagnostics

For each variant, print:

```python
n_cells = adata.obs[f"cNMF_{variant_name}_P1"].notna().sum()
print(f"  {variant_name}: {n_cells:,} / {adata.n_obs:,} cells ({n_cells/adata.n_obs*100:.1f}%) scored")
```

Expected: variants from full-dataset runs cover ~100% (modulo cells dropped during cNMF prepare). Subset variants cover ~subset-fraction-of-total. If a variant's coverage is unexpectedly low (say, 50% of expected subset size), the cell-id index probably drifted between cNMF-input and the AnnData — debug before proceeding.

## Cross-source correlation rules

`obs[program_cols].corr()` is Pearson. Correlations are computed over cells with non-NaN values for both columns; pairwise.

- **Same-variant programs** (e.g., `cNMF_full_P1` vs `cNMF_full_P2`) should have low correlation (near-zero). cNMF's consensus produces programs with low pairwise correlation by design.
- **Cross-variant same-biology programs** (e.g., `cNMF_full_P3` vs `cNMF_fullQC_P5` both representing the cycling program) should have high correlation (`r > 0.7`).
- **All-NaN row** means no cell is shared between that program's variant and any other variant's cells — a wiring bug.

Save the correlation matrix to `03_results/tables/all_programs_correlation.csv` for Stage 5.

## Pitfall — base AnnData mismatch

The `base_h5ad` should be the *same* AnnData that cNMF was prepared from (post-`10_scored.h5ad` in the reference). If you transfer onto a different version (e.g., the original `00_raw.h5ad` before annotation), the cell-id index does not match and most variants will show 0% coverage.

Verify with: `adata.obs_names[:5]` vs the `usage.index[:5]` from any pkl — they should be identical strings.

## Output structure

```
03_results/
├── checkpoints/
│   └── 13_all_programs.h5ad           # base AnnData + cNMF_<variant>_P<i> obs columns
└── tables/
    └── all_programs_correlation.csv   # program × program Pearson correlation
```

The `13_` prefix follows `scrna-pipeline-conventions`; expect downstream stages to read from these paths.
