# Metadata schemas — observed shapes across reference projects

Real-world projects do not have a stable metadata schema. The skill discovers metadata files at conventional paths and asks the user to pick a join key (Decision Pause 2). This document records the shapes observed across audited reference projects so the agent can suggest sensible candidates without locking the user in.

## Observed schemas

### <ref-scrna> (T-cell polarisation / temperature)

File: `Counts/manifest.csv` (CSV, 5 columns).

| Column      | Example values                              | Role                |
|-------------|---------------------------------------------|---------------------|
| `sample_id` | `<sample-01>`, `<sample-02>`, …            | join key            |
| `h5`        | absolute path to `*.h5`                     | source-of-truth path |
| `celltype`  | `Th1`, `Th17`                               | biological label    |
| `temp`      | `33`, `37`, `39`                            | experimental factor |
| `replicate` | `R1`, `R2`, `R3`                            | technical replicate |

Natural join key: `sample_id`.

### <ref-multi> (mouse spleen × aging × treatment)

File: `00_data/raw/SamplesMetadata.txt` (TSV, 10 columns).

| Column            | Example values                                   | Role                     |
|-------------------|--------------------------------------------------|--------------------------|
| `Project`         | `<ref-multi>`                                       | project tag (constant)   |
| `Organ`           | `Spleen`                                         | tissue (constant)        |
| `PoolMultiSample` | `<pool-1>`, `<pool-2>`, …                 | pool id (matches `pool_id`) |
| `Sample`          | `<sample-1>`, `<sample-2>`, …                   | join key (matches `sample_id`) |
| `Sex`             | `M`, `F`                                         | covariate                |
| `Age_months`      | `5`, `27`                                        | continuous covariate     |
| `MouseID`         | `1000`, `977`, …                                 | biological replicate id  |
| `Group`           | `Young`, `Old`                                   | factor (2 levels)        |
| `Treatment`       | `Control`, `Rapa`, `MetR`, `Rapa+MetR`           | factor (4 levels)        |
| `Metagroup`       | `Y_C`, `O_C`, `O_M`, `O_R`, `O_R_M`, …          | crossed factor (8 levels)|

Natural join key: `Sample`. Note that `PoolMultiSample` could also be used to verify `pool_id` derivation from the path (a useful integrity check).

## Common header patterns to suggest in Pause 2

When the agent enumerates candidate join keys, sort them by likelihood:

1. Exact match to `sample_id` (case-insensitive)
2. `Sample`, `SampleID`, `LibraryID`, `MultiSampleID` (header substring "sample" or "library" + "id")
3. The first column whose unique-value count equals `obs['sample_id'].nunique()`
4. A composite of two columns whose concatenation has the right cardinality

The agent **suggests** a candidate but never auto-selects. The user names the column.

## Validation rules to enforce before joining

After the user picks a column, run before joining:

```python
n_samples = adata.obs["sample_id"].nunique()
metadata = pd.read_csv(metadata_path, sep="\t")     # or "," for CSV
key_values = metadata[chosen_column].astype(str)

# 1. Uniqueness in metadata
assert key_values.is_unique, f"Join key '{chosen_column}' has duplicates"

# 2. Cardinality match
unmatched_in_obs = set(adata.obs["sample_id"].unique()) - set(key_values)
unmatched_in_meta = set(key_values) - set(adata.obs["sample_id"].unique())
if unmatched_in_obs:
    print(f"WARNING: {len(unmatched_in_obs)} samples in AnnData not in metadata: {sorted(unmatched_in_obs)[:5]}...")
if unmatched_in_meta:
    print(f"WARNING: {len(unmatched_in_meta)} metadata rows not in AnnData: {sorted(unmatched_in_meta)[:5]}...")
```

Surface unmatched sets to the user before joining; do not silently left-join past a cardinality mismatch.

## How the join lands in `obs`

After validation:

```python
obs_with_meta = adata.obs.merge(
    metadata, left_on="sample_id", right_on=chosen_column, how="left", validate="many_to_one"
)
adata.obs = obs_with_meta.set_index(adata.obs.index)
```

`validate="many_to_one"` is load-bearing — pandas raises if `sample_id` does not uniquely map to a metadata row. Catch the exception, surface, do not silently degrade to `validate=None`.

## Where the recorded decision lives

`analysis_config.yaml::decisions::cellranger-multi-to-anndata::metadata`:

```yaml
decisions:
  cellranger-multi-to-anndata:
    metadata_source: "00_data/raw/SamplesMetadata.txt"
    metadata_join_key: "Sample"
    species: "mouse"
```

Re-runs read this block and skip the corresponding pauses unless the input data shape (`n_samples`, `n_pools`) changed since the recorded run.
