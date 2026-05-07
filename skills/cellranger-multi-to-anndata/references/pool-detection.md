# Pool detection — glob and regex patterns

`cellranger multi` writes per-sample matrices under `per_sample_outs/<sample>/count/sample_filtered_feature_bc_matrix.h5`. The pool name is the parent directory of `per_sample_outs/`. The glob patterns below cover the layouts observed across reference projects; the skill exposes `cellranger_root` and `pool_pattern` as parameters so a project with a non-standard layout can be handled without editing the skill.

## Canonical glob

```python
import glob
h5s = sorted(glob.glob(
    f"{cellranger_root}/**/per_sample_outs/*/*/sample_filtered_feature_bc_matrix.h5",
    recursive=True,
))
```

Two `*` segments under `per_sample_outs/`:

- First `*` is the **sample name** (CellRanger's per-sample output dir).
- Second `*` is `count` (CellRanger's per-feature-type subdir; `count` for Gene Expression, `vdj` for VDJ, etc.). The fixed `count` value can be substituted for the second `*` if you want to skip non-GEX features at the glob level: `per_sample_outs/*/count/sample_filtered_feature_bc_matrix.h5`.

## Sample-id and pool-id extraction

Anchor on the directory structure, not on the filename. Two reference patterns:

```python
import re
sid  = re.search(r"per_sample_outs/([^/]+)/", fp).group(1)   # always-true on cellranger multi
pool = re.search(r"<root_segment>/([^/]+)/", fp).group(1)    # depends on layout
```

The `<root_segment>` choice depends on the project's directory naming:

| Reference project | Root segment in path | Example pool name |
|-------------------|----------------------|-------------------|
| 13403-YD          | `Counts`             | `13403-YD-P1`     |
| 14616-DM          | `CellRanger`         | `14616-DM-P1`     |

If the pool segment is not deterministic from the path, the skill's `build()` accepts a `pool_pattern` argument (a regex with one capture group). For projects where pools are flat siblings of `per_sample_outs/`, a simpler `os.path.basename(os.path.dirname(...))` walk is also valid; the regex is preferred because it tolerates extra path depth.

## Filtered vs raw barcodes

The default file is `sample_filtered_feature_bc_matrix.h5` — CellRanger's filtered output (cell-called barcodes only). Switch to `raw_feature_bc_matrix.h5` when downstream ambient-RNA correction (SoupX, CellBender) needs the unfiltered droplet distribution. The `build()` script accepts `filtered_or_raw="filtered" | "raw"`; choose at ingestion, do not flip mid-stream.

## Validation after discovery

Print and verify before reading any matrices:

```python
print(f"Found {len(h5s)} matrices")
for fp in h5s[:3] + h5s[-3:]:
    sid  = re.search(r"per_sample_outs/([^/]+)/", fp).group(1)
    pool = re.search(rf"{root_segment}/([^/]+)/", fp).group(1)
    print(f"  {pool} / {sid}")
```

Expected: `len(h5s) == n_pools × samples_per_pool`. If off by one or more, check for:

- A pool that is still mid-run (no `per_sample_outs/`)
- A sample that failed CellRanger filtering (no filtered matrix; the raw matrix may still exist)
- A symlink pointing outside `cellranger_root` (resolved by `recursive=True` only if `**` is used)
