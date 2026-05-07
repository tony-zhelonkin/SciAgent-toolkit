# Annotations autosave — `--annotations-dir` semantics

CellxGene's autosave is what lets the wet lab open a CXG instance, label some clusters, click around, close the tab, and find the labels still there tomorrow. It is *also* what fails silently when the volume mount is wrong, the directory is owned by root, or `--annotations-dir` is missing. This document is the operational reference.

## What `--annotations-dir` does

When CXG is launched with `--annotations-dir /annotations`, every label edit in the UI is persisted as a `.csv` under that directory. The naming convention is:

```
<annotations_dir>/<obs_column>-<timestamp>.csv
```

Each row maps `cell_id → label`. The same `obs` column can have multiple snapshots (one per save). On reload, CXG reads the most recent CSV per column.

## The mount contract

`/annotations` inside the container must map to a host-side, persistent, container-writable directory. The reference compose expresses this as:

```yaml
volumes:
  - <<ANNOTATIONS_HOST_PATH>>/<dataset>:/annotations
```

Three things have to be true:

1. **`<<ANNOTATIONS_HOST_PATH>>/<dataset>` exists on the host.** Docker creates parent dirs, but it creates them owned by root if the parent does not exist.
2. **It is owned by `<<UID>>:<<GID>>`.** The compose `user:` field is the container UID; if the host dir is owned by root, writes fail with `EACCES`.
3. **Read-write, not `:ro`.** A common copy-paste mistake from the data mount.

The skill's render step pre-creates the per-dataset annotation dirs on the host (with the right ownership) before `docker compose up`.

## Per-dataset isolation

Each dataset gets its own annotation directory:

```
<<ANNOTATIONS_HOST_PATH>>/
├── full/
├── t_cells/
├── b_cells/
└── myeloid/
```

This is intentional — labels written on the full dataset's CXG instance should not leak into the T-cells subset (where cell_ids are different after the per-celltype re-embed) and vice versa. The mount is per-instance, not shared.

## What goes wrong, and how to detect it

| Symptom | Likely cause | Detection |
|---------|--------------|-----------|
| Labels disappear on container restart | `--annotations-dir` flag missing or pointing inside the container's ephemeral fs | `docker compose exec cellxgene-full ls /annotations` should show CSVs after at least one save |
| Labels save but the CSVs are root-owned | Host dir was created by a previous root-running container | `ls -la <ANNOTATIONS_HOST_PATH>/<dataset>/` from the host |
| "Save" button greyed out in the UI | The container can read but not write to `/annotations` | `docker compose exec cellxgene-full touch /annotations/.test` |
| Some labels are persisted, others vanish | Multiple containers writing to the same host dir | Each `cellxgene-<dataset>` must have a distinct mount target on the host |

## Backing up annotations

Wet-lab annotations are scientifically valuable — they encode the labelling consensus. The autosave dir is a flat collection of CSVs; a daily rsync to durable storage is enough:

```bash
rsync -av --delete <<ANNOTATIONS_HOST_PATH>>/ /backup/cxg-annotations/$(date +%F)/
```

If an annotation conflicts with the AnnData's `obs` after a re-deploy (e.g., the wet lab labelled cells `Treg` but the latest QC re-embed re-numbered cells), the CSV is the source of truth — re-merge on `cell_id`. The `cell_id` index from Phase A is stable across re-deploys *only if the source `.h5ad` was not regenerated with new barcodes*. Document re-deploys carefully.

## Related cellxgene CLI flags

- `--max-category-items 5000` — categories larger than this fall back to free-text labels in the UI; setting it controls the threshold. The shipped value is fine for most projects (typical labelling: ≤ a few hundred categories per column).
- `--disable-diffexp` — *not* set in the reference. The wet lab uses CXG's diff-exp panel; disabling it removes a useful feature for marginal CPU savings.
- `--user-generated-data-dir` — alternative location flag; deprecated. Use `--annotations-dir`.
- `--experimental-annotations` — old beta flag; no longer needed.

## When the wet lab asks "where are my labels?"

The flow:

1. CSVs live at `<<ANNOTATIONS_HOST_PATH>>/<dataset>/`.
2. Each column they edited becomes a CSV; columns are picked from the CXG UI's "User annotations" panel.
3. The CSV format is `cell_id,<column_name>` — easy to merge into the AnnData's `obs` for downstream analysis.
4. If the labels look sparse (only a fraction of cells have values), that is the wet lab's choice — they only labelled what they were sure of.

To bring annotations into Python:

```python
import pandas as pd, scanpy as sc, glob, os
adata = sc.read_h5ad("03_results/objects/08_explore_full.h5ad")
for csv in sorted(glob.glob("03_results/annotation/full/*.csv")):
    col_name = os.path.basename(csv).split("-")[0]
    df = pd.read_csv(csv).set_index("cell_id")
    adata.obs[f"wetlab_{col_name}"] = df.iloc[:, 0]
adata.write_h5ad("03_results/checkpoints/10_with_wetlab_labels.h5ad", compression="gzip")
```
