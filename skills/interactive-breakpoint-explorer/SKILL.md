---
name: interactive-breakpoint-explorer
description: "The live-kernel house pattern for looking at cells at a pipeline inflection point: a Python Quarto notebook drives linked jscatter panels to brush an embedding, pulls the selection into pandas, and persists barcodes plus a labelled snapshot. Use when a human must eyeball live structure. For the static surface use decision-gate-notebook."
license: MIT
---

# Interactive Breakpoint Explorer

## Overview

At a pipeline inflection point there are calls you cannot make from a summary table,
e.g. *what is that high-mito FOXP3 corner I keep noticing, and is it one donor?* 
You have to **look**, and looking means brushing a live embedding, recoloring it, and pulling the cells you lassoed into
pandas to characterize them. 

This skill is the house pattern for that moment: 
a **live-kernel Python Quarto notebook** that drives [jscatter](https://jupyter-scatter.dev/) linked panels,
extracts a brushed selection, and persists its barcodes plus a labelled snapshot.

It is the **live-kernel sibling** of `decision-gate-notebook`, the static R review whose
`analysis_config.yaml` key a downstream stage guards on.

**Where the decision goes.** This explorer produces evidence: barcodes, a snapshot, a
characterization table. The call you make from that evidence — what you concluded, on what
numbers, and which reading you rejected — belongs in the durable record, as a topic note under
`docs/_internal/<stage-stem>/` (skill: `reasoning-trace`). A notebook cell is a working
surface, and it is overwritten on the next run.

**When to use this skill:**
- A stage just finished and a human needs to *eyeball live structure* — brush a suspicious pocket, compare a few color-bys, lasso cells and see what they are — before an expensive/irreversible next stage.
- You want a durable record of *which cells* were flagged (barcodes) plus a labelled figure, not just a chat description.
- The judgment is visual and interactive; a static PNG grid cannot answer it.

**When NOT to use this skill:**
- A downstream stage must refuse to run until sign-off → use `decision-gate-notebook`, which owns that latch.
- A static, self-contained HTML explorer with no live kernel is enough → use `bulk-rnaseq-pathway-explorer`.
- A plain results/methods write-up with nothing to brush → an ordinary `.qmd` report.

---

## Decision Tree

```
Pipeline stage just finished — a human needs to look before the next stage runs.
│
├─ The look is VISUAL + INTERACTIVE (brush a pocket, recolor, lasso cells, see what they are)
│        → INTERACTIVE BREAKPOINT EXPLORER  (this skill — live jscatter kernel)
│          → persists barcodes + a labelled snapshot
│
├─ A downstream stage must refuse to run until a human signs off
│        → decision-gate-notebook  (static R review; owns decisions.[stage].status)
│
├─ You want a shareable static HTML explorer with no live kernel (plotly, brushable UMAP)
│        → bulk-rnaseq-pathway-explorer
│
└─ You just want to write up results (no decision, no brushing)  → ordinary .qmd report
```

---

## The minimal convention — one selector, one grid, one extraction

Keep an explorer **small**. It is not a 3-notebook pipeline; it is one notebook per inflection
point with exactly three moving parts:

1. **One color-by / annotation selector** — a `CHANNELS = [...]` list you edit and re-run to
   recolor the grid (categorical obs or continuous gene/score columns).
2. **A 2-4 panel linked grid** — `PANELS = grid(df, CHANNELS)`. Brush or lasso any panel; the
   rest highlight the same cells (view/selection/hover synced).
3. **One pandas selection-extraction** — `sel = first_selection()` → `df.iloc[sel]` → a
   selected-vs-all mean table, then `save_selection` + `snapshot`.

That is the whole loop. Resist adding a second and third notebook, per-panel dashboards, or a
widget for every knob — the value is in *looking and lassoing*, not in the UI.

---

## The OOM lifecycle contract (read this first)

**This is the #1 reliability requirement.** jscatter panels are **ipywidgets**: every panel a
cell builds stays alive in the kernel until you call `.widget.close()` on it. A live explorer
session re-runs the grid cell many times while brushing — and each re-run *adds* four more live
panels on top of the old ones. Within a dozen re-runs the accumulated dead widgets **OOM the
kernel**.

The contract that prevents this: **always drive panels through `grid()` / `close_panels()` /
`live_panels()`; never a raw `_mk` loop in a cell.**

```python
PANELS = grid(df, CHANNELS)     # closes ALL prior panels first, then builds + registers the new ones
sel = first_selection()          # reads the live panels grid() registered — no variable-name traps
close_panels()                   # manual escape hatch: free every live panel (a standalone cell)
```

- `grid()` calls `close_panels()` **before** building, so re-running the grid cell never
  accumulates widgets — the count of live panels stays bounded by one grid.
- `first_selection()` with no argument reads `live_panels()`, so `df.iloc[first_selection()]`
  works no matter which panel you lassoed (no reliance on a `PANELS` variable a later cell
  shadowed).
- Keep a **standalone `close_panels()` cell** in the notebook as an escape hatch: if a grid ever
  lags, run it to free every live widget before rebuilding.

**Anti-pattern (do NOT do this):**

```python
# WRONG — leaks a new panel per channel per re-run → kernel OOM
PANELS = [_mk(c) for c in CHANNELS]
jscatter.compose([(s, c) for s, c in zip(PANELS, CHANNELS)], ...)
```

---

## Export-first discipline

The live kernel must stay light. **Never open the multi-GB `.h5ad` checkpoint in the explorer
kernel.** Instead, a project-owned `02_analysis/stages/export_explorers.py` (copy
`assets/export_explorers_skeleton.py` and fill the `# TODO(project):` stubs) materializes a
**compact parquet per inflection point** under `03_results/interactive/`:

```
03_results/interactive/
├── 01_qc_explore.parquet          # x,y + a few QC obs + a handful of markers/scores
├── 02_annotation_explore.parquet  # x,y + frozen labels + markers
└── ...
```

Each table is just 2D coords + a few obs columns + a handful of gene/score columns — indexed by
cell barcode. `load_explorer("01_qc", config=config)` reads it instantly. The exporter is a
**read-only projection** of already-computed checkpoints; it recomputes no biology. Re-run it
whenever the upstream checkpoint changes.

---

## Live-kernel caveat

This is a **live tool** — it needs a running Python kernel (VS Code + Jupyter, or JupyterLab),
not `quarto render`. Brushing/lassoing has no meaning in a batch render; the notebook is a
steering surface you drive by hand. (Its static sibling, `decision-gate-notebook`, *is*
rendered to committed GFM/HTML — that is the artifact you commit for the sign-off.)

---

## The `interactive:` config block

All project-specific lists are config-driven — never hardcoded in the helper lib. Add the
`interactive:` block to `02_analysis/config/analysis_config.yaml` (copy
`assets/interactive-config-snippet.yaml`):

| Key | Default | Governs |
|---|---|---|
| `export_stages` | `[]` | which stage checkpoints the exporter materializes |
| `load_coords_obsm` | `""` | which `obsm` embedding the exporter reads as x,y |
| `save_selection_cols` | `[x, y]` | columns persisted per selection (extend, do not shrink) |
| `categorical_obs` | `[]` | obs columns `grid()` always treats categorical |
| `marker_genes` | `[]` | marker genes the exporter adds as columns |
| `signatures` | `{}` | named gene sets the exporter scores per cell |
| `summary_stats` | `[]` | composition stats in the selection summary (`{name, column, value?}`) |

`grid_height`, `grid_rows`, and `cmap` are also read from this block with sensible floors.

---

## Prerequisites

Run from an analysis project where Scio has materialized
`02_analysis/helpers/interactive_style.py`, mounted the `interactive-style` contract library, and
provided `02_analysis/config/analysis_config.yaml`. Exporters write compact inputs under
`03_results/interactive/` from objects prepared by a committed stage.

---

## Quick Start

```bash
# 1. one explorer per inflection point, in its own folder
mkdir -p 02_analysis/notebooks/01_qc_explore
cp $SKILL/assets/explorer.qmd 02_analysis/notebooks/01_qc_explore/01_qc_explore.qmd

# 2. copy the exporter skeleton into the project and fill the # TODO(project): stubs
cp $SKILL/assets/export_explorers_skeleton.py 02_analysis/stages/export_explorers.py

# 3. add the interactive: block to config
cat $SKILL/assets/interactive-config-snippet.yaml >> 02_analysis/config/analysis_config.yaml

# 4. build the compact tables (once per checkpoint change), then open the .qmd on a LIVE kernel
python 02_analysis/stages/export_explorers.py
```

**Verify it worked:**

```bash
ls 03_results/interactive/*_explore.parquet          # compact tables exist
# in the live kernel: after brushing + save, the barcodes + snapshot landed under the notebook
ls 02_analysis/notebooks/01_qc_explore/eda/selection_*.csv
ls 02_analysis/notebooks/01_qc_explore/eda/*.png
```

---

## Progressive Depth

### Basic usage — the loop

Load a compact table, drive the grid, extract the selection:

```python
from helpers.interactive_style import (
    load_interactive_config, load_explorer, grid, first_selection, save_selection, snapshot)

config = load_interactive_config()
df = load_explorer("01_qc", config=config)

CHANNELS = ["FOXP3", "pct_counts_mt", "leiden", "tissue"]
PANELS = grid(df, CHANNELS, cat=["leiden", "tissue"])   # closes prior panels, builds + registers

# ...brush/lasso a pocket in any panel...
sel = first_selection()                                 # brushed indices from the live grid
sub = df.iloc[sel]
save_selection(df, sel, label="mt_hi_pocket", subdir="01_qc_explore", config=config)
snapshot(df, "FOXP3", name="mt_hi_pocket", subdir="01_qc_explore", indices=sel, config=config)
```

### Intermediate usage — characterize before you save

Between the grid and the save, compare the selection against everything else so the note you
write is grounded:

```python
cols = [c for c in ["FOXP3", "CTLA4", "pct_counts_mt", "n_genes_by_counts"] if c in df.columns]
import pandas as pd
pd.DataFrame({"selected": df.iloc[sel][cols].mean(), "all": df[cols].mean()}).round(3)
```

Name the pocket, save its barcodes (re-loadable later), and snapshot a labelled figure. Then
write what you concluded to `docs/_internal/<stage-stem>/<topic>.md` — the notebook holds the
evidence, the note holds the call.

### Advanced usage — reload a saved selection

Barcodes are the durable currency; re-derive indices in a later session from the saved CSV:

```python
saved = pd.read_csv("02_analysis/notebooks/01_qc_explore/eda/selection_mt_hi_pocket.csv", index_col=0)
sel_reloaded = df.index.get_indexer(saved.index)
sel_reloaded = sel_reloaded[sel_reloaded >= 0]     # keep barcodes still present in df
```

---

## Verification Checklist

After authoring an explorer, confirm:

- [ ] **Export-first:** the notebook loads a `03_results/interactive/*_explore.parquet`, never an `.h5ad` — grep the `.qmd` for `read_h5ad` / `sc.read` (should be absent).
- [ ] **Lifecycle-safe:** panels are built by `grid(...)`, not a raw `_mk` loop; a standalone `close_panels()` cell exists.
- [ ] **Selection persisted:** after a brush + save, `eda/selection_<label>.csv` (barcode-indexed) and `eda/selections_index.csv` exist.
- [ ] **Labelled snapshot:** `eda/<label>.png` + `.pdf` exist with a legend/axes/title (jscatter's own export cannot capture these).
- [ ] **The call is durable:** what you concluded from the brush is a topic note under `docs/_internal/<stage-stem>/`, surviving the next run of the notebook.
- [ ] **Config-driven:** column lists come from the `interactive:` block, not hardcoded in a cell or the helper lib.

---

## Common Pitfalls

### Pitfall: the kernel OOMs after a few grid re-runs

- **Symptom:** the kernel slows, then dies, after you re-run the grid cell a dozen times while brushing.
- **Cause:** each re-run built four new ipywidget panels *on top of* the old ones; jscatter panels stay alive until `.widget.close()`.
- **Fix:** always build via `grid()` (it calls `close_panels()` first) and keep a standalone `close_panels()` escape-hatch cell. Never hand-roll a `_mk` loop.

### Pitfall: `first_selection()` returns empty even though you lassoed

- **Symptom:** `df.iloc[first_selection()]` is empty right after brushing.
- **Cause:** you passed a stale `PANELS` variable a later cell reassigned, or you brushed a panel not in that list.
- **Fix:** call `first_selection()` with **no argument** — it reads the live panels `grid()` registered, independent of any cell's variable name.

### Pitfall: the saved figure has no legend or axes

- **Symptom:** you exported from jscatter's own viewport button and the PNG is a bare point cloud.
- **Cause:** jscatter's viewport export captures only the points — no legend, axes, title.
- **Fix:** use `snapshot(...)`, which renders a labelled matplotlib figure (and styles it through `figure-style` when available).

### Pitfall: the kernel chokes loading the data

- **Symptom:** the setup cell hangs or OOMs on load.
- **Cause:** the notebook opened the multi-GB `.h5ad` checkpoint directly.
- **Fix:** run `export_explorers.py` first and `load_explorer(...)` the compact parquet. Never open the checkpoint in the live kernel.

---

## Resources

- **jupyter-scatter (jscatter):** https://jupyter-scatter.dev/
- **Quarto with Jupyter (Python):** https://quarto.org/docs/computations/python.html
- **Bundled templates:** `assets/explorer.qmd`, `assets/export_explorers_skeleton.py`, `assets/interactive-config-snippet.yaml`

---

## When not to use

- Do not treat a notebook cell as the record. The decision you reach from brushing goes to a topic note under `docs/_internal/<stage-stem>/`; the notebook is overwritten on the next run.
- Do not open the multi-GB .h5ad in the live kernel. Always export a compact parquet under 03_results/interactive/ first (export-first discipline); the kernel loads that, not the checkpoint.
- Do not build panels with a raw _mk loop in a cell. jscatter panels are ipywidgets that leak until closed; always drive them through grid()/close_panels()/live_panels() or you OOM the kernel.
- Do not `quarto render` this as a headless artifact. It needs a live Python kernel (VS Code + Jupyter, or JupyterLab) — brushing has no meaning in a batch render.

---

## See also

- `decision-gate-notebook` — Static sibling; owns the `analysis_config.yaml` latch a downstream stage guards on
- `figure-style` — Prerequisite; the styling/saving contract the labelled matplotlib snapshot goes through before it's written to `03_results/`
- `analysis-code-conventions` — Sibling convention; narrative stages and authoritative data flow
- `reasoning-trace` — Where the decision itself is written, as a durable topic note
- `bulk-rnaseq-pathway-explorer` — Static alternative; a static, self-contained HTML explorer (plotly, no live kernel)
