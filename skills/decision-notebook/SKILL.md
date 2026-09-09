---
name: decision-notebook
description: "Router for live analysis notebooks. Pick the flavor by decision: an annotation campaign instrument (marimo + jscatter, selection manifest, rounds, cleaning; specified here), gate sign-off (decision-gate-notebook), one live look (interactive-breakpoint-explorer), or freestyle EDA. Use when work needs a notebook and the first question is which kind."
license: MIT
---

# Decision Notebook

## Overview

A notebook earns its place at a decision. The pipeline computes; a person looks, selects,
rules; the ruling lands in a durable record and the pipeline moves. This skill routes to the
notebook shape the decision calls for, and fully specifies the heaviest shape, the
annotation campaign instrument.

## Route by need

```
The analysis needs a live surface. What is the decision?
│
├─ Relabel or drop cells across many sessions, evolving label space
│    → Flavor A: annotation campaign instrument (this skill)
├─ Approve or reject what a stage just wrote (freeze a set, pick a threshold)
│    → skill: decision-gate-notebook (read-only Quarto review; verdict in analysis_config)
├─ Eyeball one inflection point once, then move on
│    → skill: interactive-breakpoint-explorer (live-kernel qmd, jscatter brush)
└─ Look around, nothing downstream waiting
     → Flavor D: freestyle EDA (below)
```

The moment a freestyle look starts producing selections or verdicts, graduate it to A or to
the gate pattern. Whatever the flavor, the call you reach goes to the durable record: the
selection manifest, `decisions.<stage>` in config, or a stage-keyed note under
`docs/_internal/` (skill: reasoning-trace).

## Flavor A: the annotation campaign instrument

One marimo notebook serves an entire campaign: rounds of lasso, rule, and marker selections
over per-lineage re-embeddings, rationale captured at save time, in-kernel evidence on
demand, and a deterministic path from saved selections to a cleaned, relabelled roster.

Distilled from one snMultiome campaign that ran it end to end: a lens compute and export
stage, a selection resolution stage, and a roster cleaning stage, plus the explorer notebook
they feed. That tree is history, not a path to open; everything below is the
dataset-agnostic contract.

### Data plane

- A compute stage owns every artifact the notebook reads. The notebook renders and selects;
  it computes embeddings for nobody.
- Per lens, per round: one explorer parquet, one row per cell. Columns: barcode, lens
  de-novo UMAP x/y, frozen-integration coordinates for cross-reference, cluster ids at
  ~3 leiden resolutions, label_fine + label_coarse, condition/sample keys, QC
  (genes, counts, mito), doublet fields, a small `expr_<gene>` marker panel, and per-cell
  activity scores (signed regulon AUCs, motif deviations) joined by barcode.
- Sidecars for runtime streaming: the full matrix as CSC npz plus an index json
  (`{"barcodes": [...], "genes": [...]}`), one per modality (RNA log1p CP10K; ATAC gene
  activity where paired). Stage-exported, guarded so reruns skip the export.
- Rounds live in config: `rounds."N"` declares a roster (null means full) and the explorer
  parquet per lens. A `ROUND` env var selects. Round refresh recipe: resolution stage,
  cleaning stage, then the lens stage with `FORCE=1`.
- Every path comes from `analysis_config.yaml`. The notebook holds zero literals.

### Panel architecture

- One creator cell builds all scatter panels; selection observers attach to the raw widgets
  and mutate traits directly, so overlay updates leave the grid unrendered and lasso stays
  fluid.
- Grid semantics: the A row shows the active lens frame on its de-novo sub-embedding, the
  B row shows the whole roster on the round's global re-embedding, and one panel shows the
  frozen integration embedding with membership opacity, so a selection reads in all three
  coordinate systems at once.
- Every panel carries an identity header rendered in HTML above the canvas: panel id, frame
  and embedding, nucleus count, active color channel, and the color key (swatches or a
  colorbar). The widget's own legend draws over the canvas; the header replaces it.
- Colors resolve from config `colors.cell_types`: label-keyed hex, one hue per coarse
  family, fine states as shades within the family, colorblind-safe anchors (Okabe-Ito,
  Tol). The same label carries the same color in every lens, full or zoomed. Categoricals
  outside the map (cluster ids, condition) cycle a fallback palette.
- One dot-size rule lives inside the single panel-creator function, scaled by point count
  so sparse zoomed panels carry the same visual weight as the dense full view.

### Selection machinery

- Selection state is a tri-state tuple (barcodes, sequence, origin); the latest selection
  wins. Origin is lasso, rule, or marker; scope is subset or full; both persist as
  `source_brush`.
- Rule selection picks a categorical column and level set, so a whole cluster or a boolean
  flag (doublet true) selects without lassoing.
- Saving writes a per-selection CSV under `03_results/interactive/selections/` and a
  manifest row: selection_id, lens, action (drop | relabel | keep), new_label, priority,
  rationale, n_cells, saved_at_utc, notebook, source_brush, status (active | superseded).
- Saves append by default (union on barcode), Replace is an explicit mode, every touch
  snapshots the prior CSV into `_history/`, and deletion is its own panel. A same-name save
  must leave earlier work recoverable.
- Priority follows specificity: the smaller selection outranks the larger one it nests
  inside; at ties, severity orders drop > relabel > keep. The resolver emits a conflicts
  table and downstream stages hold while it is non-empty.
- Labels key to barcodes alone. Cluster ids renumber on every rebuild (membership shifts,
  HVGs shift, the graph re-forms), so a cluster-keyed rule relabels the wrong cells after
  one rebuild.

### In-kernel evidence

- DE at the moment of doubt: Mann-Whitney AUC with BH correction, both directions, on CSC
  slices; contrasts are selection vs the rest of its lens, group vs group, or an arbitrary
  pair. Seconds, in the kernel, while the lasso is still live.
- Gene streaming: any gene pulls from the sidecar at runtime, joins both frames, and
  becomes a color channel.
- Marker co-detection: two runtime genes on a scatter, each axis RNA or ATAC, a quadrant
  blend (double-positive, two singles, double-negative) painted back onto the embedding,
  and a lasso on the scatter that selects like any other brush.
- Genomic tracks, paired-ATAC datasets: tabix-indexed merged fragments, per-group Tn5
  cut-site pileups over a chosen gene window, a kb-offset axis anchored at the
  strand-correct TSS, region rails (own atlas, cCREs by class, external OCR sets), and
  lane-packed gene/exon models. RNA-only datasets omit this section and keep the rest.

### Downstream resolution

- A resolution stage turns the manifest into per-barcode verdicts by severity then
  priority, and writes the conflicts table.
- A cleaning stage applies the ladder: owner selections first, then a debris rule with a
  rescue clause, then a doublet rule; then owner relabels; then `label_parents` maps every
  fine label to its coarse compartment, aborting on an unmapped label. Outputs: disposition
  table, summary, roster CSV with label_coarse + label_fine.
- The lens compute stage re-embeds each lens on the new roster with a shared recipe (HVGs
  on counts, PCA, neighbors, leiden at ~3 resolutions, fixed seed) and re-exports parquets.

### Serving

- Headless marimo on the project's assigned port, bound to 0.0.0.0; the SSH tunnel targets
  the container IP, not host loopback (skill: container-port-tunnel). Access token from the
  server log.
- Kill the server before editing the notebook .py: marimo autosave rewrites cells under
  you and patch anchors drift. Patch, syntax-check, relaunch detached (setsid), read the
  fresh token.

### Traps

- Cluster ids are nonstationary across rebuilds; barcodes are the only stable key.
- A same-name save that overwrites destroys campaign work; append-by-default plus
  `_history/` snapshots is the floor.
- Binary "activity > 0" ATAC gates pass almost everything; use quantitative fractions of a
  panel.
- snRNA-seq detects long intron-rich genes preferentially; check a novel split's markers
  against gene length before believing it.
- The widget legend draws inside the canvas; identity headers in HTML keep the canvas
  clean.

## Flavor B: gate sign-off

A stage wrote artifacts and a person must approve them before an expensive or irreversible
next stage. Use skill `decision-gate-notebook`: read-only Quarto review, verdict recorded
under `decisions.<stage>` in `analysis_config.yaml`.

## Flavor C: one live look

One inflection point, one session, evidence out. Use skill
`interactive-breakpoint-explorer`: live-kernel Python qmd, linked jscatter panels, persisted
barcodes plus a labelled snapshot.

## Flavor D: freestyle EDA

Looking around with nothing downstream waiting. A marimo notebook or live qmd over a
compact parquet export, paths from config, figures disposable. Skip the persistence
machinery. When a selection wants saving or a verdict forms, move to Flavor A or the gate
pattern and give the call a durable home.
