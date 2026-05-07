# Subset strategies — Decision Pause 1 worked examples

The five options in Decision Pause 1, with worked examples for projects that do *not* have hardcoded biology like the 13403-YD Th1/Th17 split.

## Option A — Whole dataset only (default)

```python
subsets = {"full": None}     # None means no filter
qc_variants = ["raw"]        # or ["raw", "qc"] from Pause 2
```

**When to choose:** first pass on a new dataset; no prior reason to believe a subset axis matters; compute is tight. Almost always the right starting place; expand only if Stage 5 shows the merged programs are dominated by a single celltype.

**Cost:** 1–2 cNMF runs total.

## Option B — Per major celltype

The skill enumerates `obs[celltype_column].value_counts()` and asks the user to confirm the list. Cells per subset should be ≥ 1000 for meaningful K=6..11; ≤ 1000 cells degrades cNMF's stability.

```python
# Skill output:
# celltype:
#   T_cells       45120
#   B_cells       28310
#   Myeloid       18470
#   Other          6210
#   Unannotated     580          # too small; skill flags as candidate-skip
#
# Subset list confirmed by user: ["T_cells", "B_cells", "Myeloid"]
subsets = {
    "T_cells":  "celltype == 'T_cells'",
    "B_cells":  "celltype == 'B_cells'",
    "Myeloid":  "celltype == 'Myeloid'",
}
```

**When to choose:** clear celltype annotation, multicellular dataset, the hypothesis is that programs differ within celltype.

**Cost:** N + 1 cNMF runs (N subsets + the full); 2(N+1) if Pause 2 picks "both QC".

## Option C — Per condition

```python
# Skill output:
# Group: 2 levels (Young: 50k, Old: 48k)
# Treatment: 4 levels (Control: 25k, Rapa: 25k, MetR: 24k, Rapa+MetR: 24k)
# Metagroup: 8 levels (Y_C, O_C, O_M, O_R, ...)  -- crossed
#
# User picks: Group (2 levels)
subsets = {
    "Young": "Group == 'Young'",
    "Old":   "Group == 'Old'",
}
```

**When to choose:** condition has a known effect on cell-state composition, want condition-specific programs (e.g., "what is the gene program of activated B-cells *in Young*" vs "in Old").

**Cost:** L + 1 runs where L is level count.

**Caveat.** If the factor is balanced (≈ equal cells per level), per-condition cNMF gives interpretable per-condition programs. If unbalanced (one level has 90% of cells), the small-level subset will be cell-count-starved and cNMF stability will degrade — fall back to whole-dataset and use Stage 7's ANOVA to find condition-responsive programs.

## Option D — Per celltype × condition

Cross-product. Most expensive.

```python
# Skill output: 3 celltypes × 2 conditions = 6 subsets
# Cell counts per cell:
#   T_cells × Young: 22500
#   T_cells × Old:   22620
#   B_cells × Young: 14200
#   B_cells × Old:   14110
#   Myeloid × Young:  9180
#   Myeloid × Old:    9290
subsets = {
    "T_cells_Young": "celltype == 'T_cells' & Group == 'Young'",
    "T_cells_Old":   "celltype == 'T_cells' & Group == 'Old'",
    # ... etc.
}
```

**When to choose:** both axes matter, compute budget allows, you have at least ≥ 2000 cells per cell of the cross-product (else fall back to Option B or C).

**Cost:** L × N + 1 runs.

## Option E — User-defined `obs` filter expressions

Free-form. The user names the subset name → expression:

```python
subsets = {
    "Mem_T":         "celltype.str.startswith('T_') & state == 'Memory'",
    "Naive_T":       "celltype.str.startswith('T_') & state == 'Naive'",
    "Activated_M":   "celltype == 'Macrophage' & ifn_score > 0.5",
}
```

**When to choose:** subsets that don't follow a single column. Compound logic, derived columns, score thresholds.

**Validate before running:** for each subset, evaluate the expression and report cell count. The skill's run script does this:

```python
mask = adata.obs.eval(expr)
print(f"  {subset_name}: {mask.sum()} cells")
if mask.sum() < 1000:
    warnings.warn(f"{subset_name} has only {mask.sum()} cells; cNMF stability will be degraded")
```

**Cost:** depends on subset count.

## How the recorded decision lands in config

`analysis_config.yaml::decisions::consensus-nmf-multirun::subset_strategy`:

```yaml
decisions:
  consensus-nmf-multirun:
    subset_strategy:
      mode: per_celltype            # one of: whole, per_celltype, per_condition, cross, custom
      celltype_column: celltype     # populated when mode != whole / custom
      subsets:
        - name: T_cells
          filter: "celltype == 'T_cells'"
        - name: B_cells
          filter: "celltype == 'B_cells'"
        - name: Myeloid
          filter: "celltype == 'Myeloid'"
```

Re-runs read this and skip Pause 1 unless the input AnnData's celltype distribution has changed materially since the recorded run.
