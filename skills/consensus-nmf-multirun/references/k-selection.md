# K-selection — how to read `<variant>.k_selection.png` and pick K

K selection is manual on purpose. Automated elbow-detection on the stability curve picks the right K maybe 60–70% of the time on real-world cell-state landscapes; the 30% miss case produces silently wrong programs that downstream stages amplify. Domain knowledge — "T cells in this experiment should produce ~8 transcriptional states because there are 8 known polarisations" — beats any heuristic.

## What the plot shows

cNMF's `k_selection_plot` produces two stacked panels for K = K_min..K_max:

```
   reproducibility
      │  *
      │   *
      │    * ← max here suggests stable program structure
      │     *
      │      *
      │       * *
      │           *
      └───────────────────  K
            5  6  7  8  9  10 11 12

   stability
      │   *
      │    * *
      │       *
      │        *
      │         * *
      │             *
      │              *
      └───────────────────  K
            5  6  7  8  9  10 11 12
```

- **Reproducibility** (top): how reproducible the program-gene-loadings are across cNMF iterations at this K. Higher = more reproducible.
- **Stability** (bottom): related metric quantifying how stable the program assignments are across iterations. Higher = more stable.

Both metrics fall off at high K (programs become noisy) and rise at low K (programs are too few to capture variation; "stable" here is trivial — there are only 5 of them).

## The inflection-point heuristic

Look for the **largest K** that is still on the high plateau, before the curves start to fall. Call it K*.

In the sketch above, K* ≈ 8 (last K with high-stability and high-reproducibility before the descent at 9–10).

Three failure modes of this heuristic:

- **No clear plateau** — both curves are monotonically decreasing from K_min. Likely K_min was too high; widen the range.
- **Plateau extends past K_max** — both curves are still high at K_max. Widen the range; the data supports more programs.
- **Reproducibility-stability disagreement** — reproducibility favours K=10 but stability favours K=8. Trust stability for downstream merge robustness; the higher reproducibility at K=10 typically reflects redundant programs that will collapse in Stage 5 anyway.

## Predict-before-view

Per the `mentor-mode` skill, before the user opens the plot they state their expected K. Reasoning examples:

- "T cells should have ~5 polarisations (Th1, Th2, Th17, Treg, Tfh) plus a cycling program plus a stress program → I expect K ≈ 7 ± 1"
- "Spleen has lots of cell types but the cell-type signal is in `obs['celltype']`; here I'm looking for *within-celltype* states, so K should be modest, maybe K ≈ 5 ± 1 per subset"
- "Whole-dataset cNMF on a multicellular tissue: K should reflect celltypes + states, K ≈ 10–15"

Then open the plot. If the predicted K matches the inflection, proceed. If not, the plot has new information — but inspect *why* the disagreement exists before changing K. A common cause: a stress-response program that the user did not expect adds K=1; the gene list of program 11 (or wherever it lands) makes that obvious.

## When the plot is unreadable

If both curves are flat across K, cNMF could not stably decompose. Likely causes:

- **Too few HVGs** — `num_highvar_genes=2000` is the default; if the dataset has fewer than ~5000 expressed genes (rare; usually only sparse subsets), HVGs are noisy. Lower to 1000.
- **Too few iterations** — `n_iter=100` is the default; if reproducibility is flat, try `n_iter=200`.
- **Counts not raw** — cNMF needs raw integer counts. Log-normalised input gives flat reproducibility.

The shipped script's pre-flight check verifies `adata.layers["counts"]` exists and has integer-like values; surface the warning early.

## Recording the decision

`analysis_config.yaml::decisions::consensus-nmf-multirun::k_per_variant`:

```yaml
decisions:
  consensus-nmf-multirun:
    k_per_variant:
      full:         10
      fullQC:       10
      T_cells:       8
      T_cellsQC:     8
      B_cells:       7
      B_cellsQC:     7
      Myeloid:      11
      MyeloidQC:    11
```

Re-runs read this and call `cnmf_obj.consensus(k=K, ...)` directly without re-prompting unless the input data has changed.
