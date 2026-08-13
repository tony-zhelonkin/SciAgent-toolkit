# Validation battery — does the atlas preserve biology?

A peak set is only as good as the signal it retains. After CALL→MERGE→FILTER→TIER, validate the candidate atlas on four axes before trusting it downstream. Pass/fail thresholds (provenance: `02_analysis/config/peak_filtration_config.yaml`, `validation:`) are below. Runnable gates: `checks/check_peak_atlas.R` (structural) and `checks/check_frip_retention.R` (FRiP).

## 1. FRiP retention (headline metric)

Fraction of reads in peaks, per cell (`scripts/frip.R::calculate_frip`). Compare the filtered atlas to the original.

- **Gate:** filtered median FRiP `>= 0.90` of the original median (`frip_retention_min: 0.90`).
- **Warning band:** `< 0.95` (`frip_retention_warning: 0.95`).
- Losing more than 10% of FRiP means the prune dropped real accessible signal — relax the filter.
- For speed, sample ~1000 cells (`frip_sample_size: 1000`).

## 2. TSS enrichment — a FLOOR, not a discriminator

TSS enrichment measures *library* quality (fragment density at TSSs vs flanks) and is largely **invariant across peak sets** — it does not distinguish a good atlas from a bad one. Use it only as a sanity floor.

- **Floor:** `tss_enrichment_min: 1.0`; warn below `1.5`.
- **Context:** 1.0-1.5 is normal for complex primary tissue; `> 2` is typical for pure cell lines. In the source bone-marrow multiome, TSS held at ~1.37 across unfiltered / per-sample / filtered atlases — confirming it is invariant. Do not "improve" TSS by filtering; you cannot, and trying will cost you FRiP and markers.

## 3. Peaks-per-cell and marker / biology recovery

- **Peaks per cell:** median `>= 500` (`min_peaks_per_cell: 500`). Too few means the atlas is too sparse for stable per-cell embeddings.
- **Rare-marker promoter-peak retention:** `>= 0.95` (`marker_retention_min: 0.95`). For curated rare-type markers (e.g. pDC: Tcf4/Irf8/Flt3/Siglech/Bst2; HSC: Procr/Hlf/Slamf1/Gata2/Meis1; CyclingProg: Mki67/Top2a/Pcna/Ccnb1) check that promoter peaks survive filtering, using a promoter window of `-2000 / +500` (`recovery_options.promoter_window`). Losing rare-marker promoters is the failure the rare-cell-type protection exists to prevent.

## 4. Embedding quality and external correlation

- **Embedding-structure comparison:** recompute LSI / UMAP / WNN per candidate peak set and compare the resulting structure (cluster separation, marker localization). The atlas that best resolves known populations wins — a smaller peak count is fine if structure is preserved.
- **chromVAR TF-motif activity:** TF deviations should track expected lineage TFs (a biological positive control). Route to `chromvar-motif-accessibility` for the computation.
- **Per-lineage promoter accessibility:** lineage marker promoters should be accessible in the matching lineage.

## Structural integrity (always)

`checks/check_peak_atlas.R` (reproduces `validate_peak_atlas`) asserts, exiting non-zero on any failure:
- no overlapping peaks,
- every peak exactly the fixed width (501bp),
- zero blacklist overlaps (when a blacklist is supplied),
- valid coordinates.

Run it on every candidate atlas before the biological checks — structural failures invalidate the rest.

## Suggested gate order

1. `check_peak_atlas.R` — structure (cheap, blocking).
2. `check_frip_retention.R` — FRiP `>= 0.90` (blocking).
3. Marker promoter retention `>= 0.95` and peaks-per-cell `>= 500`.
4. Embedding comparison + chromVAR (interpretive, for choosing among candidates).

## See also

- `scripts/frip.R`, `checks/check_peak_atlas.R`, `checks/check_frip_retention.R`.
- `assets/peak_filtration_config.yaml` — all thresholds, tunable.
- `chromvar-motif-accessibility` — TF activity validation.
