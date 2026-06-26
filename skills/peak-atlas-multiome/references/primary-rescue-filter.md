# Primary + Rescue: the adaptive quantify-then-prune filter

The consensus caller produces a large candidate atlas (~2.3M peaks in the source bone-marrow multiome). Most of those are real-but-rare or weak; a few are noise that slipped past the pseudo-rep filter. The FILTER stage prunes against ACTUAL cell-by-peak counts, using an adaptive per-cell-type threshold so rare populations are protected by design. Provenance: `apply_primary_rescue_filter` and `calculate_adaptive_thresholds` in `01_scripts/R/peak_utils.R`; wiring in `02_analysis/S2_6b_filter_peak_atlas.R`; tunables in `02_analysis/config/peak_filtration_config.yaml`.

This file covers (1) the quantify-then-prune phase order, and (2) the Primary+Rescue math. The first-principles framing (chase FRiP and marker coverage, not a target count) lives in `peak-atlas-framework/references/rare-celltype-protection.md`; this is the multiome-specific instantiation.

## Quantify-then-prune — the six phases

Quantifying every candidate peak across every cell is expensive, so the prune is staged: cheap metadata first, then a subsampled quantification to decide what survives, then one full quantification of only the survivors.

| Phase | Action | Why |
|---|---|---|
| 1 | **Metadata pre-filter** — keep if `n_strategies >= 2` OR (`n_strategies == 1` AND `score > quantile(score, 0.85)`) | Drops obvious singletons for free; consensus always survives, lone peaks must be top-15% |
| 2 | **Stratified cell sampling** — `n_sample = min(n_type, max(min_per_type, ceiling(0.10 * n_type)))` | A representative subsample where rare types are NOT sampled away |
| 3 | **Subsampled FeatureMatrix** — `Signac::FeatureMatrix(fragments, pre-filtered peaks, sampled cells)` | The count evidence the prune decides on, computed cheaply |
| 4 | **Primary + Rescue prune** (below) | The adaptive decision on the subsampled counts |
| 5 | **Full FeatureMatrix** over ALL cells x surviving peaks, then `CreateChromatinAssay` | The production count matrix and new assay |
| 6 | **Validation** — FRiP, marker promoter retention, TSS floor, embedding, chromVAR | Confirm biology survived (`peak-atlas-framework/references/validation-battery.md`) |

`min_per_type` is `100` in the framework template (the source script used `50`); both protect rare types — pick one and document it.

## PRIMARY — adaptive per-cluster threshold

For each cluster (cell type) of size `n`, the minimum number of accessible cells a peak must reach is:

```
threshold(cluster) = clamp(rate * n, floor, cap)
                   = max(15, min(0.02 * n, 200))
```

with `floor = 15`, `rate = 0.02` (2%), `cap = 200`. A peak is **PRIMARY-kept** if it clears its cluster's threshold in **any one** cluster:

```
cluster_count[peak, cluster] = number of that cluster's cells with count > 0
primary_keep[peak] = any_over_clusters( cluster_count[peak, c] >= threshold(c) )
```

| Cluster | n cells | threshold | as % of cluster |
|---|---|---|---|
| HSC | 147 | 15 | 10.2% |
| pDC | 180 | 15 | 8.3% |
| Neutrophil | 4500 | 90 | 2.0% |

Why the clamp:
- **Floor (15)** gives a rare type a low absolute bar it can actually reach (standard error of a proportion at n=15 is ~0.13, adequate for inference) instead of a percentage of a tiny n that demands almost-universal accessibility.
- **Cap (200)** stops an abundant cluster from imposing an absurd absolute bar (2% of 50k cells would be 1000) that would over-filter common peaks.
- **ANY-one-cluster** is the protection: a peak that is accessible in only HSCs but clears the HSC bar survives, even though it is invisible globally.

## RESCUE — replicated rare biology near the boundary

A peak that FAILS the PRIMARY filter is still rescued if it is high-confidence replicated rarity:

```
peak_presence_global = number of cells (anywhere) with count > 0
rescue_keep[peak] = (n_strategies == 4) & (peak_presence_global >= 10)
```

Found by ALL FOUR independent groupings (RNA, ATAC, WNN, CellType) AND present in at least 10 cells globally. The 4-strategy requirement is strict because it means every modality corroborated the peak; the 10-cell floor (SE ~0.16, minimally acceptable) keeps it from rescuing pure noise. This catches peaks that sit just under a per-cluster bar but are clearly real.

## Final decision and the numbers

```
final_keep = primary_keep | rescue_keep
```

On the source dataset this Primary+Rescue path produced **~544k peaks**. The simpler global-sparsity alternative (`peak_presence > 3` in the subsample, the v1 `sparsity_filter`) produced **~1.15M peaks** on the same candidate atlas. Primary+Rescue is meaningfully tighter — it strips weak globally-thin peaks the global filter keeps — yet protects rare biology the global filter would lose, because the adaptive bar reads each cluster on its own scale and the rescue arm re-adds fully-replicated rare peaks. Tighter AND safer is the whole point.

## Parameterization

All five knobs are arguments to `apply_primary_rescue_filter` (`floor`, `rate`, `cap`, `rescue_global_min`) and `calculate_adaptive_thresholds`, mirroring `peak_filtration_config.yaml::sparsity_filter_v2`. To relax for a shallow library, lower `floor` (e.g. 10) and `rescue_global_min`; to tighten, raise `rate` or lower `cap`. Re-run Phase 6 validation after any change — FRiP retention `>= 0.90` and marker promoter retention `>= 0.95` are the gates that say you did not over-prune.

## See also

- `scripts/apply_primary_rescue_filter.R` — the runnable filter (faithful to the source).
- `checks/check_marker_retention.R` — the rare-marker promoter retention gate.
- `peak-atlas-framework/references/rare-celltype-protection.md` — the shared first principle and config.
- `peak-atlas-framework/references/validation-battery.md` — the full Phase-6 battery.
