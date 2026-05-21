# Footprint precompute (Phase A.6)

## When to run

The TF Footprinting tab in ShinyMultiome.UiO reads from `pbmc@assays[["peaks"]]@motifs` and `pbmc@assays[["peaks"]]@positionEnrichment`. Both must be populated for the tab to work. Phase A.6 wraps the precompute.

The skill defaults A.6 to **OFF** because:

- It is computationally expensive (1–4 hours depending on motif count and group.by cardinality).
- The motif list and `group.by` choice need wet-lab input.
- Many wet-lab deployments do not need footprints (the Coverage tab is the main attraction).

Run Phase A.6 only when the wet-lab brief explicitly mentions TF footprinting.

## What `AddMotifs` does

```r
obj <- Signac::AddMotifs(obj, genome = bsg, pfm = pfm, assay = "peaks")
```

Steps:

1. For every PFM in `pfm`, scan every peak in `obj[["peaks"]]` for matches (using `motifmatchr::matchMotifs` under the hood).
2. Store a sparse `peaks × motifs` binary matrix in `pbmc@assays[["peaks"]]@motifs@data`.
3. Populate `@motif.names` (named character vector, PFM_id → display name).

Cost: O(n_peaks × n_PFMs). For 200k peaks × 600 PFMs (JASPAR2020 vertebrate core), ~10–30 min on a 4-core machine.

## What `Footprint` does

```r
obj <- Signac::Footprint(obj, motif.name = motifs_of_interest,
                         genome = bsg, group.by = "cell_type", assay = "peaks")
```

Steps:

1. For each motif in `motif.name`, find all peaks containing a hit.
2. Pull fragment-level coverage around each hit, stratified by `group.by`.
3. Compute per-position, per-group enrichment.
4. Store as a `SummarizedExperiment` in `pbmc@assays[["peaks"]]@positionEnrichment[[motif]]`.

Cost: O(n_hits × n_groups × motif_window_size). For 50 motifs × 20 cell types × 50k cells, ~30–90 min on a 4-core machine.

## Picking motifs

Three strategies:

### A. JASPAR2020 vertebrate core (default)

```r
pfm <- TFBSTools::getMatrixSet(
  JASPAR2020::JASPAR2020,
  list(species = 9606, all_versions = FALSE)   # 9606 = human; 10090 = mouse
)
```

~600 PFMs for human, ~200 for mouse. Comprehensive, slow. Use for an unbiased exploration.

### B. Curated marker TFs

```r
marker_tfs <- c("PAX5", "EBF1", "TCF7", "GATA3", "FOXP3", ...)
pfm <- pfm[names(pfm) %in% marker_tfs]
```

Hand-curated list of ~10–50 lineage-defining TFs. Fast, biology-driven. Default for wet-lab deliverables.

### C. Top differential TFs from upstream `chromVAR` or `SCENIC` analysis

If you've already run differential motif accessibility (chromVAR) or regulon analysis (SCENIC+) and have a ranked list of TFs of interest, footprint only those.

The skill does not bundle TF lists — the Phase A.6 wrapper accepts `motif_names` as a vector argument. Decision Pause S6 prompts the user.

## Group.by choice

`Footprint` stratifies by a metadata column. Decision: which column?

- `cell_type` — coarse (8–20 levels typical). Fast. Plots are readable.
- `leiden_<resolution>` — finer (20–80 clusters). Slower. Plots can be busy.
- `sample_id` — between-sample comparison; rarely useful for footprints (which are cell-type biology, not technical).

**Default: `cell_type`** for the curated-marker-TF strategy. Override via `footprint_group_by` argument.

## Storage cost

The `@positionEnrichment` slot adds ~50–200 MB per (motif × group.by) combination. For 50 motifs × 20 cell types: ~5–10 GB on disk. The `.rds` size grows accordingly.

To strip footprints later (e.g., re-deploy without them):

```r
obj@assays[["peaks"]]@positionEnrichment <- list()
saveRDS(obj, ...)
```

## Verification

```r
length(obj@assays[["peaks"]]@motifs@motif.names) > 0   # AddMotifs ran
length(obj@assays[["peaks"]]@positionEnrichment) > 0   # Footprint ran
names(obj@assays[["peaks"]]@positionEnrichment)        # which motifs
```

The validator (`validate_signac_rds.R --require-motifs`) tests both slots when called with that flag. Set `.env::FOOTPRINTS_ENABLED=true` to enable that mode in the standard validation step.

## Iteration pattern

```
1. Run Phase A without A.6 → deploy → wet lab uses Coverage tab for 1 week.
2. Wet lab requests footprints for TFs of interest X.
3. Re-run Phase A.6 with motif_names = X.
4. saveRDS, docker compose down, docker compose up -d.
5. Wet lab refreshes; Footprint tab populated.
```

This 4-step loop is faster than blocking the initial deploy on footprint precompute. Document this in the wet-lab onboarding so they know "footprints come second".
