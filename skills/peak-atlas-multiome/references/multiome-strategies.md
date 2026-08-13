# The four multiome strategies in depth

True 10x Multiome measures RNA and ATAC in the SAME nucleus, so one barcode can be grouped four independent ways. Peaks are discovered separately inside each grouping, and agreement across the four groupings becomes the cross-strategy confidence `n_strategies`. This file details how each grouping is built and why four. Provenance: `02_analysis/S2_6b_atac_peak_recalling_multistrategies.R` (the cluster-creation block + `call_peaks_by_strategy`).

## Why four, not one

Calling peaks on the whole dataset lets abundant populations dominate: a rare type's accessible regions never reach the pileup needed to call a peak. Grouping cells first — by transcription, by chromatin, jointly, and by expert label — means each population is peak-called in its own pseudobulk, so rare-population accessibility is preserved. And because the four groupings draw on different evidence (RNA expression, ATAC accessibility, the joint WNN graph, curated biology), a peak found by several of them is more trustworthy than one found by a single grouping. That redundancy is the point: it is both protection (rare peaks survive) and confidence (cross-modal agreement scores peaks).

## Strategy 1 — RNA clusters (`RNA_clusters`)

Cluster on transcription:

```r
DefaultAssay(obj) <- "RNA"
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30, graph.name = "RNA_snn")
obj <- FindClusters(obj, resolution = 0.5, graph.name = "RNA_snn")
obj$RNA_clusters <- obj$seurat_clusters
```

`pca` dims `1:30` is the standard RNA neighbor space; `resolution = 0.5` gives coarse, robust clusters (you want stable populations to pseudobulk, not over-split fragments). These groups capture cells that are transcriptionally alike even if their chromatin embedding smears them together.

## Strategy 2 — ATAC clusters (`ATAC_clusters`)

Cluster on chromatin:

```r
DefaultAssay(obj) <- "ATAC"
obj <- FindNeighbors(obj, reduction = "lsi", dims = 2:30, graph.name = "ATAC_snn")
obj <- FindClusters(obj, resolution = 0.5, graph.name = "ATAC_snn")
obj$ATAC_clusters <- obj$seurat_clusters
```

LSI dims `2:30` (NOT `1:30`): LSI component 1 in scATAC almost always correlates with sequencing depth rather than biology, so it is dropped from the neighbor graph (the same reason the re-embedding step runs a `DepthCor` check). These groups capture cells with similar accessibility even when RNA is sparse.

## Strategy 3 — WNN clusters (`WNN_clusters`)

The joint grouping — and the one that is only possible with paired data. Weighted-nearest-neighbor (Hao et al. 2021) learns a per-cell weighting of the RNA vs ATAC modality and builds a joint `wsnn` graph:

```r
# If a wsnn_res* column already exists from upstream WNN, reuse it:
wsnn_cols <- grep("^wsnn_res", colnames(obj@meta.data), value = TRUE)
if (length(wsnn_cols) > 0) {
  obj$WNN_clusters <- obj@meta.data[[wsnn_cols[1]]]
} else {
  obj <- FindMultiModalNeighbors(obj, reduction.list = list("pca", "lsi"),
                                 dims.list = list(1:30, 2:30))
  obj <- FindClusters(obj, graph.name = "wsnn", resolution = 0.5)
  obj$WNN_clusters <- obj$seurat_clusters
}
```

WNN groups cells that are consistent across BOTH modalities, sharpening boundaries that either modality alone blurs. It cannot be computed without paired cells — this is the strategy that distinguishes the multiome regime from the unpaired one. Never substitute RNA or ATAC clusters for WNN; doing so collapses two of the four "independent" votes into a correlated pair and inflates apparent support.

## Strategy 4 — CellType (expert annotation)

The curated biological labels (`refined_cell_type` or your equivalent annotation column). These are not derived from a single embedding; they fold in marker knowledge, doublet review, and manual rescue. They anchor the consensus in known biology and let rare, hand-curated populations (e.g. HSC, pDC) contribute their own pseudobulk even when an unsupervised clustering merged them into a neighbor.

```r
# whichever curated column exists, in priority order:
celltype_col <- c("refined_cell_type", "refined_v2d", "rough_annot_v6")
celltype_col <- celltype_col[celltype_col %in% colnames(obj@meta.data)][1]
```

## How they feed the pipeline

Each strategy is run through `call_peaks_by_strategy` (per-group pseudo-replicate calling) and merged independently into one GRanges. The four merged sets then vote: `calculate_strategy_support` counts, for each candidate peak, how many of the four overlap it (`n_strategies`, 1..4), and that count drives the score boost, the metadata pre-filter, and the rescue arm. See `scripts/call_peaks_multistrategy.R` and `peak-atlas-framework/references/support-voting.md`.

## See also

- `peak-atlas-framework/references/support-voting.md` — the compute-before-merge rule and the `adjusted_score` boost.
- `references/primary-rescue-filter.md` — how `n_strategies` gates the rescue arm.
- `scripts/call_peaks_multistrategy.R` — the runnable caller + reconciliation.
