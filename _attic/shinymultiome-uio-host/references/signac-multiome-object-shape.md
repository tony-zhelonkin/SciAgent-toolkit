# Signac multiome object shape — what ShinyMultiome.UiO expects

ShinyMultiome.UiO is a thin Shiny wrapper over Signac plotting calls. Its `server.R` makes the following hard assumptions on the loaded Seurat object (the `pbmc` global):

1. **Variable name = `pbmc`** — the loaded object lives under exactly that name. The skill's overlay `global.R` does `pbmc <- readRDS(...)`, so the .rds can be saved under any variable name; the overlay reassigns.
2. **ATAC assay named `peaks`** — `DefaultAssay(pbmc) <- "peaks"` is hard-coded for the Coverage and Footprint tabs. The Phase A.1 step `normalize_assay_names()` renames whatever the .rds calls its ChromatinAssay to `peaks` (configurable via `.env::ATAC_ASSAY` if you want a different label, but you must change it in BOTH `.env` and upstream — don't, just rename).
3. **RNA assay named `RNA` and (optionally) `SCT`** — Feature-of-Interest tab switches assays by name. The dropdown is populated from `Assays(pbmc)`; missing assays simply omit options.
4. **Default reductions present** — `umap`, `lsi`, `pca`. The Clustering tab's `dimRed` dropdown probes `Reductions(pbmc)`. Missing reductions hide options.
5. **`@meta.data` columns** — the `cellInfo` dropdown on the Clustering tab populates from `colnames(pbmc@meta.data)`. The Coverage tab's `group.by` dropdown comes from the same source. The skill's `.env::GROUP_BY_COLUMNS` is the curated subset; Phase A.4 validates each is present.

## ChromatinAssay invariants

For the Coverage tab to render anything:

```r
assay <- pbmc[[atac_assay]]               # ChromatinAssay (Signac)
length(Annotation(assay)) > 0             # gene track source
length(Fragments(assay)) > 0              # coverage source
all(seqlevels(Annotation(assay)) %in% seqnames(BSgenome))   # mm39 / hg38 alignment
```

The seqlevel alignment between `Annotation` (the gene model) and `BSgenome` (the genome sequence loaded by `global.R`) is the load-bearing invariant. Mismatch produces empty plots, no error. Phase A.2 enforces UCSC seqlevel style on `Annotation`; the BSgenome packages all use UCSC by convention.

## Fragment object internals

```r
frag <- Fragments(assay)[[1]]
frag@path        # absolute path stored in the .rds
frag@hash        # MD5 of the .tsv.gz, computed at object creation
frag@cells       # named character vector mapping fragment-file barcodes -> object cell names
```

The `@hash` is **not** validated at load time by Signac; you can rewrite `@path` to a different file with the same name and Signac will not notice. This is what makes the path-rewrite (Phase A.3) safe — we move the file (or update the path to a bind-mount) without recomputing the hash.

The `@cells` mapping is the join key between the fragment file and the Seurat object. If the fragment file has barcodes like `AAACGCT-1` but the object cells are `sample1_AAACGCT-1` (sample-prefixed), `@cells` is the lookup. Phase A.3 does **not** modify `@cells` — only `@path`.

## Multi-sample fragment handling

A multiome object may have multiple `Fragment` objects in `Fragments(assay)`, one per sample. Upstream `global.R` only handles the single-sample demo case (`frags[[1]] <- UpdatePath(...)`). The skill's overlay iterates all of them:

```r
frags <- Fragments(pbmc)
for (i in seq_along(frags)) {
  base     <- basename(GetFragmentData(frags[[i]], slot = "path"))
  new_path <- file.path(FRAGMENTS_DIR, base)
  frags[[i]] <- UpdatePath(frags[[i]], new.path = new_path)
}
Fragments(pbmc) <- NULL
Fragments(pbmc) <- frags
```

The `Fragments(pbmc) <- NULL` step is required: Signac's setter does not replace; assigning a list to `Fragments<-` appends. Clear-then-set ensures the new paths replace the old.

## Links (Peak2Gene)

Optional. `Links(assay)` returns a `GRanges` of peak–gene associations. CoveragePlot draws these as arcs at the bottom of the panel when `links = TRUE`. Empty `Links` is fine — the arcs just don't render.

Compute via `Signac::LinkPeaks(pbmc, peak.assay = atac_assay, expression.assay = rna_assay)`. Heavy: 10–60 min depending on cell count and gene count. Phase A.5 wraps this with an idempotency check — re-running on an object that already has Links is a no-op.

## Motif slot (footprint tab)

```r
pbmc@assays[["peaks"]]@motifs                # Motif object — names the PFMs
pbmc@assays[["peaks"]]@motifs@motif.names    # named character vector — drives Footprint dropdown
pbmc@assays[["peaks"]]@positionEnrichment    # list of footprint enrichment SummarizedExperiments
```

Both must be populated for the Footprint tab to work. `AddMotifs()` populates the first; `Footprint()` populates the second. Phase A.6 wraps both. Default `motif_names` is the first 20 PFMs from JASPAR2020 vertebrate core — adjust via the `motif_names` argument.

## What the skill does NOT modify

- Counts matrices (RNA, peaks)
- Reductions (umap, lsi, pca)
- Idents
- Embeddings
- Cell barcodes (object's row order is preserved)

Phase A is a *configuration overlay* on the .rds — assay names, Annotation seqlevels, fragment paths, optionally Links/motifs. The biology stays the same.
