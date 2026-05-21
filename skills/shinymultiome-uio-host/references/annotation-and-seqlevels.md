# Annotation and seqlevels — UCSC vs Ensembl

## The single most common failure mode

> Coverage tab opens, gene name is accepted, spinner runs ~5s, blank panel renders, no error.

Cause: `Annotation(obj[[atac_assay]])` has Ensembl-style seqlevels (`1, 2, ..., X, Y, MT`), while fragment files (and BSgenome) use UCSC-style (`chr1, chr2, ..., chrX, chrY, chrM`). Signac's `CoveragePlot` looks up the gene's chromosome in `Annotation`, gets `1`, asks the fragment tabix index for region `1:start-end`, gets zero rows, renders an empty panel. No error path.

## The fix

Phase A.2 always converts to UCSC:

```r
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "mm39"   # or hg38, mm10, etc.
Annotation(obj[[atac_assay]]) <- annotation
```

`seqlevelsStyle<-` is from `GenomeInfoDb`; it relies on a curated mapping from Ensembl/NCBI to UCSC for known assemblies. mm39 / GRCm39 is mapped (added to GenomeInfoDb in 1.30.x).

## Why Signac doesn't do this automatically

`Signac::GetGRangesFromEnsDb` returns Ensembl-style seqlevels by default — that's how EnsDb stores them. Signac trusts the user to convert. The default vignette (`pbmc_multiomic.html`) explicitly does:

```r
annotations <- GetGRangesFromEnsDb(EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'   # <-- the load-bearing line
genome(annotations) <- "hg38"
Annotation(pbmc[["peaks"]]) <- annotations
```

Many real-world analysis scripts skip the `seqlevelsStyle<-` line, because it works fine inside an interactive R session where `CoveragePlot` errors are swallowed silently. The bug surfaces in production when a wet-lab user opens the Coverage tab and sees nothing.

## The mismatch can be the other way around

If the analysis was done with `seqlevelsStyle = "Ensembl"` *and* the fragment files were produced with Ensembl-style barcodes (rare; Cell Ranger ARC outputs UCSC-style), then forcing UCSC here will break the *fragments* lookup instead. Phase A.2 has a defensive check: if `Annotation` is already UCSC-styled and the fragments look UCSC-styled too, leave alone. Override the conversion via `seqlevels_style = "Ensembl"` in `prepare_for_shinymultiome()` if needed.

## chrM vs chrMT vs MT

Mouse mitochondrial seqlevel naming: UCSC uses `chrM`, NCBI uses `MT`, Ensembl uses `MT`. After `seqlevelsStyle <- "UCSC"` you get `chrM`. Cell Ranger ARC fragment files use `chrM`. Most mouse multiome analyses do not load mitochondrial peaks anyway — but if the dataset contains `chrM` peaks, the alignment is correct.

## seqinfo coverage

Some EnsDb releases include scaffolds and unplaced contigs (`KI270751.1` etc) as well as the canonical chromosomes. After `seqlevelsStyle<-`, scaffolds may be dropped if they have no UCSC mapping. This is fine — wet-lab tracks are looked up by gene symbol on canonical chromosomes; scaffolds are not addressable from the UI dropdown.

## Idempotency

Phase A.2 is idempotent. If `Annotation(obj)` is already populated and UCSC-styled, the function returns the object unchanged. If the user runs Phase A multiple times (e.g., after re-saving the .rds with a new EnsDb release), only the first run rebuilds; subsequent runs are no-ops.

## Validating after Phase A

```r
ann <- Annotation(obj[["peaks"]])
stopifnot(length(ann) > 0)
stopifnot(startsWith(as.character(seqnames(ann))[1], "chr"))
stopifnot(unique(genome(ann))[1] %in% c("mm39", "hg38", "mm10", "hg19"))
```

The `validate_signac_rds.R` check encodes these assertions.
