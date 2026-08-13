#!/usr/bin/env Rscript
# check_marker_retention.R - Assert rare-marker promoter-peak retention, exit non-zero on failure
#
# Provenance:
#   Marker panels + promoter-overlap logic from
#     /data2/users/JCRLab/JBader/JBader_scHFD/02_analysis/S2_6b_filter_peak_atlas.R
#       (KEY_MARKERS, annotate_peaks_to_genes, promoter_window = c(-2000, 500))
#   Threshold (marker_retention_min = 0.95) from
#     /data2/users/JCRLab/JBader/JBader_scHFD/02_analysis/config/peak_filtration_config.yaml
#
# Gate: of the rare-marker promoter peaks present in the ORIGINAL atlas, the
# FILTERED atlas must retain >= min_ratio (default 0.95). Losing rare-marker
# promoters is exactly the failure the Primary+Rescue filter exists to prevent.
#
# Default marker panels (the source KEY_MARKERS, mouse symbols):
#   HSC        = Gata2, Procr, Hlf
#   pDC        = Tcf4, Siglech, Bst2
#   Neutrophil = S100a8, S100a9, Cebpe
# Override with a markers file (one symbol per line) via the 4th argument.
#
# Genome / annotation are parameterized via env vars:
#   MARKER_ENSDB  (default EnsDb.Mmusculus.v79)  -- gene coordinates source
#   Set them for your build (e.g. EnsDb.Hsapiens.v86 with human symbols).
#
# Usage:
#   Rscript check_marker_retention.R original.rds filtered.rds [min_ratio] [markers.txt]
#   # *.rds files contain GRanges peak atlases (or RangedSummarizedExperiment).

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript check_marker_retention.R <original> <filtered> [min_ratio] [markers.txt]")
}
orig_path <- args[[1]]
filt_path <- args[[2]]
min_ratio <- if (length(args) >= 3 && nzchar(args[[3]])) as.numeric(args[[3]]) else 0.95
markers_path <- if (length(args) >= 4 && nzchar(args[[4]])) args[[4]] else NULL

# Default rare-marker panels (source KEY_MARKERS).
DEFAULT_MARKERS <- list(
  HSC        = c("Gata2", "Procr", "Hlf"),
  pDC        = c("Tcf4", "Siglech", "Bst2"),
  Neutrophil = c("S100a8", "S100a9", "Cebpe")
)
PROMOTER_WINDOW <- c(-2000, 500)  # upstream/downstream of TSS
ENSDB_PKG <- Sys.getenv("MARKER_ENSDB", unset = "EnsDb.Mmusculus.v79")

if (!is.null(markers_path)) {
  genes <- readLines(markers_path)
  genes <- trimws(genes[nzchar(trimws(genes))])
  marker_sets <- list(markers = genes)
} else {
  marker_sets <- DEFAULT_MARKERS
}

load_peaks <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    obj <- readRDS(path)
    if (inherits(obj, "GRanges")) return(obj)
    if (inherits(obj, "RangedSummarizedExperiment")) return(SummarizedExperiment::rowRanges(obj))
    stop("RDS does not contain a GRanges / RangedSummarizedExperiment: ", path)
  }
  rtracklayer::import(path)
}

# Load the gene-coordinate EnsDb (build-parameterized).
if (!requireNamespace(ENSDB_PKG, quietly = TRUE)) {
  message("[FAIL] EnsDb package not installed: ", ENSDB_PKG,
          " (set MARKER_ENSDB to your build)")
  quit(status = 1)
}
suppressPackageStartupMessages(library(ENSDB_PKG, character.only = TRUE))
ensdb <- get(ENSDB_PKG)

#' Count peaks overlapping the promoter of any gene in `gene_names`.
#' Reproduces annotate_peaks_to_genes() (promoter_window = c(-2000, 500)).
marker_promoter_peaks <- function(peaks, gene_names) {
  g <- ensembldb::genes(ensdb)
  g <- keepStandardChromosomes(g, pruning.mode = "coarse")
  GenomeInfoDb::seqlevelsStyle(g) <- "UCSC"   # match chr-prefixed peak atlases
  g <- g[g$gene_name %in% gene_names]
  if (length(g) == 0) return(GRanges())
  proms <- promoters(g, upstream = -PROMOTER_WINDOW[1], downstream = PROMOTER_WINDOW[2])
  subsetByOverlaps(peaks, proms)
}

peaks_orig <- load_peaks(orig_path)
peaks_filt <- load_peaks(filt_path)

message("=== Rare-Marker Promoter Retention ===")
message(sprintf("Original atlas: %s peaks | Filtered atlas: %s peaks",
                format(length(peaks_orig), big.mark = ","),
                format(length(peaks_filt), big.mark = ",")))

all_markers <- unique(unlist(marker_sets))
n_orig <- length(marker_promoter_peaks(peaks_orig, all_markers))
n_filt <- length(marker_promoter_peaks(peaks_filt, all_markers))

if (n_orig == 0) {
  message("[FAIL] No marker promoter peaks in the ORIGINAL atlas - check marker symbols / genome build")
  quit(status = 1)
}
ratio <- n_filt / n_orig

issues <- character()
# Per-set check: a set dropping to zero is a hard failure.
for (set_name in names(marker_sets)) {
  o <- length(marker_promoter_peaks(peaks_orig, marker_sets[[set_name]]))
  f <- length(marker_promoter_peaks(peaks_filt, marker_sets[[set_name]]))
  r <- if (o > 0) f / o else NA_real_
  message(sprintf("  %-12s (%s): orig=%d, filt=%d (%s)",
                  set_name, paste(marker_sets[[set_name]], collapse = ", "), o, f,
                  if (is.na(r)) "n/a" else sprintf("%.1f%%", 100 * r)))
  if (!is.na(r) && f == 0) issues <- c(issues, paste0(set_name, ": ALL promoter peaks lost"))
}

message(sprintf("\nOverall retention: %.4f (threshold >= %.2f)", ratio, min_ratio))

if (ratio >= min_ratio && length(issues) == 0) {
  message("[PASS] Rare-marker promoter retention OK")
  quit(status = 0)
} else {
  if (ratio < min_ratio) issues <- c(issues, sprintf("overall %.4f below %.2f", ratio, min_ratio))
  message("[FAIL] ", paste(issues, collapse = "; "))
  quit(status = 1)
}
