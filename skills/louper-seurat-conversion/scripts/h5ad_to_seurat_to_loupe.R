#!/usr/bin/env Rscript
# h5ad_to_seurat_to_loupe.R
# Generalized reference script: AnnData h5ad → Seurat RDS → Loupe Browser .cloupe
#
# Tested against: loupeR 1.1.5 · anndataR 1.1.0 · Seurat 5.2.1 · R 4.5.0
#
# Key behaviours:
#   - Barcode auto-detection: tries existing barcodes first; falls back to
#     cellid_XXXXXXXXX format if they fail loupeR's whitelist check.
#   - Cell ID traceability: original IDs always preserved as a "cell_id" cluster
#     group in the Loupe file, regardless of barcode strategy.
#   - Uses count_mat (not counts) — loupeR v1.1.5 argument name.
#   - validate_barcodes() return value is checked as list$success, not boolean.

suppressPackageStartupMessages({
  library(anndataR)
  library(Seurat)
  library(loupeR)
})

# ── Configuration ─────────────────────────────────────────────────────────────
# Edit these paths and parameters for your project.

H5AD_IN    <- "path/to/input.h5ad"           # Input AnnData file
SEURAT_OUT <- "path/to/output_seurat.rds"    # Seurat RDS output (optional save)
LOUPE_DIR  <- "path/to/loupe_output/"        # Directory for .cloupe file
LOUPE_NAME <- "export"                        # Output: <LOUPE_NAME>.cloupe

# Layer names as they appear in AnnData (adata$layers)
COUNTS_LAYER <- "counts"            # Raw integer UMI layer (required for loupeR)
NORM_LAYER   <- "scvi_normalized"   # Optional normalized layer (saved to Seurat)

# obsm key for UMAP coordinates
UMAP_OBSM_KEY <- "X_umap"

# Cluster metadata columns from adata$obs to include in Loupe Browser.
# Columns not present in the data are skipped silently.
CLUSTER_COLS <- c("leiden", "cell_type", "tissue", "is_tumor", "batch")

# ── Setup ─────────────────────────────────────────────────────────────────────

stopifnot(file.exists(H5AD_IN))
dir.create(LOUPE_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== h5ad → Seurat → Loupe ===\n\n")

# ── Step 1: Read h5ad ─────────────────────────────────────────────────────────

cat("Step 1: Reading h5ad ...\n")
adata <- read_h5ad(H5AD_IN)

cat("  Cells x genes:", dim(adata)[1], "x", dim(adata)[2], "\n")
cat("  Layers:", paste(names(adata$layers), collapse = ", "), "\n")
cat("  obsm   :", paste(names(adata$obsm),  collapse = ", "), "\n")
cat("  obs cols (first 10):", paste(head(names(adata$obs), 10), collapse = ", "), "\n\n")

# Verify required layer exists
if (!COUNTS_LAYER %in% names(adata$layers)) {
  stop("Required layer '", COUNTS_LAYER, "' not found.\n",
       "Available layers: ", paste(names(adata$layers), collapse = ", "))
}

# ── Step 2: Convert to Seurat ─────────────────────────────────────────────────
# anndataR note: pass layer NAMES only (not "X"). X is included automatically
# as an extra layer called "X". All obsm entries are auto-converted to
# reductions, including internal scVI metadata keys — expect extra reductions.

layers_to_map <- c(COUNTS_LAYER)
if (!is.null(NORM_LAYER) && NORM_LAYER %in% names(adata$layers)) {
  layers_to_map <- c(layers_to_map, NORM_LAYER)
}

cat("Step 2: Converting to Seurat (layers:", paste(layers_to_map, collapse = ", "), ") ...\n")
seurat <- adata$as_Seurat(layers_mapping = layers_to_map)

cat("  Seurat layers   :", paste(Layers(seurat), collapse = ", "), "\n")
cat("  Reductions (obsm):", paste(Reductions(seurat), collapse = ", "), "\n")
cat("  Meta columns     :", paste(head(colnames(seurat@meta.data), 10), collapse = ", "), "\n\n")

# ── Step 3: Attach clean UMAP reduction ───────────────────────────────────────
# anndataR auto-generates reductions from obsm keys, but the names may not
# match the convention expected by downstream tools. Build a clean "umap" key.

cat("Step 3: Attaching UMAP ...\n")

if (!UMAP_OBSM_KEY %in% names(adata$obsm)) {
  stop("UMAP key '", UMAP_OBSM_KEY, "' not found in obsm.\n",
       "Available keys: ", paste(names(adata$obsm), collapse = ", "))
}

if (!"umap" %in% Reductions(seurat)) {
  umap_mat <- as.matrix(adata$obsm[[UMAP_OBSM_KEY]])
  rownames(umap_mat) <- colnames(seurat)
  colnames(umap_mat) <- c("UMAP_1", "UMAP_2")
  seurat[["umap"]] <- CreateDimReducObject(
    embeddings = umap_mat, key = "UMAP_", assay = "RNA"
  )
  cat("  UMAP added (2D)\n\n")
} else {
  cat("  UMAP already present\n\n")
}

# ── Step 4: Verify layer content ──────────────────────────────────────────────

cat("Step 4: Verifying layer content ...\n")
counts_check <- LayerData(seurat, COUNTS_LAYER)[1:3, 1:3]
cat("  First 3x3 of '", COUNTS_LAYER, "' (expect integers near 0):\n", sep = "")
print(round(counts_check, 2))
cat("\n")

# ── Step 5: Save Seurat RDS ───────────────────────────────────────────────────

if (!is.null(SEURAT_OUT) && nchar(SEURAT_OUT) > 0) {
  cat("Step 5: Saving Seurat RDS ...\n")
  dir.create(dirname(SEURAT_OUT), showWarnings = FALSE, recursive = TRUE)
  saveRDS(seurat, SEURAT_OUT)
  cat("  Saved to:", SEURAT_OUT, "\n\n")
} else {
  cat("Step 5: Skipping Seurat RDS save (SEURAT_OUT not set)\n\n")
}

# ── Step 6: Barcode strategy ───────────────────────────────────────────────────
# Try existing cell names first. If they fail loupeR's whitelist, fall back to
# cellid_XXXXXXXXX format (Loupe Browser's format for custom/non-10x data).
# Either way, preserve original IDs as a "cell_id" cluster group.

cat("Step 6: Determining barcode strategy ...\n")

original_cell_ids <- colnames(seurat)
n_cells           <- length(original_cell_ids)

cat("  Cell count:", n_cells, "\n")
cat("  Original ID example:", original_cell_ids[1], "\n")

check_original <- validate_barcodes(original_cell_ids)

if (isTRUE(check_original$success)) {
  # Existing barcodes are in a format loupeR accepts — use them directly
  loupe_barcodes   <- original_cell_ids
  barcodes_remapped <- FALSE
  cat("  Barcodes: VALID (using original cell names)\n\n")
} else {
  # Fall back to Loupe-native cellid format
  # cellid_XXXXXXXXX: 9-digit zero-padded integer, accepted by the louper binary
  cat("  Original barcodes not recognised:", check_original$msg, "\n")
  cat("  Falling back to cellid_XXXXXXXXX format ...\n")
  loupe_barcodes   <- sprintf("cellid_%09d", seq_len(n_cells))
  barcodes_remapped <- TRUE

  check_cellid <- validate_barcodes(loupe_barcodes)
  if (!isTRUE(check_cellid$success)) {
    stop("cellid barcode generation failed: ", check_cellid$msg)
  }
  cat("  Barcodes: VALID (cellid format)\n")
  cat("  Example:", loupe_barcodes[1], "\n\n")
}

# ── Step 7: Build cluster groups ──────────────────────────────────────────────

cat("Step 7: Building cluster groups ...\n")

clusters_for_loupe <- list()

for (col in CLUSTER_COLS) {
  if (col %in% colnames(seurat@meta.data)) {
    clusters_for_loupe[[col]] <- factor(seurat@meta.data[[col]])
  } else {
    cat("  Skipping '", col, "' (not in metadata)\n", sep = "")
  }
}

# Always add original cell IDs for traceability — even if barcodes weren't remapped.
# This appears as the "cell_id" category in Loupe Browser's cell info panel.
clusters_for_loupe[["cell_id"]] <- factor(original_cell_ids)

# Rename factor elements to match Loupe barcodes
clusters_for_loupe <- lapply(clusters_for_loupe, function(f) {
  names(f) <- loupe_barcodes
  f
})

cat("  Cluster groups included:\n")
for (nm in names(clusters_for_loupe)) {
  cat("    ", nm, "—", nlevels(clusters_for_loupe[[nm]]), "levels\n")
}
cat("\n")

# ── Step 8: Build projection ───────────────────────────────────────────────────

cat("Step 8: Building UMAP projection ...\n")

umap_embedding <- seurat[["umap"]]@cell.embeddings[, 1:2]
rownames(umap_embedding) <- loupe_barcodes
projections_for_loupe <- list(UMAP = umap_embedding)

cat("  Projection dimensions:", nrow(umap_embedding), "cells x 2\n\n")

# ── Step 9: Build count matrix ────────────────────────────────────────────────

cat("Step 9: Building count matrix ...\n")

counts_mat <- LayerData(seurat, COUNTS_LAYER)
colnames(counts_mat) <- loupe_barcodes

cat("  Matrix:", nrow(counts_mat), "genes x", ncol(counts_mat), "cells\n")
cat("  Output:", file.path(LOUPE_DIR, paste0(LOUPE_NAME, ".cloupe")), "\n\n")

# ── Step 10: Create Loupe Browser file ────────────────────────────────────────
# NOTE: argument is count_mat (not counts) — verified in loupeR v1.1.5.
# NOTE: force=TRUE only skips output-file overwrite; does NOT bypass the louper
#       binary's barcode whitelist check.

cat("Step 10: Creating Loupe Browser file ...\n")

create_loupe(
  count_mat   = counts_mat,
  clusters    = clusters_for_loupe,
  projections = projections_for_loupe,
  output_dir  = LOUPE_DIR,
  output_name = LOUPE_NAME
)

# ── Done ──────────────────────────────────────────────────────────────────────

loupe_file <- file.path(LOUPE_DIR, paste0(LOUPE_NAME, ".cloupe"))

if (file.exists(loupe_file)) {
  size_mb <- round(file.size(loupe_file) / 1e6, 1)
  cat("\n=== DONE ===\n")
  cat("Loupe file:", loupe_file, "(", size_mb, "MB)\n")
  if (!is.null(SEURAT_OUT) && file.exists(SEURAT_OUT)) {
    cat("Seurat RDS:", SEURAT_OUT, "\n")
  }
  if (barcodes_remapped) {
    cat("\nNote: barcodes in the .cloupe file use Loupe cellid format (not original IDs).\n")
    cat("Use the 'cell_id' cluster group in Loupe Browser to see original cell IDs.\n")
  }
} else {
  cat("\nWARNING: Expected .cloupe file not found at:", loupe_file, "\n")
  cat("Check for errors above or run create_bugreport() for diagnostics.\n")
}
