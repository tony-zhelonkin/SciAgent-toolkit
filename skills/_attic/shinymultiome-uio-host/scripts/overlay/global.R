####################################
##   shinymultiome-uio-host overlay (global.R)
####################################
# Replaces upstream EskelandLab/ShinyMultiomeUiO/global.R hard-coded values:
#   - library(BSgenome.Hsapiens.UCSC.hg38)  →  library(${BSGENOME_PKG})
#   - seuratObject = "/Users/akshay/.../*.RDS"  →  Sys.getenv("RDS_PATH")
#   - fragFilePath = "/Users/akshay/.../*.tsv.gz"  →  Sys.getenv("FRAGMENTS_DIR")
# All other behaviour (variable name `pbmc`, fragment path UpdatePath wiring) is
# preserved verbatim so upstream server.R and ui.R work without patches.
#
# Required environment variables (set via docker-compose / .env):
#   RDS_PATH          absolute path to the prepared .rds inside the container
#                     (e.g. /data/rds/S2_explore_shinymultiome.rds)
#   FRAGMENTS_DIR     directory inside the container holding *.tsv.gz + .tbi
#                     (e.g. /data/fragments). Object's existing fragment paths
#                     should already be rewritten to this dir by Phase A.3 — but
#                     this overlay re-applies UpdatePath() defensively at startup.
#   BSGENOME_PKG      e.g. BSgenome.Mmusculus.UCSC.mm39
#   ATAC_ASSAY        e.g. peaks   (must match upstream server.R hard-codes)
#   RNA_ASSAY         e.g. RNA
#
# Optional:
#   DEFAULT_REGION    gene name pre-filled in Coverage tab (UI hint, not used here)
#   GROUP_BY_COLUMNS  comma-separated metadata columns exposed in dropdowns
####################################

####################################
##   Resolve env (fail fast)
####################################
.req <- function(key) {
  val <- Sys.getenv(key, unset = NA)
  if (is.na(val) || nchar(val) == 0) {
    stop(sprintf("[shinymultiome] required env var %s is unset", key), call. = FALSE)
  }
  val
}
RDS_PATH      <- .req("RDS_PATH")
FRAGMENTS_DIR <- .req("FRAGMENTS_DIR")
BSGENOME_PKG  <- .req("BSGENOME_PKG")
ATAC_ASSAY    <- Sys.getenv("ATAC_ASSAY", unset = "peaks")
RNA_ASSAY     <- Sys.getenv("RNA_ASSAY",  unset = "RNA")

if (!file.exists(RDS_PATH)) {
  stop(sprintf("[shinymultiome] RDS_PATH does not exist: %s", RDS_PATH), call. = FALSE)
}
if (!dir.exists(FRAGMENTS_DIR)) {
  stop(sprintf("[shinymultiome] FRAGMENTS_DIR does not exist: %s", FRAGMENTS_DIR), call. = FALSE)
}

####################################
##   Load libraries
####################################
suppressPackageStartupMessages({
  library(shiny)
  library(Seurat)
  library(Signac)
  library(patchwork)
  library(ggplot2)
  library(viridis)
  library(shinybusy)
  library(shinyBS)
})

# Genome-specific BSgenome (parameterised; was hard-coded hg38 upstream).
if (!requireNamespace(BSGENOME_PKG, quietly = TRUE)) {
  stop(sprintf("[shinymultiome] BSgenome package not installed in image: %s", BSGENOME_PKG),
       call. = FALSE)
}
suppressPackageStartupMessages(library(BSGENOME_PKG, character.only = TRUE))

####################################
##   Load Seurat object
####################################
message(sprintf("[shinymultiome] loading %s ...", RDS_PATH))
pbmc <- readRDS(RDS_PATH)

# Sanity checks — fail fast at app startup with a clear error rather than
# letting upstream server.R blow up later inside a renderPlot.
if (!ATAC_ASSAY %in% Assays(pbmc)) {
  stop(sprintf("[shinymultiome] ATAC_ASSAY '%s' not found in object. Available: %s",
               ATAC_ASSAY, paste(Assays(pbmc), collapse = ", ")), call. = FALSE)
}
if (!RNA_ASSAY %in% Assays(pbmc)) {
  warning(sprintf("[shinymultiome] RNA_ASSAY '%s' not found in object. Feature tab will be limited.",
                  RNA_ASSAY), call. = FALSE)
}

####################################
##   Defensive Fragment path rewrite
####################################
# Phase A.3 (prepare_for_shinymultiome) should already have rewritten Fragment
# paths to FRAGMENTS_DIR. We re-apply UpdatePath here so the image still works
# if someone hands us an unprepared .rds (slower path: discover .tsv.gz files
# in FRAGMENTS_DIR, match by basename).
DefaultAssay(pbmc) <- ATAC_ASSAY
frags <- Fragments(pbmc)
if (length(frags) == 0) {
  warning("[shinymultiome] object has no Fragments() — Coverage tab will be empty",
          call. = FALSE)
} else {
  for (i in seq_along(frags)) {
    cur_path <- GetFragmentData(frags[[i]], slot = "path")
    base     <- basename(cur_path)
    new_path <- file.path(FRAGMENTS_DIR, base)
    if (!file.exists(new_path)) {
      stop(sprintf("[shinymultiome] fragment file not found in FRAGMENTS_DIR: %s (looked for %s)",
                   base, new_path), call. = FALSE)
    }
    if (!file.exists(paste0(new_path, ".tbi"))) {
      stop(sprintf("[shinymultiome] fragment index missing: %s.tbi", new_path), call. = FALSE)
    }
    if (!identical(cur_path, new_path)) {
      message(sprintf("[shinymultiome] UpdatePath: %s -> %s", cur_path, new_path))
      frags[[i]] <- UpdatePath(frags[[i]], new.path = new_path)
    }
  }
  Fragments(pbmc) <- NULL
  Fragments(pbmc) <- frags
}

####################################
##   Annotation sanity (UCSC seqlevels)
####################################
# Phase A.2 should have set this; warn loudly if seqlevels look Ensembl-style
# (the #1 cause of empty CoveragePlots).
ann <- Annotation(pbmc[[ATAC_ASSAY]])
if (is.null(ann) || length(ann) == 0) {
  warning("[shinymultiome] Annotation(obj[[ATAC_ASSAY]]) is empty — CoveragePlot gene tracks will be missing",
          call. = FALSE)
} else {
  s1 <- as.character(seqnames(ann))[1]
  if (!startsWith(s1, "chr")) {
    warning(sprintf("[shinymultiome] Annotation seqlevels look Ensembl-style ('%s') — convert to UCSC ('chr%s') in Phase A.2",
                    s1, s1), call. = FALSE)
  }
}

message("[shinymultiome] global.R ready. Starting Shiny ...")

####################################
##   End of overlay
####################################
