#!/usr/bin/env Rscript
####################################################################
##  prepare_for_shinymultiome.R — Phase A object preflight
####################################################################
##  Reads a Signac multiome Seurat .rds, normalises it for ShinyMultiome.UiO
##  hosting (assay names, Annotation seqlevels, Fragment paths, group.by
##  columns), optionally computes LinkPeaks and motif footprints, and writes
##  a deploy-ready .rds.
##
##  Usage from R:
##    source("scripts/prepare_for_shinymultiome.R")
##    obj <- prepare_for_shinymultiome(
##      in_rds       = "03_results/checkpoints/<input>.rds",
##      out_rds      = "03_results/objects/<output>.rds",
##      genome_label = "mm39",
##      ensdb_pkg    = "EnsDb.Mmusculus.v107",
##      atac_assay   = "peaks",
##      rna_assay    = "RNA",
##      fragments_container_dir = "/data/fragments",
##      group_by_columns = c("cell_type", "leiden_0.8", "sample_id"),
##      compute_links = TRUE,
##      compute_motif_footprints = FALSE
##    )
##
##  Usage from CLI:
##    Rscript scripts/prepare_for_shinymultiome.R \
##      --in <in.rds> --out <out.rds> --genome mm39 \
##      --ensdb EnsDb.Mmusculus.v107 --atac-assay peaks --rna-assay RNA \
##      --fragments-container-dir /data/fragments \
##      --group-by cell_type,leiden_0.8,sample_id
####################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(GenomeInfoDb)
})

# Tiny null-coalescing helper. Defined before the CLI block so the latter can
# use it for default-arg fallbacks.
`%||%` <- function(a, b) if (is.null(a)) b else a

# ------------------------------------------------------------------ #
# A.1 — Assay name normalisation
# ------------------------------------------------------------------ #
normalize_assay_names <- function(obj, rna_assay = "RNA", atac_assay = "peaks",
                                  rename_strategy = c("rename", "copy", "skip")) {
  rename_strategy <- match.arg(rename_strategy)
  current_assays <- Assays(obj)
  message(sprintf("[A.1] Current assays: %s", paste(current_assays, collapse = ", ")))

  detect_atac <- function(obj) {
    cand <- intersect(c("peaks", "ATAC", "atac", "Peaks"), Assays(obj))
    if (length(cand) == 0) {
      # last resort: any ChromatinAssay
      cls <- vapply(Assays(obj), function(a) class(obj[[a]])[1], character(1))
      cand <- names(cls)[cls == "ChromatinAssay"]
    }
    if (length(cand) == 0) {
      stop("[A.1] Could not detect a ChromatinAssay in object")
    }
    cand[1]
  }
  detect_rna <- function(obj) {
    cand <- intersect(c("RNA", "rna", "Gene Expression", "SCT"), Assays(obj))
    if (length(cand) == 0) {
      cls <- vapply(Assays(obj), function(a) class(obj[[a]])[1], character(1))
      cand <- names(cls)[cls %in% c("Assay", "Assay5")]
    }
    if (length(cand) == 0) {
      warning("[A.1] Could not detect an RNA assay; Feature tab will be limited")
      return(NA_character_)
    }
    cand[1]
  }

  src_atac <- detect_atac(obj)
  src_rna  <- detect_rna(obj)
  message(sprintf("[A.1] Detected — ATAC: '%s', RNA: '%s'", src_atac, src_rna))

  if (rename_strategy == "skip") {
    message("[A.1] Skipping rename (--rename-strategy skip)")
    return(obj)
  }

  do_rename <- function(obj, src, dst) {
    if (identical(src, dst) || is.na(src)) return(obj)
    if (dst %in% Assays(obj)) {
      stop(sprintf("[A.1] Target assay name '%s' already exists; would clobber. Pick a different name or use copy strategy.",
                   dst))
    }
    if (rename_strategy == "rename") {
      message(sprintf("[A.1] RenameAssays: %s -> %s", src, dst))
      obj <- RenameAssays(obj, assay.name = src, new.assay.name = dst)
    } else if (rename_strategy == "copy") {
      message(sprintf("[A.1] Copying assay: %s -> %s (original kept)", src, dst))
      obj[[dst]] <- obj[[src]]
    }
    obj
  }

  obj <- do_rename(obj, src_atac, atac_assay)
  if (!is.na(src_rna)) obj <- do_rename(obj, src_rna, rna_assay)

  DefaultAssay(obj) <- atac_assay
  message(sprintf("[A.1] DefaultAssay set to '%s'", atac_assay))
  obj
}

# ------------------------------------------------------------------ #
# A.2 — Annotation seqlevels
# ------------------------------------------------------------------ #
ensure_signac_annotation <- function(obj, atac_assay = "peaks",
                                     ensdb_pkg = "EnsDb.Mmusculus.v107",
                                     genome    = "mm39",
                                     seqlevels_style = "UCSC") {
  ann <- Annotation(obj[[atac_assay]])
  needs_build <- is.null(ann) || length(ann) == 0
  if (!needs_build) {
    s1 <- as.character(seqnames(ann))[1]
    needs_restyle <- !startsWith(s1, "chr") && seqlevels_style == "UCSC"
    if (needs_restyle) {
      message(sprintf("[A.2] Annotation seqlevels style is %s; converting to %s",
                      seqlevelsStyle(ann)[1], seqlevels_style))
      seqlevelsStyle(ann) <- seqlevels_style
      genome(ann) <- genome
      Annotation(obj[[atac_assay]]) <- ann
    } else {
      message("[A.2] Annotation present and seqlevels OK — leaving alone")
    }
    return(obj)
  }
  message(sprintf("[A.2] Annotation missing — building from %s", ensdb_pkg))
  if (!requireNamespace(ensdb_pkg, quietly = TRUE)) {
    stop(sprintf("[A.2] EnsDb package not installed: %s", ensdb_pkg))
  }
  ensdb <- get(ensdb_pkg, envir = asNamespace(ensdb_pkg))
  gr <- Signac::GetGRangesFromEnsDb(ensdb)
  if (seqlevels_style == "UCSC") seqlevelsStyle(gr) <- "UCSC"
  genome(gr) <- genome
  Annotation(obj[[atac_assay]]) <- gr
  message(sprintf("[A.2] Annotation built (%d ranges, seqlevels[1]=%s, genome=%s)",
                  length(gr), as.character(seqnames(gr))[1], genome))
  obj
}

# ------------------------------------------------------------------ #
# A.3 — Fragment file paths
# ------------------------------------------------------------------ #
rewrite_fragment_paths <- function(obj, atac_assay = "peaks",
                                   container_dir = "/data/fragments",
                                   host_dir      = NULL) {
  # `container_dir` is the path INSIDE the container (what gets written into
  # the .rds for the deployed app). `host_dir` is the path on the host machine
  # where the same files live now (what we use to verify existence). If
  # host_dir is NULL the function uses container_dir for verification too —
  # only correct when running on the deploy machine.
  verify_dir <- if (is.null(host_dir)) container_dir else host_dir
  frags <- Fragments(obj[[atac_assay]])
  if (length(frags) == 0) {
    warning("[A.3] No Fragments() on object; skipping rewrite")
    return(obj)
  }
  for (i in seq_along(frags)) {
    cur <- GetFragmentData(frags[[i]], slot = "path")
    base <- basename(cur)
    new_container <- file.path(container_dir, base)
    new_verify    <- file.path(verify_dir,    base)
    if (!file.exists(new_verify)) {
      stop(sprintf("[A.3] Fragment file missing on host (%s). Move it to %s before deploy.",
                   new_verify, verify_dir))
    }
    if (!file.exists(paste0(new_verify, ".tbi"))) {
      stop(sprintf("[A.3] Fragment index missing on host: %s.tbi", new_verify))
    }
    if (!identical(cur, new_container)) {
      message(sprintf("[A.3] UpdatePath: %s -> %s", cur, new_container))
      frags[[i]] <- UpdatePath(frags[[i]], new.path = new_container)
    }
  }
  Fragments(obj[[atac_assay]]) <- NULL
  Fragments(obj[[atac_assay]]) <- frags
  obj
}

# ------------------------------------------------------------------ #
# A.4 — group.by column validation
# ------------------------------------------------------------------ #
validate_group_columns <- function(obj, required, max_values_per_column = 60) {
  md <- obj@meta.data
  missing <- setdiff(required, colnames(md))
  if (length(missing) > 0) {
    stop(sprintf("[A.4] Required group.by columns missing: %s",
                 paste(missing, collapse = ", ")))
  }
  for (col in required) {
    n <- length(unique(md[[col]]))
    msg <- sprintf("[A.4]   %s: %d distinct values", col, n)
    if (n > max_values_per_column) {
      warning(sprintf("%s — exceeds max_values_per_column=%d. Coverage plots will be unreadable; consider coarsening.",
                      msg, max_values_per_column))
    } else {
      message(msg)
    }
  }
  invisible(TRUE)
}

# ------------------------------------------------------------------ #
# A.5 — Optional: LinkPeaks
# ------------------------------------------------------------------ #
maybe_compute_links <- function(obj, atac_assay = "peaks", rna_assay = "RNA",
                                compute_links = FALSE, genes.use = NULL) {
  if (!compute_links) {
    message("[A.5] compute_links=FALSE — skipping")
    return(obj)
  }
  existing <- tryCatch(length(Links(obj[[atac_assay]])), error = function(e) 0)
  if (existing > 0) {
    message(sprintf("[A.5] Links already populated (%d) — skipping", existing))
    return(obj)
  }
  if (!rna_assay %in% Assays(obj)) {
    warning(sprintf("[A.5] RNA assay '%s' not present — cannot compute Links. Skipping.",
                    rna_assay))
    return(obj)
  }
  message(sprintf("[A.5] Running Signac::LinkPeaks (atac=%s, rna=%s, genes=%s)",
                  atac_assay, rna_assay,
                  if (is.null(genes.use)) "NULL (all)" else sprintf("%d", length(genes.use))))
  obj <- Signac::LinkPeaks(
    obj, peak.assay = atac_assay, expression.assay = rna_assay, genes.use = genes.use
  )
  message(sprintf("[A.5] Done — %d Links populated", length(Links(obj[[atac_assay]]))))
  obj
}

# ------------------------------------------------------------------ #
# A.6 — Optional: motif scan + footprint precompute
# ------------------------------------------------------------------ #
maybe_compute_footprints <- function(obj, atac_assay = "peaks",
                                     compute_motif_footprints = FALSE,
                                     bsgenome_pkg  = "BSgenome.Mmusculus.UCSC.mm39",
                                     pfm = NULL,
                                     motif_names = NULL,
                                     group.by = NULL) {
  if (!compute_motif_footprints) {
    message("[A.6] compute_motif_footprints=FALSE — skipping")
    return(obj)
  }
  if (!requireNamespace(bsgenome_pkg, quietly = TRUE)) {
    stop(sprintf("[A.6] BSgenome package not installed: %s", bsgenome_pkg))
  }
  bsg <- get(bsgenome_pkg, envir = asNamespace(bsgenome_pkg))
  if (is.null(pfm)) {
    if (!requireNamespace("JASPAR2020", quietly = TRUE) ||
        !requireNamespace("TFBSTools", quietly = TRUE)) {
      stop("[A.6] pfm=NULL requires JASPAR2020 + TFBSTools installed")
    }
    pfm <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020,
                                   list(species = 9606, all_versions = FALSE))
  }
  message(sprintf("[A.6] AddMotifs (assay=%s, %d PFMs)", atac_assay, length(pfm)))
  obj <- Signac::AddMotifs(obj, genome = bsg, pfm = pfm, assay = atac_assay)
  if (is.null(motif_names)) {
    motif_names <- names(pfm)[1:min(length(pfm), 20)]
    message(sprintf("[A.6] motif_names=NULL — using first %d", length(motif_names)))
  }
  message(sprintf("[A.6] Footprint (group.by=%s, %d motifs)",
                  if (is.null(group.by)) "<default Idents>" else group.by,
                  length(motif_names)))
  obj <- Signac::Footprint(obj, motif.name = motif_names,
                           genome = bsg, group.by = group.by, assay = atac_assay)
  obj
}

# ------------------------------------------------------------------ #
# Top-level driver
# ------------------------------------------------------------------ #
prepare_for_shinymultiome <- function(
    in_rds, out_rds,
    genome_label,
    ensdb_pkg,
    atac_assay = "peaks",
    rna_assay  = "RNA",
    rename_strategy = c("rename", "copy", "skip"),
    fragments_container_dir = "/data/fragments",
    fragments_host_dir      = NULL,
    group_by_columns        = c("cell_type", "leiden_0.8", "sample_id"),
    max_values_per_column   = 60,
    compute_links           = TRUE,
    compute_motif_footprints = FALSE,
    bsgenome_pkg            = NULL,
    pfm                     = NULL,
    motif_names             = NULL,
    footprint_group_by      = NULL
) {
  rename_strategy <- match.arg(rename_strategy)
  if (!file.exists(in_rds)) stop(sprintf("in_rds does not exist: %s", in_rds))
  message(sprintf("[A] Reading %s", in_rds))
  obj <- readRDS(in_rds)

  obj <- normalize_assay_names(obj, rna_assay = rna_assay, atac_assay = atac_assay,
                               rename_strategy = rename_strategy)
  obj <- ensure_signac_annotation(obj, atac_assay = atac_assay,
                                  ensdb_pkg = ensdb_pkg,
                                  genome    = genome_label,
                                  seqlevels_style = "UCSC")
  obj <- rewrite_fragment_paths(obj, atac_assay = atac_assay,
                                container_dir = fragments_container_dir,
                                host_dir      = fragments_host_dir)
  validate_group_columns(obj, required = group_by_columns,
                         max_values_per_column = max_values_per_column)
  obj <- maybe_compute_links(obj, atac_assay = atac_assay, rna_assay = rna_assay,
                             compute_links = compute_links)
  if (compute_motif_footprints) {
    if (is.null(bsgenome_pkg)) stop("compute_motif_footprints=TRUE requires bsgenome_pkg")
    obj <- maybe_compute_footprints(obj, atac_assay = atac_assay,
                                    compute_motif_footprints = TRUE,
                                    bsgenome_pkg = bsgenome_pkg,
                                    pfm = pfm, motif_names = motif_names,
                                    group.by = footprint_group_by)
  }

  if (!dir.exists(dirname(out_rds))) dir.create(dirname(out_rds), recursive = TRUE)
  message(sprintf("[A] Saving %s", out_rds))
  saveRDS(obj, out_rds)
  message("[A] Done. Run: Rscript checks/validate_signac_rds.R ", out_rds)
  invisible(obj)
}

# ------------------------------------------------------------------ #
# CLI entry point
# ------------------------------------------------------------------ #
.parse_args <- function(args) {
  pairs <- list()
  i <- 1
  while (i <= length(args)) {
    k <- args[i]
    if (!startsWith(k, "--")) stop(sprintf("Unexpected arg: %s", k))
    v <- if ((i + 1) <= length(args) && !startsWith(args[i + 1], "--")) args[i + 1] else "TRUE"
    pairs[[sub("^--", "", k)]] <- v
    i <- if (v == "TRUE" && (i + 1) > length(args)) i + 1 else i + 2
  }
  pairs
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0L) {
    cat("Usage: Rscript prepare_for_shinymultiome.R --in <in.rds> --out <out.rds> ",
        "--genome <label> --ensdb <pkg> [--atac-assay peaks] [--rna-assay RNA] ",
        "[--fragments-container-dir /data/fragments] [--fragments-host-dir <path>] ",
        "[--group-by col1,col2,col3] [--no-links] [--footprints] [--bsgenome <pkg>]\n",
        sep = "")
    quit(status = 1)
  }
  p <- .parse_args(args)
  prepare_for_shinymultiome(
    in_rds       = p$`in`,
    out_rds      = p$out,
    genome_label = p$genome,
    ensdb_pkg    = p$ensdb,
    atac_assay   = p$`atac-assay` %||% "peaks",
    rna_assay    = p$`rna-assay` %||% "RNA",
    fragments_container_dir = p$`fragments-container-dir` %||% "/data/fragments",
    fragments_host_dir      = p$`fragments-host-dir`,
    group_by_columns        = strsplit(p$`group-by` %||% "cell_type,leiden_0.8,sample_id", ",")[[1]],
    compute_links           = !isTRUE(p$`no-links` == "TRUE"),
    compute_motif_footprints = isTRUE(p$footprints == "TRUE"),
    bsgenome_pkg            = p$bsgenome
  )
}
