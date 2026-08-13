#!/usr/bin/env Rscript
####################################################################
##  validate_signac_rds.R — pre-deploy schema check
####################################################################
##  Run before `docker compose up` to fail fast if Phase A was skipped or
##  incomplete. Usage:
##
##      Rscript checks/validate_signac_rds.R 03_results/objects/<dataset>.rds \
##          [--atac-assay peaks] [--rna-assay RNA] \
##          [--fragments-host-dir /scratch/<project>/03_results/fragments] \
##          [--group-by cell_type,leiden_0.8,sample_id] \
##          [--require-links] [--require-motifs]
##
##  Exits 0 on success, 1 on any failure. Prints a per-check verdict so the
##  agent (and the user) can see exactly what failed.
####################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(GenomeInfoDb)
})

`%||%` <- function(a, b) if (is.null(a) || (is.character(a) && nchar(a) == 0)) b else a

.parse_args <- function(args) {
  if (length(args) < 1L) stop("Usage: validate_signac_rds.R <path> [--key val ...]")
  rds <- args[1]
  pairs <- list()
  i <- 2
  while (i <= length(args)) {
    k <- args[i]
    if (!startsWith(k, "--")) stop(sprintf("Unexpected arg: %s", k))
    is_flag <- (i + 1) > length(args) || startsWith(args[i + 1], "--")
    v <- if (is_flag) "TRUE" else args[i + 1]
    pairs[[sub("^--", "", k)]] <- v
    i <- if (is_flag) i + 1 else i + 2
  }
  list(rds = rds, opts = pairs)
}

.report <- function(name, ok, detail = "") {
  tag <- if (ok) "OK  " else "FAIL"
  suffix <- if (nchar(detail) > 0) sprintf(" — %s", detail) else ""
  cat(sprintf("  [%s] %s%s\n", tag, name, suffix))
  ok
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  parsed <- .parse_args(args)
  rds <- parsed$rds
  opts <- parsed$opts

  atac_assay <- opts$`atac-assay` %||% "peaks"
  rna_assay  <- opts$`rna-assay` %||% "RNA"
  frag_host_dir <- opts$`fragments-host-dir` %||% NULL
  group_by <- strsplit(opts$`group-by` %||% "cell_type,leiden_0.8,sample_id", ",")[[1]]
  require_links  <- isTRUE(opts$`require-links`  == "TRUE")
  require_motifs <- isTRUE(opts$`require-motifs` == "TRUE")

  if (!file.exists(rds)) {
    cat(sprintf("  [FAIL] file does not exist: %s\n", rds))
    quit(status = 1)
  }

  cat(sprintf("Validating Signac multiome schema for %s\n", rds))
  obj <- tryCatch(readRDS(rds), error = function(e) {
    cat(sprintf("  [FAIL] could not readRDS: %s\n", conditionMessage(e)))
    quit(status = 1)
  })
  cat(sprintf("  shape: %d cells, assays = [%s]\n",
              ncol(obj), paste(Assays(obj), collapse = ", ")))

  results <- c()

  # 1. assay names
  results <- c(results, .report(
    sprintf("ATAC assay '%s' present", atac_assay),
    atac_assay %in% Assays(obj),
    sprintf("got [%s]", paste(Assays(obj), collapse = ", "))
  ))
  results <- c(results, .report(
    sprintf("RNA assay '%s' present", rna_assay),
    rna_assay %in% Assays(obj)
  ))

  # 2. ChromatinAssay class
  if (atac_assay %in% Assays(obj)) {
    cls <- class(obj[[atac_assay]])[1]
    results <- c(results, .report(
      sprintf("'%s' is a ChromatinAssay", atac_assay),
      cls == "ChromatinAssay",
      sprintf("got %s", cls)
    ))
  }

  # 3. Annotation present + UCSC seqlevels
  ann_ok_present <- FALSE
  ann_ok_ucsc    <- FALSE
  if (atac_assay %in% Assays(obj) && class(obj[[atac_assay]])[1] == "ChromatinAssay") {
    ann <- Annotation(obj[[atac_assay]])
    ann_ok_present <- !is.null(ann) && length(ann) > 0
    if (ann_ok_present) {
      s1 <- as.character(seqnames(ann))[1]
      ann_ok_ucsc <- startsWith(s1, "chr")
    }
  }
  results <- c(results, .report("Annotation() present and non-empty", ann_ok_present))
  results <- c(results, .report("Annotation seqlevels are UCSC-styled", ann_ok_ucsc))

  # 4. Fragments accessible
  frag_results <- c()
  if (atac_assay %in% Assays(obj) && class(obj[[atac_assay]])[1] == "ChromatinAssay") {
    frags <- Fragments(obj[[atac_assay]])
    if (length(frags) == 0) {
      frag_results <- c(frag_results, .report("Fragments() non-empty", FALSE,
                                              "no Fragments objects on assay"))
    } else {
      for (i in seq_along(frags)) {
        cur <- GetFragmentData(frags[[i]], slot = "path")
        # If --fragments-host-dir given, validate by basename inside that dir
        # (the path stored in the .rds is a CONTAINER path; verifier runs on host)
        if (!is.null(frag_host_dir)) {
          base <- basename(cur)
          host_path <- file.path(frag_host_dir, base)
        } else {
          host_path <- cur  # validator running inside container or host==container
        }
        frag_results <- c(frag_results, .report(
          sprintf("fragment file accessible: %s", basename(cur)),
          file.exists(host_path),
          sprintf("checked %s", host_path)
        ))
        frag_results <- c(frag_results, .report(
          sprintf("fragment .tbi accessible: %s.tbi", basename(cur)),
          file.exists(paste0(host_path, ".tbi"))
        ))
      }
    }
  }
  results <- c(results, frag_results)

  # 5. group.by columns
  for (col in group_by) {
    results <- c(results, .report(
      sprintf("meta.data['%s'] present", col),
      col %in% colnames(obj@meta.data),
      sprintf("%d distinct values",
              if (col %in% colnames(obj@meta.data)) length(unique(obj@meta.data[[col]])) else 0)
    ))
  }

  # 6. Optional: Links
  if (require_links && atac_assay %in% Assays(obj)) {
    n <- tryCatch(length(Links(obj[[atac_assay]])), error = function(e) 0)
    results <- c(results, .report("Links() populated (--require-links)",
                                  n > 0, sprintf("%d links", n)))
  }

  # 7. Optional: Motifs + footprint enrichment
  if (require_motifs && atac_assay %in% Assays(obj)) {
    has_motifs <- tryCatch(
      !is.null(obj@assays[[atac_assay]]@motifs) &&
        length(obj@assays[[atac_assay]]@motifs@motif.names) > 0,
      error = function(e) FALSE
    )
    has_pe <- tryCatch(
      length(obj@assays[[atac_assay]]@positionEnrichment) > 0,
      error = function(e) FALSE
    )
    results <- c(results, .report("Motif slot populated (--require-motifs)", has_motifs))
    results <- c(results, .report("positionEnrichment populated (--require-motifs)", has_pe))
  }

  n_pass <- sum(results)
  n_total <- length(results)
  cat(sprintf("\n  %d/%d checks passed\n", n_pass, n_total))
  if (n_pass != n_total) quit(status = 1)
}

main()
