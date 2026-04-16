# extract_gene_loadings.R -- CoReSh-to-GSEA bridge
#
# Given a top-ranking CORESH hit (GSE), compute gene-level loadings onto the
# query direction and emit a symbol-keyed gene set for downstream GSEA.
#
# Usage (inside R):
#   source("extract_gene_loadings.R")
#   gmt <- build_coresh_gmt(
#     top_hits = top_ranking,                 # data.table from coresh_batch.R
#     queries  = queries,                     # named list of integer Entrez queries
#     chunk_dir = "preprocessed_chunks/hsa/",
#     species  = "human",
#     n_top    = 50,
#     min_size = 15,
#     max_size = 500,
#     jaccard_threshold = 0.8
#   )
#   writeLines(gmt, "coresh_derived_sets.gmt")
#
# Requires scripts/symbols_to_entrez.R in the same directory.

suppressPackageStartupMessages({
  library(data.table)
  library(qs2)
})

# Caller must source("symbols_to_entrez.R") before sourcing this file --
# keeping the two scripts decoupled avoids brittle path detection.
if (!exists("sym2ent", mode = "function") || !exists("ent2sym", mode = "function")) {
  stop("source('symbols_to_entrez.R') before sourcing extract_gene_loadings.R")
}

.index_chunks <- function(chunk_dir) {
  # one-time index: GSE -> chunk path. Memoized via options().
  key <- paste0("coresh_index_", digest::digest(chunk_dir))
  idx <- getOption(key)
  if (!is.null(idx)) return(idx)
  chunk_paths <- list.files(chunk_dir, pattern = "full_objects\\.qs2$",
                            full.names = TRUE)
  idx <- rbindlist(lapply(chunk_paths, function(p) {
    ch <- qs_read(p)
    data.table(gse = vapply(ch, function(o) o$gseId, character(1)),
               chunk = p)
  }))
  do.call(options, setNames(list(idx), key))
  idx
}

extract_gene_loadings <- function(chunk_path, gse_id, query, n_top = 50) {
  ch <- qs_read(chunk_path)
  hit <- Filter(function(o) o$gseId == gse_id, ch)
  if (length(hit) == 0L) stop(sprintf("GSE %s not found in %s", gse_id, chunk_path))
  obj <- hit[[1]]

  query_idxs <- na.omit(match(query, obj$rownames))
  if (length(query_idxs) < 3L) stop(sprintf(
    "GSE %s: only %d/%d query genes present; aborting loading extraction.",
    gse_id, length(query_idxs), length(query)))

  E <- obj$E1024 / 1024
  profile <- colSums(E[query_idxs, , drop = FALSE])
  profile <- profile / sqrt(sum(profile^2))
  gene_loadings <- as.numeric(E %*% profile)

  ord <- order(abs(gene_loadings), decreasing = TRUE)
  keep <- head(ord, n_top)
  data.table(
    entrez  = obj$rownames[keep],
    loading = gene_loadings[keep],
    rank    = seq_along(keep)
  )
}

dedupe_jaccard <- function(sets, threshold = 0.8) {
  keep <- rep(TRUE, length(sets))
  for (i in seq_along(sets)) {
    if (!keep[i]) next
    for (j in seq_len(i - 1)) {
      if (!keep[j]) next
      a <- sets[[i]]$genes; b <- sets[[j]]$genes
      jacc <- length(intersect(a, b)) / length(union(a, b))
      if (jacc > threshold) { keep[i] <- FALSE; break }
    }
  }
  sets[keep]
}

build_coresh_gmt <- function(top_hits, queries, chunk_dir, species = "human",
                             n_top = 50, min_size = 15, max_size = 500,
                             jaccard_threshold = 0.8) {
  idx <- .index_chunks(chunk_dir)
  sets <- vector("list", nrow(top_hits))
  for (i in seq_len(nrow(top_hits))) {
    row <- top_hits[i]
    chunk_path <- idx[gse == row$gse]$chunk[1]
    if (is.na(chunk_path)) { message("skip: ", row$gse, " not indexed"); next }
    tryCatch({
      loadings <- extract_gene_loadings(chunk_path, row$gse,
                                        queries[[row$query_name]], n_top = n_top)
      symbols <- ent2sym(loadings$entrez, species = species)
      genes   <- unique(stats::na.omit(as.character(symbols)))
      if (length(genes) >= min_size && length(genes) <= max_size) {
        sets[[i]] <- list(
          name  = sprintf("CORESH_%s_%s", row$query_name, row$gse),
          genes = genes)
      }
    }, error = function(e) message("skip: ", row$gse, " -- ", conditionMessage(e)))
  }
  sets <- Filter(Negate(is.null), sets)
  sets <- dedupe_jaccard(sets, threshold = jaccard_threshold)

  vapply(sets, function(s) paste(c(s$name, "-", s$genes), collapse = "\t"),
         character(1))
}
