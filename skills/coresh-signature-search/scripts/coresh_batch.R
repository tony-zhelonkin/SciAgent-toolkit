# coresh_batch.R -- run multiple CORESH queries across all chunks
#
# Usage (inside R):
#   source("coresh_batch.R")
#   queries <- list(
#     iron_uptake  = c(7037L, 4891L, 55240L),
#     iron_storage = c(2495L, 2512L)
#   )
#   results <- coresh_batch(queries,
#                           chunk_dir = "preprocessed_chunks/hsa/",
#                           n_cores = 8, pvalues = FALSE)
#   data.table::fwrite(results, "coresh_results.csv")
#
# Returns a tidy data.table with columns:
#   query_name, gse, gpl, pctVar, pval, size, rank
#
# `rank` is the 1-based rank within each query_name (by -pctVar if pvalues=FALSE,
# by pval if pvalues=TRUE).
#
# Pre-flight checks:
#   * query is integer vector of length >= 3
#   * chunk_dir exists and contains *_full_objects.qs2 files
#   * first chunk matches expected structure ($E1024, $rownames, $gseId, $totalVar)

suppressPackageStartupMessages({
  library(data.table)
  library(qs2)
  library(BiocParallel)
  library(fgsea)
})

coreshMatch <- function(obj, query, calculatePvalues = FALSE) {
  E <- obj$E1024 / 1024
  queryIdxs <- na.omit(match(query, obj$rownames))
  k <- length(queryIdxs)
  if (k == 0L) {
    return(data.table(gse = obj$gseId, gpl = obj$gplId,
                      pctVar = 0, pval = NA_real_, size = 0L))
  }
  curProfile <- colSums(E[queryIdxs, , drop = FALSE])
  queryVar <- sum(curProfile^2)
  pval <- if (calculatePvalues) {
    tryCatch(
      fgsea:::gesecaCpp(E, queryVar, k,
                        sampleSize = 21, seed = 1, eps = 1e-300)[[1]]$pval,
      error = function(e) NA_real_)
  } else NA_real_
  data.table(gse = obj$gseId, gpl = obj$gplId,
             pctVar = queryVar / k / obj$totalVar * 100,
             pval = pval, size = k)
}

coresh_batch <- function(queries, chunk_dir, n_cores = 4, pvalues = FALSE) {
  stopifnot(is.list(queries), !is.null(names(queries)),
            all(nzchar(names(queries))))
  for (q in queries) stopifnot(is.integer(q), length(q) >= 3L)

  chunk_dir <- normalizePath(chunk_dir, mustWork = TRUE)
  chunk_paths <- list.files(chunk_dir, pattern = "full_objects\\.qs2$",
                            full.names = TRUE)
  if (length(chunk_paths) == 0L) stop(sprintf(
    "No *_full_objects.qs2 files in %s. See references/synapse-data-setup.md.",
    chunk_dir))

  # sanity-check first chunk
  first <- qs_read(chunk_paths[1])
  expected <- c("gseId", "gplId", "E1024", "rownames", "totalVar")
  missing <- setdiff(expected, names(first[[1]]))
  if (length(missing) > 0L) stop(sprintf(
    "Chunk %s missing fields: %s", chunk_paths[1], paste(missing, collapse = ", ")))
  stopifnot(is.integer(first[[1]]$rownames))

  bpparam <- MulticoreParam(n_cores, progressbar = TRUE)

  per_query <- lapply(names(queries), function(qname) {
    message(sprintf("[coresh_batch] query '%s' (k=%d, pvalues=%s)",
                    qname, length(queries[[qname]]), pvalues))
    ranking <- rbindlist(bplapply(chunk_paths, function(path) {
      ch <- qs_read(path)
      rbindlist(lapply(ch, coreshMatch,
                       query = queries[[qname]],
                       calculatePvalues = pvalues))
    }, BPPARAM = bpparam))
    ranking[, query_name := qname]
    if (pvalues) {
      ranking <- ranking[order(pval)]
    } else {
      ranking <- ranking[order(-pctVar)]
    }
    ranking[, rank := seq_len(.N)]
    ranking
  })

  result <- rbindlist(per_query)
  setcolorder(result, c("query_name", "gse", "gpl", "pctVar", "pval", "size", "rank"))
  result
}
