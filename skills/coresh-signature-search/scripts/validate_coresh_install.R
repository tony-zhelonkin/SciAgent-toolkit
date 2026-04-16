# validate_coresh_install.R -- pre-flight check for CORESH setup
#
# Usage:
#   Rscript validate_coresh_install.R [chunk_dir]
#
# Default chunk_dir: $CORESH_CHUNKS env var, else "preprocessed_chunks".
#
# Checks (stops on first failure):
#   1. coresh package loads
#   2. BiocParallel, qs2, fgsea, data.table load
#   3. chunk_dir exists
#   4. hsa/ or mmu/ subdir has >= 1 *_full_objects.qs2 file
#   5. first chunk reads and contains expected per-object attributes
#   6. rownames on first object are integer (Entrez), not character

args <- commandArgs(trailingOnly = TRUE)
chunk_dir <- if (length(args) >= 1) args[1] else {
  Sys.getenv("CORESH_CHUNKS", unset = "preprocessed_chunks")
}

ok <- function(msg) cat("[ OK  ]", msg, "\n")
fail <- function(msg) { cat("[FAIL ]", msg, "\n"); quit(status = 1) }

# 1-2: packages
for (pkg in c("coresh", "BiocParallel", "qs2", "fgsea", "data.table")) {
  if (!requireNamespace(pkg, quietly = TRUE)) fail(
    sprintf("%s not installed. Install via remotes::install_github('alserglab/coresh') and BiocManager.", pkg))
}
ok("R packages loaded: coresh, BiocParallel, qs2, fgsea, data.table")

# 3: chunk_dir exists
if (!dir.exists(chunk_dir)) fail(sprintf(
  "chunk_dir '%s' not found. Run the Synapse download (references/synapse-data-setup.md).",
  chunk_dir))
ok(sprintf("chunk_dir exists: %s", normalizePath(chunk_dir)))

# 4: at least one species subdir with chunks
for (sp in c("hsa", "mmu")) {
  sp_dir <- file.path(chunk_dir, sp)
  if (!dir.exists(sp_dir)) { cat("[WARN ]", sp, "subdir missing\n"); next }
  chunk_paths <- list.files(sp_dir, pattern = "full_objects\\.qs2$", full.names = TRUE)
  if (length(chunk_paths) == 0L) {
    cat("[WARN ]", sp, "has 0 chunks\n")
  } else {
    ok(sprintf("%s/: %d chunks", sp, length(chunk_paths)))
  }
}

# 5-6: load first available chunk and verify structure
candidates <- unlist(lapply(c("hsa", "mmu"), function(sp) {
  list.files(file.path(chunk_dir, sp), pattern = "full_objects\\.qs2$", full.names = TRUE)
}))
if (length(candidates) == 0L) fail("No chunks in any subdir.")

suppressPackageStartupMessages(library(qs2))
first <- qs2::qs_read(candidates[1])
if (!is.list(first) || length(first) == 0L) fail(sprintf(
  "Chunk %s did not deserialize to a non-empty list.", candidates[1]))
obj <- first[[1]]
expected <- c("gseId", "gplId", "E1024", "rownames", "totalVar")
missing <- setdiff(expected, names(obj))
if (length(missing) > 0L) fail(sprintf(
  "First object in %s missing fields: %s", candidates[1], paste(missing, collapse = ", ")))
if (!is.integer(obj$rownames)) fail(
  "obj$rownames is not integer -- queries passed as symbols will silently return zero hits.")

ok(sprintf("First chunk verified: %d objects, first GSE=%s, rownames integer (Entrez).",
           length(first), obj$gseId))
cat("\n[ READY ] CORESH install looks good. Ready to run coresh_batch.R.\n")
