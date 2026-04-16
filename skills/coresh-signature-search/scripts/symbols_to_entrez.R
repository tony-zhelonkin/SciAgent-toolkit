# symbols_to_entrez.R -- gene symbol -> integer Entrez ID helper for CORESH
#
# Usage (inside R):
#   source("symbols_to_entrez.R")
#   query <- sym2ent(c("TFRC", "SLC11A2"), species = "human")
#
# Usage (CLI):
#   Rscript symbols_to_entrez.R human TFRC SLC11A2 STEAP3
#
# Warns on unmapped symbols. Returns integer vector (possibly empty).
# Covers human (org.Hs.eg.db) and mouse (org.Mm.eg.db).

sym2ent <- function(symbols, species = c("human", "mouse")) {
  species <- match.arg(species)
  db <- switch(species,
               human = { requireNamespace("org.Hs.eg.db", quietly = TRUE);
                         org.Hs.eg.db::org.Hs.eg.db },
               mouse = { requireNamespace("org.Mm.eg.db", quietly = TRUE);
                         org.Mm.eg.db::org.Mm.eg.db })
  if (is.null(db)) stop(sprintf(
    "org.%s.eg.db not installed. Install via BiocManager::install('org.%s.eg.db').",
    ifelse(species == "human", "Hs", "Mm"),
    ifelse(species == "human", "Hs", "Mm")))

  mapped <- AnnotationDbi::mapIds(
    db, keys = symbols, keytype = "SYMBOL",
    column = "ENTREZID", multiVals = "first")

  unmapped <- symbols[is.na(mapped)]
  if (length(unmapped) > 0) {
    warning(sprintf("%d/%d symbols failed to map: %s",
                    length(unmapped), length(symbols),
                    paste(utils::head(unmapped, 20), collapse = ", ")),
            call. = FALSE)
  }
  as.integer(stats::na.omit(mapped))
}

ent2sym <- function(entrez_ids, species = c("human", "mouse")) {
  species <- match.arg(species)
  db <- switch(species,
               human = org.Hs.eg.db::org.Hs.eg.db,
               mouse = org.Mm.eg.db::org.Mm.eg.db)
  AnnotationDbi::mapIds(
    db, keys = as.character(entrez_ids),
    keytype = "ENTREZID", column = "SYMBOL", multiVals = "first")
}

# CLI entrypoint
if (!interactive() && !is.null(sys.call(0))) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 2) {
    species <- args[1]
    symbols <- args[-1]
    result <- sym2ent(symbols, species = species)
    cat(paste(result, collapse = ","), "\n")
  }
}
