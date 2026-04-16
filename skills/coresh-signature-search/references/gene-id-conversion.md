# Gene ID conversion -- most common failure mode

CORESH chunks store `rownames` as **integer Entrez Gene IDs** (e.g. `7037L` for TFRC). If you pass symbols, `match("TFRC", c(7037L, 4891L, ...))` returns `NA` silently -- R does not warn on type mismatches like this. Your `size` column comes out zero, `pctVar` is noise, and nothing complains.

Always convert symbols -> Entrez integers BEFORE calling `coreshMatch()`, and assert on the result.

## R recipes

### Human

```r
library(org.Hs.eg.db)

sym2ent_human <- function(symbols) {
  mapped <- mapIds(org.Hs.eg.db,
                   keys = symbols,
                   keytype = "SYMBOL",
                   column = "ENTREZID",
                   multiVals = "first")
  # mapped is a named character vector; NA for unmapped symbols
  unmapped <- symbols[is.na(mapped)]
  if (length(unmapped) > 0) {
    warning(sprintf("%d/%d symbols failed to map: %s",
                    length(unmapped), length(symbols),
                    paste(unmapped, collapse = ", ")))
  }
  as.integer(na.omit(mapped))
}

query <- sym2ent_human(c("TFRC", "SLC11A2", "SLC11A1", "STEAP3", "STEAP4"))
# c(7037L, 4891L, 6556L, 55240L, 147339L)
```

### Mouse

```r
library(org.Mm.eg.db)

sym2ent_mouse <- function(symbols) {
  mapped <- mapIds(org.Mm.eg.db, keys = symbols,
                   keytype = "SYMBOL", column = "ENTREZID",
                   multiVals = "first")
  as.integer(na.omit(mapped))
}
```

Both are packaged together in `scripts/symbols_to_entrez.R` with a `species` argument.

## Cross-species orthologs

For cross-species validation (human signature -> query mouse compendium or vice versa), use an ortholog table before the Entrez conversion. Three reasonable options:

### Option A: `homologene` (Bioconductor, simple, static)

```r
library(homologene)

mouse_orthologs <- homologene(
  c("TFRC", "SLC11A2", "SLC11A1", "STEAP3", "STEAP4"),
  inTax = 9606, outTax = 10090    # human -> mouse
)
mouse_symbols <- mouse_orthologs$`10090`
mouse_query   <- sym2ent_mouse(mouse_symbols)
```

### Option B: `babelgene` (NCBI + Ensembl, more coverage)

```r
library(babelgene)

map <- orthologs(genes = c("TFRC", "SLC11A2"), species = "mouse")
mouse_query <- sym2ent_mouse(map$symbol)
```

### Option C: Ensembl BioMart (most thorough, slowest)

Use when `homologene`/`babelgene` coverage is insufficient for niche genes. See the `biomaRt` package vignette.

## Verifying the conversion on a real chunk

Always sanity-check against an actual chunk before the full sweep:

```r
chunk <- qs_read("preprocessed_chunks/hsa/001_full_objects.qs2")
stopifnot(is.integer(chunk[[1]]$rownames))     # integer, not character
stopifnot(length(query) ==
          length(na.omit(match(query, chunk[[1]]$rownames))))
```

If the second assertion fails, some query genes are not represented on that platform -- usually because the dataset is microarray and the platform doesn't probe the gene. Platform coverage varies; `size` on the final ranking tells you how consistently your genes appear across the compendium.

## Python (for signature preparation outside R)

```python
import mygene
mg = mygene.MyGeneInfo()

def sym2ent(symbols, species="human"):
    res = mg.querymany(symbols, scopes="symbol",
                       fields="entrezgene", species=species)
    return [int(r["entrezgene"]) for r in res if "entrezgene" in r]

query = sym2ent(["TFRC", "SLC11A2"], species="human")
# [7037, 4891]
```

`mygene` queries the MyGene.info API online, so it's slightly heavier than `org.Hs.eg.db` but up-to-date with NCBI.

## Reverse direction: Entrez -> symbols for output

When extracting gene loadings from top hits (the CoReSh-to-GSEA bridge), you end up with integer Entrez IDs and need symbols for downstream GSEA with MSigDB-style gene sets. Same `mapIds` call, reversed:

```r
ent2sym_human <- function(entrez_ids) {
  mapIds(org.Hs.eg.db,
         keys = as.character(entrez_ids),     # keys must be character!
         keytype = "ENTREZID",
         column = "SYMBOL",
         multiVals = "first")
}
```

Note the `as.character()` -- `mapIds` accepts character keys only, even for numeric ID types.
