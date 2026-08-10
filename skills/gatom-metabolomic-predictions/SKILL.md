---
name: gatom-metabolomic-predictions
description: "Finds maximally-regulated metabolic subnetworks from differential expression using atom-transition graphs and BUM-scored SGMWCS. Use when interpreting transcriptomic or metabolomic DE through KEGG/Rhea pathway structure, or linking enzyme expression to metabolite flow. R-based; needs raw p-values. For pathway enrichment use fgsea instead."
license: MIT
---

# GATOM - Active Metabolic Module Discovery

## Purpose
Finds maximally-regulated metabolic subnetworks from DE data using atom-transition networks. Solves SGMWCS via BUM model scoring—**not permutation testing**.

## Critical: Inverted Graph Structure
**Nodes = Metabolites** (atoms): Have `label`, `score`. **NO log2FC.**
**Edges = Reactions** (genes): Have `log2FC`, `pval`, `Symbol`, `enzyme`.

Extract genes from edges: `E(module)$Symbol`, not vertices.

## Core Workflow
```r
# 1. Prepare DE - EXACT FORMAT
gene.de <- data.frame(
 ID = symbols,        # Symbol, RefSeq, or Entrez
 pval = P.Value,      # RAW p-value (NOT padj!)
 log2FC = logFC,
 baseMean = 2^AveExpr # LINEAR scale
)

# 2. Build graph - met.db REQUIRED even without metabolite data
g <- makeMetabolicGraph(
 network = network, topology = "atoms",
 org.gatom.anno = org.anno, gene.de = gene.de,
 met.db = met.db,  # REQUIRED
 met.de = NULL,
 gene2reaction.extra = gene2reaction.extra  # Required for Rhea/Combined
)

# 3. Score and solve
gs <- scoreGraph(g, k.gene = 50, k.met = NULL)
m <- solve_mwcsp(rnc_solver(), gs)$graph

# 4. Extract genes from EDGES
module_genes <- igraph::as_data_frame(m, "edges")$Symbol
```

## Critical Parameters
| Parameter | Impact |
|-----------|--------|
| `k.gene` | Module size (heuristic). 25=stringent, 50=default, 75=exploratory. **Run sensitivity analysis.** |
| `gene2reaction.extra` | **REQUIRED for Rhea/Combined**—silent failure without it |
| `met.db` | **REQUIRED even with met.de=NULL**—topology needs it |

## Input Format Gotchas
| Field | Requirement | Common Error |
|-------|-------------|--------------|
| `pval` | **Raw p-value** (both `gene.de` and `met.de`) | Using padj → wrong BUM scoring |
| `baseMean` | **Linear scale, `gene.de` ONLY** — `met.de` has no `baseMean` requirement (verified against `ctlab/gatom` source, `R/de.R::getMetDEMeta`: its columns are `ID`/`pval`/`logPval`/`log2FC` plus an optional `signalColumn`, no `baseMean`) | Assuming `met.de` needs `2^AveExpr` too and blocking on a dataset that lacks it — it doesn't need it |
| `ID` | Must match org.anno (`gene.de`) or `met.db`'s base ID / one of its `mapFrom` ID types (`met.de`) | Wrong ID type → empty graph (check console output: `"Found DE table for ... with %s IDs"`) |

Handle duplicates: keep lowest p-value per gene.

## `met.de` ID space — verified per network, and where it breaks in practice

`met.de$ID` is auto-matched against `met.db`'s base ID or any of `met.db$mapFrom`'s alternate ID
types (`getMetDEMeta` / `findIdColumn`, `R/de.R`). Confirmed by inspecting `ctlab/gatom`'s own
tutorial vignette and its linked example DE files (not just this SKILL's prior wording):

| Network | `met.db` file | Base metabolite ID | `met.de$ID` confirmed to accept |
|---|---|---|---|
| KEGG | `met.kegg.db.rds` | KEGG compound | KEGG compound; also **HMDB** and ChEBI via `mapFrom` — confirmed from a real published example (`Ctrl.vs.MandLPSandIFNg.met.de.tsv.gz`) using HMDB IDs directly as `ID` |
| Rhea | `met.rhea.db.rds` | ChEBI | ChEBI |
| Combined | `met.combined.db` | KEGG + ChEBI | KEGG or ChEBI |
| Rhea-lipid subnetwork (`topology="metabolites"`) | `met.rhea.lipids.db.rds` | — | **Full lipid/fatty-acid common names** (e.g. `15S-HETE`, `10Z-Heptadecenoic acid`) — confirmed from GATOM's own worked lipidomics example (`Ctrl.vs.HighFat.lipid.de.csv`, real ST001289 mouse data). **NOT** LIPID MAPS shorthand class notation (`DG 18:0_20:4`, `PC 38:4`) — the default output style of LipidMatch/LipidAnnotator/MS-DIAL. Shorthand IDs will silently fail to match (empty/near-empty graph) unless translated to common names or ChEBI first.

**Real-world ID-format inconsistency warning — do not hardcode an assumed format.** ID systems
drift version-to-version and a given `met.db` snapshot may not match what your DE table uses even
when the *system* (e.g. "HMDB") is nominally the same:

- HMDB accessions changed format over time: older 5-digit (`HMDB00634`, seen in the published KEGG
  example above) vs. current 7-digit, zero-padded (`HMDB0000008`). A `met.db` built against one
  format will not silently match IDs in the other — `findIdColumn`'s fuzzy sampling can fail
  quietly and route you straight into the "Empty graph" pitfall below with no error, only a console
  message stating which ID type it *did* detect.
- **Always inspect the actual loaded `met.db$mapFrom[[<idType>]]` (or `met.db$metabolites`) IDs
  after downloading/loading it, and compare a handful against your own `met.de$ID` values by eye,
  before assuming a network/ID-system pairing from this table (or from any past project's
  experience) will resolve.** This table records what was *verified to work in one example
  dataset*, not a guarantee for every `met.db` build or every metabolomics/lipidomics core
  facility's ID export.

## Network Files
| Network | Required Files |
|---------|----------------|
| KEGG | `network.kegg.rds`, `met.kegg.db.rds` |
| Combined | Above + `gene2reaction.combined.*.tsv` |
| Rhea-lipid | `network.rhea.lipids.rds`, `met.rhea.lipids.db.rds` (verified filename per `ctlab/gatom` tutorial vignette — some references shorten this to `met.lipids.db.rds`; check what actually downloads), gene2reaction file; use `topology="metabolites"` |

Download (manual): `http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/`

**Recommended download method (RNAseq-toolkit helper):**
```r
source(file.path(toolkit_dir, "scripts/GSEA/GSEA_processing/load_reference_db.R"))
download_gatom_references(dest_dir = "00_data/references/gatom")
# Downloads: network.kegg.rds, network.combined.rds, met.*.rds, org.Mm.eg.gatom.anno.rds, gene2reaction TSVs
```

## Sensitivity Analysis
k.gene is heuristic—run multiple values:
```r
core_genes <- Reduce(intersect, lapply(c(25, 50, 75), function(k) {
 gs <- scoreGraph(g, k.gene = k)
 unique(igraph::as_data_frame(solve_mwcsp(rnc_solver(), gs)$graph, "edges")$Symbol)
}))
```

## Common Pitfalls
| Problem | Fix |
|---------|-----|
| Empty graph | Check ID type—console shows "Found DE table for genes with X IDs" |
| `met.de` matches far fewer metabolites than expected (not fully empty, easy to miss) | Usually an ID-format-version mismatch, not a real absence of signal — see "`met.de` ID space" table above (old vs current HMDB format; lipid shorthand vs common name) |
| No genes in module | Extract from `E(m)$Symbol`, not vertices |
| Rhea network empty | Download organism-specific gene2reaction.*.tsv |
| Wrong viz colors | Color **edges** by log2FC, not nodes |

## Export
```r
write_graph(m, "module.graphml", format = "graphml")  # Cytoscape
saveModuleToHtml(m, "module.html", name = "Module")   # Interactive
```

---
*Install*: `BiocManager::install("gatom")`
*Test*: `data("networkEx")` then verify graph has edges

---

## When not to use

- Do not use adjusted p-values (padj) — GATOM's BUM scoring requires raw p-values.
- Do not use for pathway-level enrichment statistics. Use fgsea or clusterProfiler instead.

---

## See also

- `genenmf-metaprogram-discovery`

Upstream docs: https://artyomovlab.wustl.edu/publications/supp_materials/GATOM/
