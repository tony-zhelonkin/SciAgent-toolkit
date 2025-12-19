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
| `pval` | **Raw p-value** | Using padj → wrong BUM scoring |
| `baseMean` | **Linear scale** | Forgetting `2^AveExpr` from limma |
| `ID` | Must match org.anno | Wrong ID type → empty graph (check console output) |

Handle duplicates: keep lowest p-value per gene.

## Network Files
| Network | Required Files |
|---------|----------------|
| KEGG | `network.kegg.rds`, `met.kegg.db.rds` |
| Combined | Above + `gene2reaction.combined.*.tsv` |
| Rhea-lipid | `network.rhea.lipids.rds`, `met.lipids.db.rds`, gene2reaction file; use `topology="metabolites"` |

Download: `http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/`

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