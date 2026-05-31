# mLLMCelltype — R / Seurat Workflow

Install:

```r
install.packages("mLLMCelltype")                                   # CRAN
# or: devtools::install_github("cafferychen777/mLLMCelltype", subdir = "R")
```

The R `interactive_consensus_annotation()` takes `FindAllMarkers()` output directly
as its `input` (it expects the `cluster`/`gene` columns Seurat produces), so no
manual dict construction is needed.

```r
library(mLLMCelltype)
library(Seurat)
library(dplyr)

pbmc <- readRDS("your_seurat_object.rds")     # already clustered

# Marker genes per cluster (positive markers only)
pbmc_markers <- FindAllMarkers(pbmc, only.pos = TRUE,
                               min.pct = 0.25, logfc.threshold = 0.25)

# Cache speeds up identical re-runs
cache_dir <- "./mllmcelltype_cache"
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

consensus_results <- interactive_consensus_annotation(
  input       = pbmc_markers,
  tissue_name = "human PBMC",                  # NOTE: R uses tissue_name, Python uses tissue
  models      = c("claude-sonnet-4-6",         # Anthropic
                  "gpt-5.5",                    # OpenAI
                  "gemini-3.1-pro-preview",     # Google
                  "qwen3.6-plus"),              # Alibaba
  api_keys    = list(
    anthropic = "your-anthropic-key",
    openai    = "your-openai-key",
    gemini    = "your-gemini-key",
    qwen      = "your-qwen-key"
  ),
  consensus_threshold   = 0.7,
  max_discussion_rounds = 3,
  cache_dir             = cache_dir
)

# Map consensus labels back onto the Seurat object
final <- consensus_results$consensus                   # named vector: cluster -> label
pbmc$cell_type <- final[as.character(Idents(pbmc))]
```

## Key R/Python differences

| | Python | R |
|---|---|---|
| Marker input | `marker_genes` dict `{cluster: [genes]}` | `input` = `FindAllMarkers()` data frame |
| Tissue arg | `tissue=` | `tissue_name=` |
| API keys | env vars or `api_keys` dict | `api_keys = list(provider = key)` |
| Result access | `res["consensus"]` | `consensus_results$consensus` |

Gene symbols (not Ensembl IDs), cluster-ID round-tripping, and the uncertainty
caveats from the main SKILL.md apply identically in R.

R docs: https://cafferyang.com/mLLMCelltype/
