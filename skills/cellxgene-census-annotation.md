# CELLxGENE Census - Cell Type Annotation via Reference Data skill

## Purpose
Census provides cloud-hosted access to >33M annotated single cells with pre-computed embeddings (scVI, Geneformer), enabling annotation transfer from reference to query data without downloading massive datasets. Primary use: leverage existing annotations from similar cells to annotate new data.

## Core Workflow: Annotation Transfer

```python
import cellxgene_census

# CRITICAL: Always open with version for reproducibility
census = cellxgene_census.open_soma(census_version="2023-12-15")

# Query with value_filter (TileDB syntax, NOT pandas)
adata_ref = cellxgene_census.get_anndata(
    census=census,
    organism="Homo sapiens",
    obs_value_filter="tissue == 'lung' and is_primary_data == True",  # Remove duplicates
    obs_embeddings=["scvi"],  # Pre-computed embeddings
    column_names={"obs": ["cell_type", "tissue"]}  # Minimize metadata
)

census.close()  # Always close
```

## Critical Non-Obvious Patterns

### 1. Duplicate Cell Filter (ESSENTIAL)
```python
# Census contains duplicate cells across datasets
# ALWAYS include: is_primary_data == True
obs_value_filter="cell_type == 'neuron' and is_primary_data == True"
```

### 2. Query Syntax Gotchas
```python
# Use Python operators, NOT pandas syntax
# ✓ "sex == 'female' and cell_type in ['B cell', 'neuron']"
# ✗ "sex.isin(['female'])" - TileDB QueryCondition, not pandas
```

### 3. scVI Model Integration (Non-Standard)
```python
import scvi

# Prepare query data for Census model
scvi.model.SCVI.prepare_query_anndata(query_adata, "model_path")
vae_q = scvi.model.SCVI.load_query_data(query_adata, "model_path")

# KEY TRICK: Bypass training requirement
vae_q.is_trained = True  
query_adata.obsm["scvi"] = vae_q.get_latent_representation()
```

### 4. Vector Search for Annotation Transfer
```python
# EXPERIMENTAL API - requires cellxgene_census[experimental]
neighbors = cellxgene_census.experimental.find_nearest_obs(
    embedding_name="scvi",
    organism="homo_sapiens", 
    census_version="2023-12-15",
    query=query_adata,  # Must have matching obsm layer
    k=30,  # More neighbors = better confidence
    nprobe=20  # Search accuracy vs speed tradeoff
)

# Returns NeighborObs with distances and neighbor soma_joinids
# Use predict_obs_metadata for majority-vote annotation
predictions = cellxgene_census.experimental.predict_obs_metadata(
    organism="homo_sapiens",
    census_version="2023-12-15", 
    neighbors=neighbors,
    attributes=["cell_type", "tissue_general"]  # Returns confidence scores
)
```

### 5. Geneformer Tokenization Pattern
```python
from geneformer import TranscriptomeTokenizer

# REQUIRED format: ENSEMBL IDs WITHOUT version (e.g., ENSG00000139618 not ENSG00000139618.2)
adata.var["ensembl_id"] = adata.var.index.str.split(".").str[0]

# Custom attributes for join tracking
tokenizer = TranscriptomeTokenizer(custom_attr_name_dict={"joinid": "joinid"})
tokenizer.tokenize_data(
    data_directory="h5ad_dir/",
    output_directory="token_dir/",
    file_format="h5ad"
)
```

## Domain Context

**Census Hierarchy:**
```
census["census_data"]["homo_sapiens"]
    .obs                    # Cell metadata (SOMADataFrame)
    .ms["RNA"].var          # Gene metadata  
    .ms["RNA"].X("raw")     # Raw counts matrix
    .obsm["scvi"]          # scVI embeddings (if requested)
```

**Embedding Models:**
- `scvi`: Continuous latent space, good for fine-grained cell states
- `geneformer`: Token-based, requires special preprocessing, better for cell classes

**Confidence Interpretation:**
- Prediction confidence <0.5: Unreliable, cell may be rare/novel type
- High k (30-50) neighbors: More robust but slower, better for ambiguous cells
- `is_primary_data==True` removes technical replicates but may exclude biologica replicates

## Common Pitfalls

**Memory Management:**
```python
# Census queries load into memory - calculate size first
# Rough estimate: 500MB per 100K cells with full gene set
# Use column_names to limit metadata: {"obs": ["cell_type"]}
# Use var_value_filter to limit genes: "feature_id in ['ENSG00000161798']"
```

**Model Version Mismatch:**
```python
# scVI model genes ≠ query genes → prepare_query_anndata handles this
# Check coverage: typically need >50% overlap
# Model reports: "Found X% reference vars in query data"
```

**soma_joinid vs Index:**
```python
# soma_joinid is Census internal ID for joins
# NOT sequential, NOT preserved across versions
# Use obs_coords parameter for specific cell retrieval:
obs_coords=sorted(neighbors.neighbor_ids[:, 0].tolist())
```

**Census Version Stability:**
```python
# Weekly releases are transient, use long-term support versions
# LTS versions: "2023-07-25", "2023-12-15", "2024-07-01"
# Always specify census_version for reproducibility
```

## Integration with Standard Workflows

**After Census Annotation Transfer:**
1. Validate predictions with marker gene expression
2. Refine ambiguous calls with sub-clustering
3. Cross-reference confidence scores with QC metrics
4. Use Census predictions as prior for probabilistic models

**Typical Pipeline:**
```
Query Data → Embed with Census Model → Vector Search → 
→ Predict Labels → Validate with Markers → Integrate to AnnData
```

**Alternative: Direct Random Forest (No Vector Search):**
```python
from sklearn.ensemble import RandomForestClassifier

# Train on Census embeddings
rfc = RandomForestClassifier()
rfc.fit(census_adata.obsm["scvi"], census_adata.obs["cell_type"])

# Predict on query
query_adata.obs["predicted_cell_type"] = rfc.predict(query_adata.obsm["scvi"])

# Get confidence (probability of predicted class)
probs = rfc.predict_proba(query_adata.obsm["scvi"])
conf = [probs[i][rfc.classes_ == pred[i]][0] 
        for i, pred in enumerate(query_adata.obs["predicted_cell_type"])]
```