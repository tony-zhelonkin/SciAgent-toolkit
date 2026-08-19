# Embedding Deliverables

Use this guide when a pipeline produces both a de novo discovery embedding and a
supervised or batch-integrated embedding.

## Preserve both analytical roles

Package both embeddings under role-explicit names:

| AnnData key | Basis | Primary question |
|---|---|---|
| `X_umap_unsupervised` | de novo representation | What cell states and continua are present? |
| `X_umap_integrated` | supervised or corrected representation | How do cell types compare across samples or conditions? |

The exact upstream representation belongs in metadata. The key names describe
the role a collaborator can rely on.

Supervised and corrected spaces can sharpen shared cell-type structure while
compressing within-label variation. The unsupervised space retains a discovery
view of sub-states and continua. Keeping both prevents one analytical purpose
from silently standing in for the other.

## Packaging

Before writing the canonical object:

```python
required = ("X_umap_unsupervised", "X_umap_integrated")
for key in required:
    coords = adata.obsm[key]
    assert coords.shape[0] == adata.n_obs
    assert np.isfinite(coords).all()
```

Record in `uns["embeddings"]` for each embedding:

- source representation and method;
- cells included during fitting;
- whether query cells were projected or jointly fit;
- batch, label, or reference information used; and
- the scientific role the embedding supports.

## Cross-format names

Carry the roles into derived objects with stable, documented names:

- Seurat: `umap.unsup` and `umap.integrated`;
- Loupe: `umap_unsupervised` and `umap_integrated`; and
- CellxGene: schema-compatible names whose labels preserve the two roles.

If a target format imposes a default embedding, choose it explicitly for the
collaborator workflow and retain the other as an available alternative.

## Completeness checks

Validate all cells, including treated, query, rare, and reference subsets:

- coordinate matrices have one row per observation;
- row order matches `obs_names` after save and reload;
- coordinates are finite;
- no subset was silently excluded during projection;
- role labels and method metadata survive conversion; and
- representative biological structures are visible in the expected space.

## Interpretation handoff

The objects README should tell collaborators which embedding to open for cell
state discovery and which to use for integrated cell-type comparison. It should
also state that the spaces answer different questions and name any known
projection or reference limitations.

## Verification

Open each derived format independently. Compare cell order, coordinates for a
sample of cells, embedding names, and default-view selection against the
canonical AnnData object.
