# mLLMCelltype — Uncertainty Metrics, Triage, and Hierarchical Annotation

## The two uncertainty metrics

For every cluster, `interactive_consensus_annotation` returns:

- **Consensus proportion** (`consensus_proportion`, 0–1): the fraction of models that
  agree on the final label. High = models converge; low = disagreement.
- **Shannon entropy** (`entropy`, ≥ 0): spread of the vote distribution across
  candidate labels. Low (≈0) = unanimous; high = votes scattered over many types.

A cluster is treated as **controversial** (and sent to discussion rounds) when its
proportion falls below `consensus_threshold` (default 0.7) OR its entropy exceeds
`entropy_threshold` (default 1.0).

## Triage pattern

```python
review = {
    c for c in res["consensus"]
    if res["consensus_proportion"][c] < 0.7 or res["entropy"][c] > 1.0
}
print("Manually validate:", review)

# Confirm with marker expression for the flagged clusters
import scanpy as sc
sc.pl.dotplot(adata, var_names=top_markers_for(review), groupby="leiden")
```

Interpretation guide:

| Proportion | Entropy | Action |
|---|---|---|
| high (≥0.8) | low (≈0) | Trust the label; spot-check only |
| moderate | moderate | Accept but verify markers |
| low (<0.7) | high (>1.0) | Likely mixed cluster, rare/novel type, or ambiguous markers — re-examine |

Low confidence often signals a **doublet/mixed cluster** or **insufficient marker
specificity** rather than a hard biological question — increasing `max_discussion_rounds`
rarely fixes a genuinely mixed cluster.

## Uncertainty-driven re-clustering

When a cluster is persistently high-entropy, sub-cluster it and re-annotate:

```python
sub = adata[adata.obs["leiden"] == "7"].copy()
sc.pp.neighbors(sub); sc.tl.leiden(sub, resolution=1.0, key_added="subleiden")
sc.tl.rank_genes_groups(sub, "subleiden", method="wilcoxon")
sub_markers = {c: sub.uns["rank_genes_groups"]["names"][c][:10].tolist()
               for c in sub.obs["subleiden"].cat.categories}
sub_res = interactive_consensus_annotation(sub_markers, species="human", tissue="blood",
                                           models=[...])
```

## Hierarchical / multi-resolution annotation

mLLMCelltype supports multi-resolution analysis with consistency checks: annotate at a
coarse resolution to get major lineages, then at a finer resolution for subtypes, and
pass coarse calls as `additional_context` to the fine pass so the models keep
parent/child labels consistent (e.g. fine "CD8+ effector T" must sit under coarse "T cell").

```python
coarse = interactive_consensus_annotation(coarse_markers, species="human",
                                           tissue="blood", models=[...])
fine = interactive_consensus_annotation(
    fine_markers, species="human", tissue="blood", models=[...],
    additional_context=f"Coarse-resolution labels for the parent clusters: {coarse['consensus']}. "
                       f"Keep fine subtypes consistent with these lineages.",
)
```

## Reporting

`format_discussion_report(res)` renders the full deliberation (per-model votes,
discussion turns, final rationale) as Markdown — attach it to methods/supplement for
transparency, since recording the reasoning process is a design goal of the tool.
