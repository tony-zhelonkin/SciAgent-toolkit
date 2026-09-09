# Gene-Identifier Packaging

Use this guide when finalizing the active variable names in a collaborator-facing
AnnData object.

## Working identity and delivery identity

Keep stable Ensembl identifiers as `var_names` during build and analysis stages.
They support unambiguous joins, gene-union construction, annotation lookup, and
gene-set matching.

At the packaging milestone:

1. copy the active Ensembl identifiers into `var["gene_id"]`;
2. make the approved gene symbols the active `var_names`;
3. make duplicate symbols unique with a stable, documented suffix rule; and
4. record the round-trip semantics in `uns` and the objects README.

```python
adata.var["gene_id"] = adata.var_names.astype(str)
adata.var_names = adata.var["gene_name"].astype(str)
adata.var_names_make_unique(join="-")
adata.uns["var_names_note"] = (
    "var_names are unique gene symbols; Ensembl identifiers are in "
    "var['gene_id']"
)
```

The switch occurs once during finalization. Intermediate objects retain the
identifier system their upstream joins and models expect.

## Packaging contract

Before the switch, validate:

- every variable has a non-empty stable identifier;
- symbol provenance and genome build are recorded;
- missing symbols have an explicit policy;
- duplicate symbols are understood before uniquification; and
- any feature filtering is complete.

After the switch, validate:

- `var_names` are unique;
- `var["gene_id"]` is complete and preserves pre-switch order;
- the symbol-to-identifier mapping survives save and reload;
- named-gene queries used by collaborators resolve as expected; and
- downstream conversions preserve both fields.

## Derived formats

The Seurat export keeps symbols as feature names and carries `gene_id` in feature
metadata. Loupe and CellxGene exports should preserve both identifiers to the
extent their schemas permit. Document any format limitation rather than silently
dropping the stable identifier.

## Why packaging is the boundary

Scientists commonly query delivered objects by familiar gene symbol. Build and
analysis code benefit from stable identifiers. A single explicit packaging step
serves both audiences while avoiding mid-pipeline identifier changes that can
silently alter joins, matrices, or gene-set membership.

## Verification

Sample genes with unique, duplicated, and missing symbols. Trace each from the
last checkpoint through the canonical object and every derivative. Confirm row
order, matrix dimensions, active names, and preserved identifiers remain aligned.
