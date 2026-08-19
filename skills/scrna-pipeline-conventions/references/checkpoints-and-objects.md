# Checkpoints and Objects

Use this guide for AnnData recovery boundaries and finalized scRNA deliverables.
Read `analysis-code-conventions/references/restartability-and-dataflow.md` for
the general dependency, invalidation, and safe-write doctrine.

## Two object lifecycles

Keep intermediate recovery state distinct from the reviewed milestone:

```text
03_results/
├── checkpoints/    # intermediate .h5ad state consumed by later stages
└── objects/        # finalized milestone object and its documentation
```

A checkpoint is keyed to the analysis state it represents and may be regenerated.
The canonical object is named for the project or deliverable and is the source
from which collaborator formats are derived.

## Checkpoint boundaries

Write a checkpoint after an expensive or failure-prone scRNA operation when a
later run can validate and resume from that state. Common boundaries include:

- ingestion and gene-union construction;
- cell and gene QC;
- ambient-RNA or doublet handling;
- normalized or latent representations needed by later stages; and
- annotation state before final packaging.

Each stage reads one declared upstream object and writes a new object. Record the
input identity, consequential configuration, tool versions, shape, and validation
summary with the checkpoint or in its sidecar metadata.

Do not overwrite the input object in place. Write a temporary object, reopen it,
validate it, and then promote it to the intended checkpoint path.

## Checkpoint validation

Before reuse, verify the invariants the next stage assumes, including as relevant:

- unique observation and variable names;
- expected cell and gene counts;
- required `obs`, `var`, `layers`, `obsm`, and `uns` keys;
- sparse/dense representation expectations;
- absence or documented meaning of missing coordinates and labels; and
- agreement between recorded dependencies and current project state.

An existing `.h5ad` is evidence that a write occurred. Validation establishes
whether it is a usable recovery point.

## Finalized milestone object

Package one canonical annotated AnnData object under `03_results/objects/` and
document it in the sibling `README.md`. The package review covers:

- active and preserved gene identifiers;
- both required embedding roles when the pipeline produces them;
- cleaned and documented observation metadata;
- layers and raw-count semantics;
- cell order and identifier uniqueness; and
- the provenance needed to regenerate it from checkpoints.

The final filename should identify the project or dataset rather than an
intermediate stage. Treat this object as the source for every derived format.

## Derivation chain

Generate Seurat, Loupe, and CellxGene forms from the canonical AnnData object:

```text
canonical .h5ad → Seurat .rds → Loupe .cloupe
        └───────→ CellxGene-ready .h5ad
```

Preserve cell order, gene identifiers, embeddings, and documented metadata across
each conversion. Make corrections in the canonical object and regenerate its
derivatives so the chain has one authoritative source.

## Verification

1. Resume from every supported checkpoint and confirm its dependency record.
2. Reopen the finalized `.h5ad` in a fresh process and run the package invariants.
3. Compare observation order and required fields across every derivative.
4. Confirm the objects README names the source object, derivation commands, schema,
   embedding roles, and known limitations.
