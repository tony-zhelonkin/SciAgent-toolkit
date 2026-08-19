# Metadata Packaging

Use this guide for the deliberate final cleanup of AnnData `obs` and for
preserving uncertainty in consensus annotations.

## One packaging pass

Keep working columns stable while stages still join, filter, and model against
them. Perform one reviewed cleanup when creating the canonical deliverable:

1. deduplicate identity columns into approved canonical names;
2. remove derivable working columns after confirming their sources remain;
3. move truly object-level, single-valued metadata into `uns`;
4. apply collaborator-approved labels and factor order; and
5. document every delivered field in the objects README.

Implement the pass as one named operation whose mapping is inspectable and
tested. Preserve source columns until the final call so upstream stages retain
their expected schemas.

## Canonical identity fields

Choose one field for each identity role used by the project, such as sample,
individual, pool, condition, batch, and final cell label. Before dropping aliases,
verify they agree on every observation or record the reconciliation rule.

Label changes are data transformations. Store the mapping in configuration or a
reviewed table, including factor order and the meaning of missing or unknown
values.

## Derivable and object-level fields

Remove a derived working field only when:

- its source values and transformation are preserved;
- no downstream consumer declares it as part of the contract; and
- regenerating it produces the same value and type.

Move a single-valued field to `uns` when it describes the object as a whole.
Retain it in `obs` when row-level consumers require it, even if the current object
happens to contain one value.

## Annotation consensus: preserve abstention

When combining label-transfer methods, preserve confidently unplaceable cells as
an explicit uncertain or novel state. A conservative policy protects treated or
query populations whose biology is absent from the training labels.

Read confidence thresholds and the consensus policy from project configuration.
The package should record:

- component methods and reference versions;
- confidence fields and their meanings;
- the consensus rule and abstention label;
- final label provenance; and
- uncertain or novel fractions by relevant condition.

Reference recall on known coarse types establishes agreement with the reference
categories. Also inspect whether condition-specific or rare states accumulate in
the abstention set, because forced assignment can erase the signal the experiment
was designed to find.

## Schema verification

After cleanup and reload, confirm:

- canonical identity columns are present, typed, and complete as intended;
- aliases slated for removal truly agree with retained fields;
- category levels and label mappings match the approved table;
- object-level fields appear in `uns` with clear names;
- uncertainty states and confidence values survive every conversion; and
- the objects README describes each retained `obs` field and label set.

Compare cell order and values across AnnData, Seurat, Loupe, and CellxGene
derivatives. Document target-format limitations and retain the canonical AnnData
object as the full-fidelity source.
