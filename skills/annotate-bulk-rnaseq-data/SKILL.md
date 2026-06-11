---
name: annotate-bulk-rnaseq-data
description: 'Router for annotating bulk RNA-seq featureCounts matrices before edgeR/limma DE. Branches into two paths: regular gene-symbol annotation (Ensembl/biomaRt, the usual case) and transposable-element annotation (parse TE IDs into Subfamily:Family:Class, build TE DGEList), with a combined gene+TE matrix when both are needed. R-based; enforces the rule "annotate before filtering." Use when preparing featureCounts gene (and TE) outputs for differential expression. Routes to references/gene-annotation.md (gene path) and references/te-annotation.md (TE path). For upstream STAR/featureCounts TE preprocessing use star-te-preprocessing. For single-cell count matrices use single-cell-rna-qc instead.'
license: MIT
metadata:
  scope: implementation
  requires: []
  skill-author: SciAgent-toolkit
  last-reviewed: 2026-06-08
  # CHANGELOG: 1.1.0 — split into thin router + references/{gene,te}-annotation.md;
  #            fixed TE-ID label to Subfamily:Family:Class; added star-te-preprocessing back-edge.
  version: 1.1.0
  category: workflow
  tier: standard
  tags: []
  complementary-skills:
  - star-te-preprocessing
  - single-cell-rna-qc
  contraindications:
  - Do not filter low-count rows before annotation — you lose Ensembl IDs irreversibly.
  - Do not use for single-cell scRNA-seq. Use single-cell-rna-qc and scanpy/Seurat workflows instead.
---

# Annotate RNAseq Data (Genes + TEs)

## Overview

This skill annotates already-counted bulk RNA-seq matrices (featureCounts outputs) and assembles edgeR `DGEList` objects for downstream limma/edgeR DE. It is one logical step **between counting and DE**, and it splits into two paths:

- **Gene-symbol annotation (the usual case)** — every preprocessed dataset needs this: map Ensembl IDs to gene symbols via Ensembl/biomaRt and build the gene `DGEList`.
- **Transposable-element annotation (the rare case)** — parse TE IDs into `Subfamily:Family:Class`, build the TE `DGEList`, and optionally row-bind genes + TEs into a unified annotated matrix.

A gene-only agent reads this router plus `references/gene-annotation.md` and stops there; it never has to read TE-specific material.

**When to use this skill:**
- Annotate Ensembl IDs with gene symbols before filtering (gene path)
- Process transposable-element count matrices into family/subfamily annotations (TE path)
- Build gene and/or TE `DGEList` objects for downstream DE
- Combine gene + TE counts into a unified annotated matrix

**When NOT to use this skill:**
- Single-cell scRNA-seq count matrices → use `single-cell-rna-qc`
- Generating the count matrices themselves (STAR alignment, SAF building, featureCounts) → use `star-te-preprocessing`

---

## Decision Tree / Routing

```
Annotating a bulk RNA-seq count matrix?
│
├─ Annotating regular gene symbols (the usual case)?
│     →  references/gene-annotation.md
│        Ensembl→Symbol via biomaRt/org.db, gene DGEList assembly.
│
├─ Annotating transposable elements?
│     →  references/te-annotation.md
│        Parse TE IDs (Subfamily:Family:Class), TE annotation, TE DGEList.
│
└─ Both, for a combined gene+TE matrix?
      →  read BOTH references/gene-annotation.md and references/te-annotation.md
         (gene + TE blocks are row-bound into one annotated matrix).
```

---

## Shared Principle

**Always annotate BEFORE filtering.** Never drop low-count rows first — filtering before annotation loses Ensembl/TE IDs irreversibly, and the mapping back to symbols/families cannot be recovered. This rule holds on both paths.

---

## Complementary Skills

| When you need... | Use skill | Relationship |
|---|---|---|
| Produce the integer TE + gene count matrices this skill annotates (STAR Random-One, SAF, featureCounts) | `star-te-preprocessing` | Prerequisite (upstream step) |
| Single-cell scRNA-seq QC and annotation | `single-cell-rna-qc` | Alternative (different modality) |
| Downstream DE → GSEA on the annotated DGEList | `bulk-rnaseq-gsea` | Next step |

The canonical handoff chain is `star-te-preprocessing → annotate-bulk-rnaseq-data → bulk-rnaseq-gsea`.

---

## Resources

- **Gene path detail:** `references/gene-annotation.md`
- **TE path detail:** `references/te-annotation.md`
- **Gene helpers:** RNAseq-toolkit **v2.0.0** — `scripts/General/{io_helpers,annotate_genes,dge_helpers}.R`
- **TE helpers:** TE-RNAseq-toolkit **v2.0.0** — `scripts/te_utils.R`, TE-ID parser `R/te_utils.R::parse_te_id`
