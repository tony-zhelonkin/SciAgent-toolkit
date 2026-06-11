---
name: annotate-bulk-rnaseq-data
description: 'Router for annotating bulk RNA-seq featureCounts matrices before edgeR/limma DE. Branches into two paths: regular gene-symbol annotation (Ensembl/biomaRt, the usual case) and transposable-element annotation (parse TE IDs into Subfamily:Family:Class, build TE DGEList), with a combined gene+TE matrix when both are needed. R-based; enforces the rule "annotate before filtering." Use when preparing featureCounts gene (and TE) outputs for differential expression. Routes to references/gene-annotation.md (gene path) and references/te-annotation.md (TE path). For upstream STAR/featureCounts TE preprocessing use star-te-preprocessing. For single-cell count matrices use single-cell-rna-qc instead.'
license: MIT
metadata:
  scope: implementation
  requires: []
  skill-author: SciAgent-toolkit
  last-reviewed: 2026-06-11
  # CHANGELOG: 1.1.1 — added joint gene+TE analysis caveats at the combined-matrix handoff
  #            (genes-only size factors, within-feature-type-only validity, no gene-vs-TE
  #            magnitude comparison, no TPM for TEs).
  #            1.1.0 — split into thin router + references/{gene,te}-annotation.md;
  #            fixed TE-ID label to Subfamily:Family:Class; added star-te-preprocessing back-edge.
  version: 1.1.1
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

## Joint gene+TE matrix — analysis caveats (read before combined-mode DE)

A row-bound gene+TE matrix is valid only because exonic TE loci were subtracted upstream (no double-counting) — but mutual exclusivity is **necessary, not sufficient**. When you take a combined matrix into joint normalization/DE, the rules below are load-bearing — but they are **graded options, not mandates** (grade scale + key claims + gaps live in `te-gene-featurecounts/SKILL.md` "Evidence & open questions", authoritative source: the evidence-graded reconciliation, note 13). The *joint matrix itself* is a reviewed construct (**grade A**).

- **Size factors from genes ONLY — grade B / CONTESTED.** Estimate DESeq2 size factors on the gene submatrix (`estimateSizeFactors(dds, controlGenes = which(feature_type == "gene"))`). TE-Seq advocates this (the long-tailed, multimapper-inflated (`-M`) TE minority can violate the "most features unchanged" assumption and drag *gene* fold-changes); but the dominant tool **TEtranscripts pools** genes+TEs. Reasonable, not universal — sanity-check against pooled size factors and confirm gene LFCs are stable.
- **Within-feature-type, across-sample DE ONLY — grade C / inference.** The combined object is valid for gene-vs-sample and TE-vs-sample comparison; the per-sample basis offset and multimapper bias cancel across samples within a feature type. Sound mechanistic inference, not stated in a TE primary source.
- **NEVER compare gene-vs-TE magnitude within a sample — grade C / inference.** Genes (unique-only, possibly Salmon/EM/length-modeled) and TEs (`-M` integer, no length model, possibly `-s 0` both-strand) sit on different measurement bases. "This TE is expressed like that gene" and "TE % of transcriptome" are not interpretable as biology.
- **No TPM/FPKM for TE meta-features — grade C / inference.** A summed multi-locus subfamily has no single length, so length-normalized units are undefined. Use model-normalized counts / logCPM / DESeq2 LFCs only.

If TE rows were counted stranded with a sense/antisense split upstream (`--te-strand sense_antisense`), keep sense and antisense as **separate** TE features — the more principled best-practice (**grade B / SQuIRE-specific**) for preserving bidirectional TE biology on a stranded library without breaking gene-comparability. The field is split (TEtranscripts defaults `--stranded no`); `-s 0` standalone matrices remain a valid choice, and mode-switching is not required.

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
