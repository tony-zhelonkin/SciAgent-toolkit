# Toolkit lanes

Lanes are a compact map from a user's situation to related skills. Skill descriptions provide the
actual routing triggers, and each skill body carries its own boundaries and related-skill links.

| Lane | What it is for | Composes with |
|---|---|---|
| `base` | Day-to-day bioinformatics: scverse and scRNA-seq foundations, bulk RNA-seq interpretation, scATAC, multiomics, GRNs, figures, pipeline conventions, and helpers for docs, biology, exploration, captions, curation, review, and handoff. | Any specialist analysis lane; `science-architect` for orchestrated work. |
| `architect` | Human-gated, architecture-first software changes: map, expert review, synthesis, design, phased planning, implementation handoff, verification, portfolio decisions, and retrospective audits whose artifacts compound downstream. | `base` when software work accompanies analysis; `pathway-signature` when interpretation accompanies architecture work. |
| `science-architect` | Multi-phase analysis plans: Opus decomposes and reviews, Sonnet implements one-script phases, and durable plan/research artifacts, figure variants, interpretation campaigns, and human decision gates carry the work. | `base` and the specialist lane that owns the scientific work. |
| `multiome` | Paired 10x Multiome RNA+ATAC from shared barcodes: primary processing, joint embeddings and WNN, consensus peaks, motif enrichment, peak-to-gene linkage, and coverage tracks. | `base`; `scatac-regulatory` for deeper chromatin analysis; `pathway-signature` for downstream functional interpretation. |
| `pathway-signature` | Functional interpretation of pseudobulk differential expression: multi-database GSEA, TF and pathway activity, public-data signature search, and interactive dashboards. | `base` and any assay lane that produces differential-expression results. |
| `rnaseq-fastq-preprocessing` | Bulk RNA-seq FASTQs through nf-core/rnaseq to BAMs and gene counts, with the Random-One STAR and grouped-SAF featureCounts path for TE-compatible integer counts. | `base` for downstream annotation and DE; `pathway-signature` for enrichment and interpretation. |
| `scatac-regulatory` | Chromatin-accessibility-first analysis: differential accessibility, sub-peak CREs, motif activity and enrichment, TF footprinting, peak-to-gene linkage, and ATAC-driven GRNs. | `base`; `multiome` when RNA and ATAC share barcodes; `pathway-signature` for downstream functional interpretation. |
