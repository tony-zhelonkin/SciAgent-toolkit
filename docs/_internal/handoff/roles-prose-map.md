# Role prose redistribution map

This map covers the 493 lines in the seven `roles/*.yaml` files. Step 4a leaves those files in
place. It relocates durable reader knowledge so their later deletion removes only the role,
activation, and provenance mechanism.

## Shared disposition

- Durable lane purpose and lane-to-lane composition live in `docs/lanes.md`.
- User-situation triggers live in the frontmatter descriptions of the skills named below. Existing
  descriptions already carried most of these triggers; the map identifies those existing homes as
  well as the descriptions changed in this step.
- Boundaries live in each relevant skill's single `## When not to use` section. Existing sections
  were extended where the role header supplied a missing boundary.
- Related skill links remain in `## See also`. New links in `delegate-cli` are plain skill names.
- The role `description:` fields condense their Purpose headers and now live in the matching
  `docs/lanes.md` row.
- The `skills:`, `agents:`, and `commands:` lists are provenance membership for the retiring role
  mechanism. Their leaf files remain in `skills/`, `agents/`, and `commands/`; the useful inline
  skill summaries are already represented by skill descriptions. The membership associations are
  deliberately dropped.

## Description capacity

`delegate-cli` was 344/350 characters and could not absorb the science-orchestration trigger as
written. Its description was compacted to 325 characters while retaining its direct codex/agy
requests, tool scope, and disambiguation, then adding the four orchestration triggers and durable
artifact/review terms. `architecture-first-dev` expanded from 318 to 349 characters. Other
applicable near-cap descriptions already stated their lane trigger, so they were left intact;
`decision-gate-notebook`, for example, is 346 characters and already names pipeline inflection
points, read-only re-plotting, and the durable approval gate.

## `architect`

| Header section | Destination |
|---|---|
| Purpose | `docs/lanes.md` (`architect`) records the design-before-code discipline and compounding map/review/design/plan artifacts. `skills/architecture-first-dev/SKILL.md` remains the routing spine. |
| Use when | All three triggers now appear in the rewritten `architecture-first-dev` description: non-trivial features/refactors/redesigns, cross-cutting expert review, and work where scientific reasoning matters alongside structure. |
| Do NOT use when | New `architecture-first-dev` `## When not to use` bullets cover one-line fixes, throwaway scripts, exploratory data analysis, and running bioinformatics pipelines. |
| Pipeline | `architecture-first-dev` contains the routing tree and stage contract; `docs/workflows/architect/00-quickstart.md` contains the readable end-to-end sequence and composable reviewer syntax. |
| Composes with | The `architect` row in `docs/lanes.md` records composition with `base` and `pathway-signature`. |

The YAML description maps to the same lane row and router description. Its agent and command
provenance is preserved functionally in the router references and architect workflow docs; role
membership itself is dropped.

## `base`

| Header section | Destination |
|---|---|
| Purpose | The `base` row in `docs/lanes.md` records the day-to-day scverse, bulk RNA-seq, scATAC, multiomics, GRN, figure, convention, and helper scope. |
| Use when | General single-cell and exploratory work is already routed by the `scanpy` description; bulk interpretation by `bulk-rnaseq-gsea` and `bulk-rnaseq-activity-inference`; new-project workflow setup by `scrna-pipeline-conventions`. |
| Do NOT use when | `scrna-pipeline-conventions` already redirects software-design work to `architecture-first-dev` and excludes throwaway exploration. `scanpy` now redirects pure pathway/TF interpretation to `bulk-rnaseq-gsea` or `bulk-rnaseq-activity-inference`. |

The YAML description maps to the `base` lane row. The seven helper-agent categories in Purpose and
their provenance list remain represented by the leaf files under `agents/`; the association with a
default role is dropped.

## `science-architect`

| Header section | Destination |
|---|---|
| Purpose | The `science-architect` row in `docs/lanes.md` records multi-phase, model-tiered execution, review gates, figure variants, interpretation campaigns, and durable artifacts. |
| Use when | The compacted `delegate-cli` description now names multi-phase worker fan-out and the `pipeline-plan`, `explore-and-plan`, `add-figure-variant`, and `interpret-storm` triggers. `decision-gate-notebook` and `interactive-breakpoint-explorer` descriptions already route pipeline inflection points and persistent review evidence. |
| Do NOT use when | `decision-gate-notebook` already excludes un-gated keyboard exploration and now redirects software architecture and pure pathway/TF interpretation to concrete skills. Its existing section also keeps the notebook read-only and scoped to decisions. |
| Activation | Dropped: `sciagent activate base science-architect` documents the retiring activation mechanism. |
| Composes with | The `science-architect` row in `docs/lanes.md` records composition with `base` and the scientific lane that owns the work. Skill-level companions `decision-gate-notebook` and `interactive-breakpoint-explorer` were added to `delegate-cli` `## See also`. |

The YAML description maps to the lane row and the three skill descriptions. The empty `agents:`
provenance and four command attributions are role-status data and are dropped; the commands and
skills remain as leaf artifacts.

## `multiome`

| Header section | Destination |
|---|---|
| Purpose | The `multiome` row in `docs/lanes.md` records paired shared-barcode RNA+ATAC, joint embedding/WNN, consensus peaks, motifs, peak-to-gene links, and tracks. |
| Use when | `cellranger-arc-multiome` routes 10x ARC primary processing; `muon-multimodal-analysis`, `signac-chromatin-analysis`, `seurat-multimodal-analysis`, and `scvi-multivi` route joint embeddings/WNN; `peak-atlas-multiome` routes paired consensus peaks; `pycistarget-motif-enrichment` and `pyranges-peak-gene-linkage` route the motif/linkage endpoints. |
| Do NOT use when | `peak-atlas-multiome`, `pyranges-peak-gene-linkage`, and `scvi-multivi` already redirect unpaired data. `cellranger-arc-multiome` already excludes scATAC-only and scRNA-only data. `muon-multimodal-analysis` now carries all three lane boundaries and points to `scglue-unpaired-multiomics-integration`, `snapatac2-atac-preprocessing`, and `scanpy`/scvi skills. |
| Activate as | Dropped: `si activate base multiome` is retired activation syntax. |

The YAML description maps to the lane row. Its leaf provenance remains in the catalog; the overlay
membership is dropped.

## `pathway-signature`

| Header section | Destination |
|---|---|
| Purpose | The `pathway-signature` row in `docs/lanes.md` records pseudobulk functional interpretation through GSEA, TF/pathway activity, signature search, and interactive dashboards. |
| Use when | The triggers are already explicit in `bulk-rnaseq-gsea` (multi-database/custom GSEA), `bulk-rnaseq-activity-inference` (CollecTRI/PROGENy), `coresh-signature-search` (public-data coregulation), and `bulk-rnaseq-pathway-explorer` (`master_unified.csv` to shareable HTML). |
| Do NOT use when | `bulk-rnaseq-gsea`, the lane router, now excludes upstream scRNA-seq work, chromatin/accessibility analysis, and RNA velocity, with concrete skill redirects. Existing `## When not to use` sections on all four interpretation skills retain their narrower method boundaries. |

The YAML description maps to the lane row. Agent, skill, and command provenance membership is
dropped while every leaf artifact remains available.

## `rnaseq-fastq-preprocessing`

| Header section | Destination |
|---|---|
| Purpose | The `rnaseq-fastq-preprocessing` row in `docs/lanes.md` records the FASTQ-to-BAM/counts lane and the TE-compatible Random-One/grouped-SAF integer-count path. |
| Use when | `nfcore-rnaseq-execution` already routes FASTQs, samplesheet/container/disk operations, STAR BAMs, and gene counts; `star-te-preprocessing` routes TE-preserving alignment; `te-reference-saf-build` routes shared reference construction; `te-gene-featurecounts` routes gene+TE matrices. |
| Do NOT use when | `nfcore-rnaseq-execution` now covers downstream annotation/DE/enrichment/interpretation, 10x velocity, and chromatin/ATAC exclusions. The existing sections in `star-te-preprocessing`, `te-gene-featurecounts`, and `annotate-bulk-rnaseq-data` preserve the finer handoff boundaries. |
| Stack | Dropped: `sciagent activate base rnaseq-fastq-preprocessing` documents the retiring activation mechanism. The workflow composition survives in `docs/lanes.md`. |

The YAML description maps to the lane row. The leaf provenance remains in the flat catalog; overlay
membership is dropped.

## `scatac-regulatory`

| Header section | Destination |
|---|---|
| Purpose | The `scatac-regulatory` row in `docs/lanes.md` records accessibility-first differential analysis, CREs, motif activity/enrichment, footprinting, linkage, and ATAC-driven GRNs. |
| Use when | The relevant descriptions already route these situations: `crescendo-scatac-cre-analysis` for sub-peak CREs, `chromvar-motif-accessibility` for motif activity, `tf-footprint-differential-analysis` for occupancy, `pycistopic-atac-topic-modeling` and `pycistarget-motif-enrichment` for topics/motifs, and `scenic-grn-inference` for regulatory networks. |
| Do NOT use when | `snapatac2-atac-preprocessing`, the lane's preprocessing entry point, now excludes scRNA-seq QC/normalization/clustering and RNA velocity with concrete redirects. Existing method-specific sections on the regulatory leaf skills remain in place. |

The YAML description maps to the lane row. Its agent, skill, and command provenance membership is
dropped while the leaf artifacts remain.

## Deliberately dropped content

| Content | Reason |
|---|---|
| Three activation/stack snippets | They teach `sciagent activate base <overlay>`, the mechanism being retired in steps 4b/5. |
| Role/overlay attribution and status-reporting language | It describes provenance ownership and effective-stack reporting, both properties of the retiring role mechanism. |
| Seven sets of `skills:`, `agents:`, and `commands:` membership | Flat catalog discovery routes directly on leaf descriptions. The files and their content survive; the role-to-leaf association has no reader-facing behavior after role deletion. |
| Lane-choice wording such as “use base” | Concrete skill redirects replace role selection. Lane-to-lane workflow relationships remain in `docs/lanes.md`. |

No Purpose fact, user-situation trigger, durable usage boundary, pipeline description, or composition
relationship was dropped.
