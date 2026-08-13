# Compatibility declaration rescind map

Owner ruling 1 on 2026-08-12 rescinded the `compatibility:` frontmatter mechanism and its coupling
drift guard. This map records every declaration that was removed and the reader-facing disposition
of its content. Project paths already used throughout a skill remain part of that skill's workflow
contract; genuine module, version, helper, and companion-skill requirements are stated in body
prose.

| Skill | Removed declaration | Disposition |
|---|---|---|
| `annotate-bulk-rnaseq-data` | `sciagent-scaffold: 00_data/metadata/, 00_data/processed/, 02_analysis/config/analysis_config.yaml, 03_results/annotated_outputs/; external-module: RNAseq-toolkit, TE-RNAseq-toolkit` | `## Prerequisites` states RNAseq-toolkit v0.2.0 and TE-RNAseq-toolkit v0.1.0 with their helper paths. The project paths remain in the routing and reference workflows; the declaration carried no additional reader fact about them. |
| `bulk-rnaseq-activity-inference` | `sciagent-scaffold: 02_analysis/config/config.R, 02_analysis/config/color_config.R, 03_results/checkpoints/, 03_results/tables/, 03_results/plots/` | These are workflow inputs and outputs already shown in `## Quick Start`, `## Pipeline Architecture`, and the output table. The declaration carried no additional reader fact. |
| `bulk-rnaseq-gsea` | `sciagent-scaffold: 02_analysis/config/config.R, 02_analysis/config/color_config.R, 02_analysis/config/pipeline.yaml, 02_analysis/helpers/normalize_gsea.R, 03_results/checkpoints/, 03_results/tables/, 03_results/plots/GSEA/, 03_results/interactive/; external-module: RNAseq-toolkit` | `## Prerequisites` now states the project-vendored RNAseq-toolkit and the required layout; `## Resources` retains its workflow-doc path. The individual paths remain in Quick Start, pipeline, and reference prose. |
| `bulk-rnaseq-pathway-explorer` | `sciagent-scaffold: 03_results/tables/master_unified.csv, 03_results/interactive/; external-module: pathway-explorer` | New `## Prerequisites` states pathway-explorer v2.0.0, the input table, and the output directory. |
| `cellranger-multi-to-anndata` | `sciagent-scaffold: 00_data/raw/, 03_results/checkpoints/` | Both paths remain explicit in `## Quick Start`, `## Standard Workflow`, and the metadata decision pause. The declaration carried no additional reader fact. |
| `consensus-nmf-multirun` | `sciagent-scaffold: 03_results/checkpoints/, 03_results/cnmf/, 03_results/tables/, 03_results/plots/` | The complete path contract remains in `## Quick Start`, the seven-stage workflow, and its output descriptions. The declaration carried no additional reader fact. |
| `decision-gate-notebook` | `sciagent-scaffold: 02_analysis/config/analysis_config.yaml, 02_analysis/helpers/figure_style.R, 02_analysis/notebooks/, 03_results/` | The config, figure shim, notebook location, and read-only results contract remain explicit in `## The three parts`, `## Quick Start`, and verification. The declaration carried no additional reader fact. |
| `delegate-cli` | `sciagent-scaffold: docs/_internal/reasoning/` | `## Fan-out orchestration defaults` already directs stateful research output to `docs/_internal/reasoning/`. The declaration carried no additional reader fact. |
| `figure-style` | `sciagent-toolkit: figure-style; sciagent-scaffold: 02_analysis/helpers/figure_style.R, 02_analysis/helpers/figure_style.py, 02_analysis/config/analysis_config.yaml, 03_results/` | `## Prerequisites — use the project shim` states the helper and contract-library relationship; the following config and results-placement sections state the remaining paths. |
| `interactive-breakpoint-explorer` | `sciagent-toolkit: interactive-style; sciagent-scaffold: 02_analysis/config/analysis_config.yaml, 02_analysis/helpers/interactive_style.py, 02_analysis/notebooks/, 02_analysis/stages/export_explorers.py, 03_results/interactive/, 03_results/objects/` | New `## Prerequisites` states the mounted contract library, shim, config, committed-stage input, and interactive output. Quick Start retains the notebook and exporter paths. |
| `iterative-peak-merging` | `sciagent-scaffold: 01_scripts/R_scripts/createIterativeOverlapPeakSet.R` | `## Local Script Location`, `## Requirements`, and every command example already state the required project-local script. The declaration carried no additional reader fact. |
| `peak-atlas-multiome` | `sibling-skill: peak-atlas-framework` | New `## Prerequisites` names `peak-atlas-framework` as the required companion and links its references, scripts, and checks. |
| `reasoning-trace` | `sciagent-scaffold: docs/_internal/reasoning/, docs/_internal/research/, 02_analysis/stages/, 03_results/` | `## Where traces live`, `## No-ephemeral discipline`, and `## Done when` already state all four locations and their relationship. The declaration carried no additional reader fact. |
| `scrna-cxg-host` | `sciagent-scaffold: 03_results/checkpoints/, 03_results/objects/, 03_results/annotation/` | The Phase A and Phase B quick starts, server-path decision, and verification checklist already state all three locations. The declaration carried no additional reader fact. |
| `te-geneset-gsea` | `sciagent-scaffold: 02_analysis/config/analysis_config.yaml, 03_results/objects/, 03_results/master/; external-module: TE-RNAseq-toolkit` | New `## Prerequisites` states TE-RNAseq-toolkit v0.1.0 and the analysis-config layout; Quick Start and advanced workflow retain the objects and master-output paths. |

The three real defects surfaced by the audit remain fixed: the interactive-style shim/import/root
walk, the peak-atlas framework path resolution and fatal missing-framework diagnostic, and their
regression tests. Manifest schema v2, honest `block_write` I/O failure propagation, and verb help
also remain in place.
