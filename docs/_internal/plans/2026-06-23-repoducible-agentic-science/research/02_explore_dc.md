# Findings: DC_mouse_cancer Repo Exploration

> Trace produced 2026-06-23 by the `explorer-dc` (Explore) agent. Read-only agent; persisted by orchestrator from returned findings.

**Target:** `/scratch/current/antonz/projects/DC-nexus/DC_mouse_cancer` (single-cell: scanpy/scVI, mm39, ~85K cells, 5 pipeline stages)

## Verdict Table — Five Concerns

| # | Concern | Status | Evidence |
|---|---------|--------|----------|
| 1 | FIGURE LEGIBILITY | **Documented & enforced** | AGENTS.md §"Figure standards"; `plot_utils.py` `set_paper_style()` + `save_figure()` + `POINT_SIZE` from config; `analysis_config.yaml:visualization` fully parameterized |
| 2 | RESULTS PLACEMENT | **Documented & enforced** | AGENTS.md schema; `config.py` `PATHS.stage_dir()`; strict `.gitignore`; `objects/` checkpoints |
| 3 | README ADJACENCY | **Documented & enforced (one gap)** | Every stage 00–05 has README.md; caption format specified; `interactive/`+`master/` lack READMEs |
| 4 | PLANNING DECOMPOSITION | **Followed but NOT in docs/_internal/plans/** | Plans landed in `reasoning/` (ADR-style) + `sessions/` (handoffs); `plans/` only `.gitkeep`. Date-slug/Opus-decompose convention NOT enforced; orchestration ad-hoc Sonnet+Opus pairs per stage |
| 5 | REPRODUCIBILITY / NO EPHEMERAL | **Documented & enforced** | All scripts durable, logged, checkpointed; `requirements-optional.txt`; `logs/rerun/` full pipeline logs; `.pid`/`.flag` sentinels |

## Top-level Layout

`00_data/{raw,processed,references}` (read-only) · `01_modules/{SciAgent-toolkit,RNAseq-toolkit,dev-env}` · `02_analysis/{config,helpers,scripts,notebooks}` · `03_results/{00_build,01_qc,02_eda,03_integrated,04_convert,05_annotation,objects,master,interactive,_scratch}` · `.agents/`+`.claude/` (symlinks → toolkit) · `docs/_internal/` (gitignored) · `logs/`.

`.sciagent/manifest.json`: stack `["base","_injected"]`, 38 skill symlinks, 7 agent symlinks, 1 injected skill (`mllmcelltype-consensus-annotation`).

## AGENTS.md / CLAUDE.md

CLAUDE.md = `@AGENTS.md` (11 bytes). Umbrella CLAUDE.md at `/scratch/current/antonz/projects/DC-nexus/CLAUDE.md` explains submodule vs umbrella model.

AGENTS.md (165 lines incl. injected block):
- **Critical rules:** Config not hardcoding; Normalize then visualize (checkpoint under objects/); Read-only input; Cache expensive ops.
- **Figure standards (54–68):** "Figures are judged as they will be seen: shrunk in a journal column or projected to a conference auditorium. Favour bigger, fewer, bolder marks." One saver two files (`save_figure()` → png+pdf); uniform dots/fonts from config; equal embeddings (pca-based/scvi-based parallel subdirs); cross-dataset parity (colon/ln share stem, colours, ordering).
- **Artifact captions (70–76):** every file in `03_results/<stage>/` needs caption in README.md: `## <filename>`, 1-sentence finding, then `Script | Function | Config | Input` table.
- **Documentation namespace (78–97):** session handoff → sessions/; research note → research/; decision log → reasoning/; public phased plan → docs/plan/.

ABSENT: no `docs/_internal/plans/{date-slug}/` planning decomposition pattern.

## .claude/ and .agents/

Symlink farms into toolkit. Agents (7): bio-interpreter, captions, code-reviewer, doc-curator, docs-librarian, handoff, insight-explorer. Commands (1): commit. Skills (38): full scRNA/scATAC/scVI stack + bulk-RNA + metabolomics + injected mllmcelltype-consensus-annotation. Several installed-but-never-invoked (shinymultiome-uio-host, treearches, scenic-grn-inference) — base role installs regardless of per-project scope.

## 02_analysis — Compute/Viz

Explicit `NN_<slug>.py` (compute) / `NN_<slug>_viz.py` (viz). Stage 05 annotation: 05_markers.py, 05_annotate.py, 05_reconcile.py, 05_reconcile_rds_mapback.R, 05_annotate_viz.py.

Centralized `analysis_config.yaml` (~45KB, 640+ lines) consumed by Python (`config.py` → CONFIG/PATHS/PARAMS/OBJECTS) and R (yaml::read_yaml). `PATHS.stage_dir()` resolves all paths; no hardcoding. Helpers: plot_utils, config, anndata_utils, metadata_utils, markers, qc_utils. `requirements-optional.txt` documents non-installable packages + workarounds (reproducibility artifact). `annotation_profiles/` = second config layer for LLM annotation.

## 03_results — Stages

| Stage | figures/ | tables/ | README.md | Notes |
|-------|----------|---------|-----------|-------|
| 00_build | empty | 2 CSVs | yes | captioned |
| 01_qc | 16 png+pdf | yes | yes | 24-artifact captions |
| 02_eda | pca-based/+scvi-based/ | 4 CSVs | yes | embedding-parity note |
| 03_integrated | yes | yes | yes | |
| 04_convert | none | yes | yes | README notes no figures by design |
| 05_annotation | pca-based/{umap,dotplot}/+scvi-based/... | 9 CSVs | yes (long) | passA/passB/reconcile/_preview subdirs |
| interactive | — | — | no | empty |
| master | — | — | no | accumulator |
| objects | — | h5ad+rds (gitignored) | no | checkpoints incl. scvi_{colon,ln,integrated}/ |

Source-table adjacency followed; LLM trace dirs (passA/passB/trace) kept inside 05_annotation. READMEs detailed (function name, config keys w/ values).

## logs/ Pattern (reveals actual workflow)

Files: `02_embed_colon.log` vs `_viz.log` (compute/viz separate processes); `rerun/MASTER_*.log` (chained 12-stage background pipeline, triggered by UMAP rotation-drift bug from unseeded umap, fixed + rerun); `idem_rc*_ln.log` (RECLUSTER_ONLY mode re-runs Leiden on frozen embeddings, verifies obsm MD5 hashes); `rename_*_viz.log` (viz re-run after stem renames WITHOUT recompute); `05_annotate_execute.log` (0 bytes = dry run) vs full annotation logs (cost-discovery before paid LLM calls); `.pid`/`.flag` sentinels. **No throwaway scripts** — all log invocations map to files in scripts/.

## docs/_internal

README (namespace + naming), scientific-context.md, handoffs/ (dated kickoffs), sessions/ (8 dated handoffs), reasoning/ (ADR-style decision logs incl. `2026-06-03_object-build-plan.md` = "Decision log & orchestration plan" with Wave A/B/C/D table + D1–D10 decisions; annotation strategy/normalization/palette/reconciliation), research/ (marker reference), plans/ (EMPTY), reports/ (EMPTY), personal-notes/, .ref/_scratch/ (archived old scripts).

Effective planning pattern: reasoning note → per-session handoffs. The closest plan (`2026-06-03_object-build-plan.md`) states "Orchestration — Sonnet implementer + Opus reviewer per stage/object" and "Opus audits code AND the produced artifacts" — closest existing match to the 5-concern model, but NOT in plans/, no phase-sizing or 2-3-cadence.

## Figure Theme Code (plot_utils.py — single SSOT)

- `set_paper_style()` — sets rcParams from CONFIG['visualization'] (font.size, axes.titlesize bold, labelsize bold, tick sizes, legend; removes top/right spines; `pdf.fonttype=42` Illustrator-editable; savefig.dpi 300).
- `save_figure()` — png@300 + pdf@600, extension-agnostic, closes fig.
- `POINT_SIZE` = 12.0 single constant.
- `embedding_dir()` — pca-based/scvi-based from config.
- `rasterize_axes()`, `purge_figures()` (deletes stale {prefix}*.png/pdf at script start — prevents orphans), `add_centroid_labels()` (white-halo), `annotation_palette()` (deterministic, loud warn on unmapped), `grouped_dotplot()` (data-driven sizing).

Config viz block: `fig_save_dpi 300, fig_rasterized_dpi 600, umap_point_size 12, base_size 13, title_size 15, label_size 12, tick_size 11, annotation_onplot_labels true, onplot_label_size 9, halo_width 2.5, embedding_subdirs {pca: pca-based, scvi: scvi-based}`. Comment in config: "# Legibility for print + conference slides (read from a distance)". Gap: no explicit axis-truncation guard in set_paper_style.

## Single-cell-specific Deviations

1. Dual-basis (pca-based/scvi-based) parallel figure subdirs — unique; enforced by `embedding_dir()`.
2. Checkpoints are h5ad/rds + scVI model dirs in objects/ (gitignored, never deleted).
3. LLM-assisted annotation (mllmcelltype) with passA/passB/trace audit trail (discussion JSON, token counts, git HEAD) preserved as reproducibility artifact.
4. Config.yaml serves Python + R (h5ad→Seurat conversion stage).

## Dense Summary

DC is the strongest repo for figure legibility + results placement (both documented in AGENTS.md and mechanically enforced via plot_utils.py + config.py, png+pdf always, purge_figures prevents stale, rasterized dense layers). README adjacency documented + followed in all numbered stages (minor: interactive/master lack READMEs). Biggest gap = planning decomposition (plans/ empty; planning in reasoning/ ADRs + handoffs; Sonnet+Opus pairs per stage but no phase-sizing/2-3 cadence). Reproducibility excellent (full chained reruns in logs/rerun/, RECLUSTER_ONLY obsm-hash idempotency checks, UMAP drift incident root-caused + rerun, requirements-optional.txt). Logs reveal compute/viz strictly separate; viz re-run after renames without recompute; dry-run cost discovery before paid LLM. SC adaptations: dual-basis subdirs, LLM annotation audit trail, Python+R config bridge.
