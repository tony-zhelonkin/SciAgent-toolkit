# Findings: STING-cGAS-GSE329522 Repo Exploration

> Trace produced 2026-06-23 by the `explorer-sting` (Explore) agent. The agent runs read-only; this file was persisted by the orchestrator from the agent's returned findings.

## Verdict Table — 5 Concerns

| Concern | Status | Evidence |
|---|---|---|
| 1. FIGURE LEGIBILITY | (b) documented but not fully enforced | AGENTS.md documents multi-panel prohibition and one-claim-per-figure rule, but no explicit font-size floor, line-width floor, conference-room/journal-column legibility spec, or no-truncated-label mandate. `base_size=11` in primary viz script; `custom_minimal_theme` uses `base_size=12`. No systematic constraint preventing small labels. |
| 2. RESULTS PLACEMENT | (a) documented & enforced (mostly) | AGENTS.md §"Normalize then visualize" mandates `03_results/<stage>/{figures,tables}`, objects to `03_results/objects/`, master to `03_results/master/`. `config.R` `stage_dir()` enforces this programmatically. Actual files match. One gap: `03_results/03_de/` has no README.md despite AGENTS.md requirement. |
| 3. README ADJACENCY | (b) documented but not enforced | AGENTS.md line 71: "Each `03_results/<phase>/` carries a `README.md`" — explicit mandate. 01_qc, 02_eda, 04_tf all have README.md. **03_de is missing its README.md**. No automated check enforces this. |
| 4. PLANNING DECOMPOSITION | (d) absent | `docs/_internal/plans/` exists but is **empty** (only `.gitkeep`). No date-slug plan directories with phase decomposition. No Opus/Sonnet split documented. No plan artifact format defined. |
| 5. REPRODUCIBILITY / NO EPHEMERAL SCRIPTS | (c) followed by convention but undocumented | All outputs in durable numbered scripts. `03_results/_scratch/` empty. Reasoning artifacts in `docs/_internal/reasoning/` & `sessions/`. But no explicit policy against ephemeral scripts in AGENTS.md; `docs/_internal/_scratch/notes.md` is the kind of chat-reasoning the toolkit wants to avoid. |

## Top-Level Layout

| Directory | Role |
|---|---|
| `00_data/` | raw/, processed/, references/ — input data |
| `01_modules/` | Git submodules: `RNAseq-toolkit/`, `SciAgent-toolkit/`, `.ref/` (gitignored) |
| `02_analysis/` | scripts/, config/, helpers/ (empty), notebooks/ (empty) |
| `03_results/` | stage dirs (01_qc, 02_eda, 03_de, 04_tf) + objects/, master/, interactive/, _scratch/, 99_fake/ |
| `.agents/`, `.claude/` | symlink farms into SciAgent-toolkit (agents/, commands/, skills/) |
| `.sciagent/manifest.json` | symlink registry; stack `["base","pathway-signature"]`, 97 symlinks, `injected: []` |
| `AGENTS.md` | primary agent instruction doc (~5KB) |
| `CLAUDE.md` | 11 bytes: `@AGENTS.md` redirect |
| `docs/` | docs/_internal/ (reasoning, plans(empty), sessions, research, inbox, reply-package), docs/{plan,stages,reference}/ (empty) |

## Agent Instruction Surface (AGENTS.md)

Effective instruction file is AGENTS.md (CLAUDE.md forwards). Contains (1) auto-generated SCIAGENT:ROLES manifest block (stack base + pathway-signature; 44 skills, 7 sub-agents, 1 command), and (2) project-specific "Analysis conventions" block.

**Figure & conclusion discipline (lines 68–73):**
- "One claim = one dedicated, captioned figure. No conclusions that live only in chat/handoff."
- "Multi-panel ONLY when sub-panels share the same axis and are directly comparable."
- "Every figure documents how it was generated — auditable, never a black box."
- "Each `03_results/<phase>/` carries a `README.md` captioning, laconically: (1) the STATEMENT the artifacts make, and (2) the MECHANISM behind that statement."

ABSENT from AGENTS.md: font-size floors, line-width minimums, journal-column/conference-room legibility requirements, no-truncated-axis-label rule, ink-to-info ratio.

**Results placement (lines 75–77):** explicit COMPUTE script `NN_<name>.R` (no ggplot/ggsave) + VIZ script `NN_<name>_viz.R` (no statistics). "Viz must run standalone after compute and must never recompute statistics."

**Planning decomposition:** absent. **Reproducibility:** partially present ("no conclusions in chat/handoff") but no explicit anti-ephemeral policy.

## 02_analysis — Scripts and Config

Compute+viz pairs throughout: `00_setup_metadata.R`, `01_mapping_qc.R`+`_viz.R`, `02_de_limma_trend.R`+`_viz.R`, `03_decoupler_tf.R`+`_viz.R`, plus 03b–03g pairs. Helpers/notebooks empty.

Centralized config:
- `analysis_config.yaml` — single YAML: metadata, design, thresholds, palettes, viz settings (`width=10, height=8, width_narrow=6, width_wide=14, dpi=300, base_size=12`), stage registry, path roots, schemas.
- `config.R` — loads YAML; exposes `DIR_RESULTS/OBJECTS/MASTER/SCRATCH`, `stage_dir()`, `load_or_compute()`, color constants (`DIVERGING_COLORS`, `GROUP_COLORS`, `AXIS_COLORS`, `MODULE_COLORS`), gene vectors, `provisional_caption()`.
- `env/` — key_packages.txt, package_manifest.csv, sessionInfo.txt (reproducible env snapshot).

Scripts pick output via `stage_dir("04_tf","figures")` reading YAML. No hardcoded paths.

## 03_results — Stage Folders

| Stage | figures/ | tables/ | README.md | Source-table-adjacent? |
|---|---|---|---|---|
| `01_qc/` | 4 PDFs | 4 CSVs | YES | YES |
| `02_eda/` | 1 PDF | 2 CSVs | YES (template only — unfilled) | YES |
| `03_de/` | 4 PDFs | 2 CSVs | **ABSENT** | YES |
| `04_tf/` | 18 PDFs + deck_assets/ | 30 CSVs | YES (743 lines, exemplary) | YES |

Special dirs: `_scratch/` empty (gitkeep); `99_fake/` empty; `interactive/` empty; `master/` cross-stage tables; `objects/` checkpoint RDS.

04_tf README is exemplary — per-figure STATEMENT + MECHANISM + script/function/input table for all 18 figures. 02_eda README is unfilled template. 03_de README missing entirely.

## docs/_internal

`README.md` defines layout: scientific-context.md, reasoning/, sessions/, research/. One-way reference rule (internal may reference public; public must not link back).

- `plans/` — EMPTY (gitkeep). No plan format/example.
- `reasoning/` — 5 dated bundles (e.g. 2026-06-08-hif1a-rank-instability with 00-SYNTHESIS/01-method/02-validation/03-adversarial-redteam; figure-audit; deck-review; serendipity-resolution). SYNTHESIS + advocate + adversary structure.
- `sessions/` — README only.
- `research/` — 2 dated notes.
- `inbox/` — email thread + reply-reasoning/ (11 numbered deliberation files).
- `reply-package/` — dated packages with build_narrative_pptx.py + pptx + REPLY_DRAFT.md.
- `_scratch/` — preview PNG subdirs (previews, previews_adv, previews_audit, previews_impl, previews_fix, previews_tript) + notes.md (unstructured scratchpad).

## Figure Styling / Theme Code

Primary viz theme (`03_decoupler_tf_viz.R:75-79`): `theme_minimal(base_size=11)` + titles 12pt, subtitles 9pt, captions 7pt. Linewidths 0.7 (primary), 0.3–0.6 (secondary). Toolkit theme `custom_minimal_theme.R` = `theme_classic(base_size=12)`. Config YAML viz block: base_size 12, title 14, label 11, width_narrow 6, width_wide 14, dpi 300. `save_fig()` writes stamped PDF + optional 300-dpi PNG. Colors centralized in config.R (Okabe-Ito based). No legibility mandate ties sizes to floors; no truncation guard.

## Ephemeral vs Durable

Durable: numbered committed R scripts; `load_or_compute()` checkpointing; structured dated reasoning bundles; reply-package build scripts. Quasi-ephemeral: `docs/_internal/_scratch/notes.md` (unstructured); preview PNG iteration dirs. `03_results/_scratch/` empty (policy honored).

## Surprises

- CLAUDE.md is 11 bytes (pure redirect).
- `omnipathr-log/` at repo root holds 5 runtime logs NOT gitignored (`.gitignore` only covers `logs/*.log`).
- `.venv_pkgs/` at repo root not gitignored (only `.venv/` is).
- `01_modules/.ref/` convention for reference codebases (gitignored).

## Summary for Orchestrator

Strongest: rigorous compute/viz split (12 matched pairs, `stage_dir()` routing, viz never recomputes); centralized `analysis_config.yaml`+`config.R`; durable dated reasoning bundles with adversarial structure; single-sourced color palettes.

Biggest gaps: (1) figure legibility has no font/linewidth floor, no conference-room spec, no truncation rule (base_size 11, 7pt captions would fail projection); (2) results placement minor (03_de README missing; stage_dir warns not blocks); (3) README adjacency not enforced (03_de missing, 02_eda unfilled); (4) planning decomposition completely absent (empty plans/); (5) reproducibility followed but undocumented (no explicit anti-ephemeral policy).
