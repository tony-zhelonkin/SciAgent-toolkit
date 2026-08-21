# SciAgent-toolkit — Authoritative Map (2026-06-23)

> Trace produced by `explorer-toolkit` (Explore, read-only). Persisted by orchestrator. Target = the canonical toolkit at `/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit`. This is the backbone reference for the synthesis.

## 1. Purpose, Philosophy, Architecture

**Purpose:** `sciagent` is a per-project AI-harness context manager. It manages "role" bundles (skills + sub-agents + slash-commands + output-style) and injects them into a project's `.claude/` and `.agents/` via symlinks, natively visible to Claude Code and Pi.

**Philosophy (README + 2026-06-04 mind-palace doc):** "In an agent-operated repository, the layout, the docs, and the context files are not documentation of the work — they are the interface through which the work is done." Tidiness as a structural property, not a reminder. The 2026-06-04 plan notes craft rules live in 3 uncoordinated homes (per-project AGENTS.md body, docs/guidelines/, scaffold template) and recommends consolidating into a toolkit-owned managed block `<!-- BEGIN SCIAGENT:CRAFT -->` (NOT yet implemented).

**Architecture:** single bash CLI `bin/sciagent`; logic in `lib/sciagent/*.sh` (no deps; grep/awk YAML; jq optional). Stack depth capped at 2 (base + overlay). Dual-track symlinks into `.claude/` + `.agents/`. Mutations tracked in `.sciagent/manifest.json` (v1). Managed block in AGENTS.md delimited by `<!-- BEGIN SCIAGENT:ROLES v1 hash=... -->` with SHA1 drift detection.

**tags.yaml — 21-tag vocabulary:** trajectory, integration, annotation, reference-mapping, de, pathway, grn, factor-analysis, metaprogram, multimodal, chromatin, motif, qc, preprocessing, conversion, viz, report, hosting, architecture, tooling. Drives `sciagent inject --tag <name>`.

## 2. lib/ and bin/ — Machinery

Verbs: `activate`, `deactivate`, `inject`, `eject`, `validate`, `status`, `list`, `new`, `gitignore`.

Modules: block.sh (managed-block r/w/hash), frontmatter.sh, roles.sh, symlinks.sh (dual-track + manifest), skill_deps.sh (transitive requires graph), collisions.sh (cross-namespace name collisions; allowlist), stack.sh (last-wins shadowing + render_block_body), claude_settings.sh (outputStyle apply/revert), activate.sh (Phase A stack walk → B requires closure → B.5 complementary soft-warns → symlinks → block write → manifest), inject.sh (`_injected` synthetic overlay), eject.sh, deactivate.sh, validate.sh (tag-vocab, requires cycle, collision, docs layout), status.sh (`--json`).

Injection: `sciagent inject <name>` auto-detects kind; ambiguous hard-fails; `--tag` pulls all skills with tag. Symlinks: `symlink_create_dual` makes `ln -sfn` in both `.claude/` and `.agents/`. Stack: `stack_walk base [overlay]` last-wins; `render_block_body` renders Skills (effective/injected/inherited), Sub-agents, Slash commands, Output style.

## 3. Skills Inventory (~72)

Notable: `anndata`, `scanpy` (concept; plots umap/dotplot/heatmap/violin), `single-cell-rna-qc`, `anndatar-seurat-scanpy-conversion`, `scrna-pipeline-conventions` (concept; "house style: numbered scripts, 03_results/ layout, multi-checkpoint" — BUT uses a FLAT `03_results/tables/`+`03_results/plots/` layout that CONTRADICTS the stage-based template), `bulk-rnaseq-gsea`, `bulk-rnaseq-activity-inference`, `bulk-rnaseq-pathway-explorer` (HTML dashboards), full scvi-* ecosystem, mofa-*, ATAC/motif/grn, `architecture-first-dev` (concept), `architecture-treemap`, `skill-creator` (meta), `annotate-bulk-rnaseq-data`, `te-geneset-gsea` (newer).

**No dedicated figure-legibility/aesthetics skill. No `viz`-tagged skill** despite the `viz` tag existing. No reasoning-trace/decision-log skill.

Skill frontmatter format:
```yaml
name, description (no angle brackets), license: MIT
metadata: { scope: concept|implementation, requires: [], skill-author, last-reviewed,
  category: analysis|practice|meta|tooling, tier: simple|standard|advanced, version,
  tags: [...from tags.yaml], complementary-skills: [], contraindications: [] }
```
Scope caps linted: concept ≤500 body lines, implementation ≤350.

## 4. Agents Inventory (19)

**analysis-base/ (7, in base role, all sonnet):**
- `captions` — fire-and-forget; writes README.md to `03_results/<stage>/` with publication captions (Script/Function/Config/Input table); traces generating function.
- `doc-curator` — audits C1 (file-size), C2 (one-way reference), C3 (uncaptioned artifacts). Never deletes without confirmation.
- `handoff` — writes dated `YYYY-MM-DD_<slug>.md` to sessions/; **Step 2 scans 03_results/ for uncaptioned artifacts**.
- `bio-interpreter` — web → mechanism; writes dated note to research/.
- `insight-explorer` — stats exploration of RDS/CSV; recommends viz.
- `code-reviewer` — before/after refactor review vs project docs.
- `docs-librarian` — web search for tool docs/params/issues.

**architect/ (12):** bioinf, wetlab, **graphic (Tufte information-graphic reviewer — data-ink, channel economy, perceptual encoding, legibility, affordances)**, stat, divergent (saboteur), ml, mapper (read-only cartographer→map.md), **slicer (opus)**, architect, synth, status-reporter, meta-architect, feature-reviser.

Figure-touching: graphic (architect-role ONLY), captions, doc-curator (C3), handoff (Step 2).

## 5. Commands Inventory

**architect/ (18):** `/map` (cartography→map.md), `/review` (parallel expert panel→review/<name>.md; `--iterate` cross-informed), `/synthesize`, `/design`, `/architect`, **`/plan`** (decompose approved design into phases→plan/README.md + phase-NN.md; each phase ≤3-5 files, testable, ordered data-models-first), **`/implement`** (single-phase default or `--auto` end-to-end; stamps frontmatter status:shipped, commit:, files_touched:, verification:), **`/verify`** (mechanical drift check plan vs code vs ADR; CLEAN/INCOMPLETE/DRIFT/NEEDS REVIEW; frontmatter-first), `/status`, `/diagram`, `/audit-slice`, `/components-extract`, `/synthesize-audit`, `/architecture-treemap`, `/meta-map`, `/meta-design`, `/meta-apply`, **`/meta-plan`** (cross-feature sequencing: dependency graph, portfolio-phase table, collision check).

**universal (1):** `/commit` (atomic, imperative, no AI attribution).

Existing pipeline: `/map`→`/review`→`/synthesize`→`/design`→`/plan`→`/implement`→`/verify`. **Software-architecture oriented** (data models→logic→storage→API). No bio-specific plan template. **No Opus-decompose/Sonnet-implement split at command level. No every-N-phases review cadence. `--auto` runs in one session, not detached background. No `--resume`.**

## 6. Roles, System-Prompts, Templates

Roles: base.yaml (7 agents + 36 skills + commit), architect.yaml (12 agents + 2 skills + 18 commands + output_style architect-mentor), pathway-signature.yaml, rnaseq-fastq-preprocessing.yaml, scatac-regulatory.yaml.

System-prompts: architect-mentor.md (output style for architect role).

**templates/project/analysis/** (what `sciagent new project` writes): AGENTS.md template, CLAUDE.md (1-line @AGENTS.md), analysis_config.yaml (full), `03_results/{01_qc,02_eda}/{tables,figures}/` with README.md.template per stage + objects/ + master/ + interactive/ + _scratch/, docs/_internal/{scientific-context, sessions/README, reasoning/.gitkeep, research/README}, .gitignore-seed.

**CRITICAL:** template uses `03_results/<stage>/{tables,figures}/` with README at `03_results/<stage>/README.md` — matches the user's desired pattern. But `scrna-pipeline-conventions` SKILL.md uses a DIFFERENT flat layout — unresolved contradiction.

Config template viz: dpi 300, width 10, height 8, base_size 12, title 14, label 11. README templates define caption format (`## <filename>` / **Finding** / Script|Function|Config|Input). .gitignore-seed hard-ignores objects/*.{h5ad,rds,...}; **figure/table gitignore rules COMMENTED OUT**; _scratch/* ignored.

## 7. docs/ — Prior Direction

`docs/architecture.md` — canonical design spec (two-role stack, dual symlinks, managed block, manifest schema, Step-0 output-path protocol).
`docs/guidelines/` — visualization.md (theme_publication() R fn, Okabe-Ito, 300 dpi, Nature widths 89/183mm, save_publication_plot(), quality checklist incl. "Font sizes readable at target size", "No overlapping text", "Axis labels include units"), code_style.md, core_architecture.md, data_processing.md, gsea_analysis.md, master_tables.md. **These are SSOT islands — AGENTS.md templates never reference them; base_size=12 is too small; not enforced by any agent.**
`docs/_internal/plans/2026-06-04_mindpalace...` — craft lives in 3 uncoordinated places; rules frozen per-project; recommends SCIAGENT:CRAFT managed block + wire code-reviewer/doc-curator to gate against code_style.md.
`docs/_internal/plans/2026-06-17_annotate-skill-refactor.md` — example of using architect pipeline on the toolkit itself.
`docs/proposals/` — sciagent-extension-design-spec.md; ai-research/01-trace-recording-claude-code.md (trace-recording proposal, not implemented).

## 8. Tests = Contracts

activate/deactivate/inject/eject semantics; managed-block roundtrip/drift/markers; collision allowlist; manifest ownership; no-duplicate-basenames; no-exit-in-libs; new_project_types; **test_results_gitignore_layout** (objects/ hard-ignored, _scratch/* ignored — does NOT test {figures,tables} structure); skill_frontmatter_valid; requires graph integrity; **skill_scope_lint** (concept≤500/impl≤350); stack shadowing; status; **test_tags_vocabulary** (closed vocab); **test_validate_docs_layout** (docs/_internal must be gitignored; .md in 03_results warns; non-standard handoff names warn). **No figure/table gitignore enforcement (deliberately off).**

## Pain Points — Already-Provides vs Missing

| | Already provides | Missing |
|---|---|---|
| **A. Figure legibility** | docs/guidelines/visualization.md (theme_publication, 300dpi, Okabe-Ito, Nature widths, checklist); config template sizes; `graphic` agent (Tufte reviewer, architect-only) | No enforcement path (guidelines never loaded into agent context); no linewidth standard; no label-truncation rule; graphic agent inaccessible in analysis role; no `viz`-tagged skill; no conference-room sizing variant; base_size 12 too small |
| **B. Results placement** | template `03_results/<stage>/{tables,figures}/README.md`; config stage subdirs; AGENTS.md template specifies layout | `scrna-pipeline-conventions` SKILL.md uses contradictory flat layout; no validate check for stage-dir naming; no enforcement that source tables accompany figures |
| **C. README adjacency** | `captions` agent; `doc-curator` C3; `handoff` Step 2; README templates; AGENTS.md "captions written as artifact produced" | No auto-trigger on figure creation; no pre-commit/CI caption check; captioning remains manual; no session-loop mandate |
| **D. Planning decomposition** | `/plan` (≤3-5 file testable phases); `/implement --auto`; `/verify`; `/meta-plan`; phase status machine | No model-tier split (Opus-decompose/Sonnet-implement); no every-N-phase Opus review cadence; no actual runnable+artifacts gate (only anchor grep); no detached background exec; no bio phase template; no `--resume` |
| **E. Reproducibility** | scrna-pipeline-conventions (numbered scripts, multi-checkpoint, config.py, decisions: replay); AGENTS.md rules; handoff; dated research notes; git-SHA phase frontmatter; gsea_seed | No reasoning-trace agent populating reasoning/; no CoT capture; no script-durability audit agent; no session-replay skill; biomart pinning planned not shipped; trace-recording proposal unimplemented |

## Dense Synthesis

Mechanism: bash-only harness-agnostic context manager; `sciagent activate base [overlay]` reads role YAML, resolves requires closure, dual-symlinks into `.claude/`+`.agents/`, stamps AGENTS.md managed block w/ SHA1 drift detection, writes manifest. 72 skills, 19 agents, 19 commands.

The toolkit already ships a competent software-architecture planning pipeline (`/map`→`/review`→`/design`→`/plan`→`/implement`→`/verify`) and a strong caption/curation triad (captions + doc-curator + handoff). The gaps are exactly the user's pain points: (A) figure legibility exists only as dead-end guidelines + an architect-only `graphic` agent — never injected, never enforced, too-small defaults; (B) the correct stage-based results layout is in the template but contradicted by a skill and not validated; (C) captioning is manual, no trigger/hook; (D) the planning pipeline has no model-tiering, no review cadence, no real artifact gate, no detached execution; (E) reasoning-trace persistence is proposed but unimplemented. The fix surface: a `viz`/figure-style skill + pre-save checklist agent; reconcile layout SSOT + add validate check; auto-trigger captions via hook; an Opus-decompose/Sonnet-implement/Opus-review orchestration command set; a reasoning-trace agent + anti-ephemeral rule. Hooks (deterministic) are the right enforcement substrate for C and E.
