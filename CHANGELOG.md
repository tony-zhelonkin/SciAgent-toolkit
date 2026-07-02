# Changelog

All notable changes to the SciAgent-toolkit project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

Reproducibility: portable, pin-respecting symlinks + a toolkit-locality guard,
so a dropped-in agent activates against the copy the project actually ships.

### Fixed
- **`sciagent activate` now reliably materializes `.claude/statusline.sh` (chmod +x) and `.claude/settings.json` (statusLine + gold defaults: `editorMode: vim`, `effortLevel: xhigh`, `alwaysThinkingEnabled`, `autoMemoryEnabled: false`, hooks) for EVERY project, not just ones scaffolded via `sciagent new project`.** Root cause: those files were only ever rendered from `templates/project/_common/.claude/{settings.json,statusline.sh}.template` inside `new.sh`'s `_render_tree`, which runs once at scaffold time; `activate.sh`'s `cmd_activate` never touched them at all, so any project that vendored the toolkit without running `new project` (the common case) never got a status line. New `claude_settings_ensure_statusline`/`claude_settings_ensure_project_defaults` (`lib/sciagent/claude_settings.sh`) are now called unconditionally from `cmd_activate` (and therefore from `sciagent update`, which re-runs activate): statusline.sh is written only if absent and the executable bit is always re-asserted; settings.json is written verbatim if absent, or — if present — has any *missing* top-level keys backfilled from the template via a reverse `jq` merge that always lets the existing file's values win, so a user's own settings are never clobbered. `new.sh`'s `_render_tree` also lost the executable bit when rendering `statusline.sh.template` (plain `sed > dst`, no mode copy) — `_new_project` now `chmod +x`s the rendered file.
- **Skill/agent/command symlinks are now RELATIVE when the toolkit lives inside the project tree.** `symlink_create_dual` (used by `activate`) and the `inject` symlink sites previously wrote every `.claude/*`/`.agents/*` link as an absolute path into whatever `$SCIAGENT_TOOLKIT` happened to be at activation time — breaking portability across host/container/machine and escaping the submodule pin. A new `symlink_target_for` helper (mirroring `symlink_create_helper_lib`) relativizes each link to its own directory via `realpath -ms --relative-to` when the toolkit is at/below the activation CWD, and falls back to the absolute path for external/global checkouts (where a relative `../../../…` chain would be fragile). Applies to all three namespaces (skills, agents, commands) and both mirrors (`.claude` and `.agents`); the manifest still records link *paths*, so `deactivate` teardown round-trips unchanged.

### Added
- **Toolkit-locality guard for mutating verbs.** `activate`/`inject`/`eject` now refuse to run when the project ships its own `./01_modules/SciAgent-toolkit` but the resolved `$SCIAGENT_TOOLKIT` is a *different* (external/global) toolkit — the footgun where a global `sciagent` silently symlinks a project against, and pins it to, the wrong copy. In-repo vs external is decided by real-path equality against `./01_modules/SciAgent-toolkit`. Override with the `--allow-external-toolkit` flag or `SCIAGENT_ALLOW_EXTERNAL_TOOLKIT=1`. `deactivate` is deliberately unguarded (it only removes what its own manifest owns). New tests: `tests/test_relative_symlinks.sh`, `tests/test_toolkit_locality_guard.sh`.

Unified figure style: the figure-style contract drops the dual print/screen
variant model in favor of ONE legible tier emitted in BOTH formats. Promoted
from the in-project Wave-0 prototype that proved it across a full results tree.

### Changed
- **`lib/figure-style/figure_helpers.R` — unified single-variant, dual-format contract.** `project_theme()` is now ONE legible tier (no per-variant font floors; plain non-bold axis titles, bold title/legend-title/strip, right legend with inter-row air, `axis.line` 0.4, `plot.title.position` "plot", config-driven margins) read from the `figures:` block. `save_figure()` emits exactly `<name>.pdf` (cairo, Unicode glyphs) + `<name>.png` (no `.print`/`.screen` suffix), shares ONE geometry, purges stale same-stem files, and NEVER re-themes (the caller owns all theming); `name` may carry a subdir. `save_overview()` keys its README caption on `<name>.png`. `style_series()` is the alignment-safe running-sum normalizer (via the `grs_restyle` closure: ES-clamp, single top-justified right legend, bottom-only xticks, config `running_sum_heights`), with `style_running_sum` as an alias. The `variant` parameter is retained on `project_theme`/`set_paper_style`/`save_figure` but IGNORED (drop-in compat); `void=TRUE` borderless-panel behavior preserved.
- **`lib/figure-style/figure_helpers.py` — contract parity.** `set_paper_style()`/`project_theme()` single tier (plain axis titles, `variant` accepted-but-ignored); `save_figure()` emits `<name>.pdf` + `<name>.png` (no suffix), one shared geometry (`wide`/`width`/`height` overrides), purges stale same-stem files; `save_overview()` keys the caption on `<name>.png`. Removed the per-variant font-bump/geometry/normalize helpers.
- **`analysis_config.yaml:figures` template + `figure_style.{R,py}` shim templates** rewritten to the unified keys (`base_size` 14, `title_size`, `subtitle_size`, `axis_title_size`, `axis_text_size`, `strip_size`, `legend_text_size`, `caption_size`, `label_size`, `cue_size`, `line_width`, `point_size`, `width`/`height`, `width_wide`/`width_narrow`, `dpi`, `formats: [pdf, png]`, `top_n`, `running_sum_ylim`/`running_sum_top`/`running_sum_heights`, `nes_cap`, `caption_wrap_column`, sub-layout dirs). Dropped `base_size_column`/`width_column`/`height_column`/`variants`.
- **`validate --check figure-style`** base-font floor lowered 16 → 14 to match the unified single tier (the new template default no longer trips the toolkit's own guardrail).
- **`skills/figure-style/SKILL.md`** rewritten for the unified single-variant, dual-format style (removed `.print`/`.screen` language).
- **`craft.yaml` (SCIAGENT:CRAFT single source)** — the always-on Figures craft standard rendered into every repo's `AGENTS.md` now states ONE legible tier: `figure_base_size` floor 16 → 14, the `figure_print_base_size` floor removed, and the body line reworded to "one legible tier, base >= 14pt; emit a vector PDF + raster PNG from one plot object" (no print/screen split). Propagated via `sciagent update`.
- **Craft / contract prose swept to the unified model** — the `figure-audit` and `captions` agents, the `add-figure-variant` and `interpret-storm` commands, the plan + project `AGENTS.md` templates, `docs/guidelines/visualization.md`, and the `base` role comment now reference `<stem>.pdf` + `<stem>.png` (no `.print`/`.screen`), a single 14 pt tier, and `save_figure()`/`save_overview()` (dropped `save_figure(variant="both")`, `base_size_column`, dual-variant language).
- **`tests/test_validate_figure_style.sh`** conformant fixture uses the unified `<name>.png` naming.

### Added
- **Okabe-Ito palette helpers**: `scale_color_okabe()`/`scale_fill_okabe()` (R) and `okabe_palette()` (Python), sourced from `colors.okabe_ito` with a canonical 8-colour fallback.

## [3.3.0] - 2026-06-24

Skill lifecycle: a natural, judgment-driven path from `experimental` through
`stable` and `deprecated` to a reference-only attic — surfaced softly, never
enforced. No hooks, no fail-closed `validate` check; the toolkit only shows the
state and leaves the call to a human during practice.

### Added
- **`metadata.status:` skill field** (`experimental | stable | deprecated`, default `stable`). Absent/empty means `stable`, so only non-stable skills carry the field. `sciagent list skills` tags non-stable skills `[status]` (stable shown plain); a deprecated *active* skill earns a one-line migrate-off nudge in the `status` Notes section. Soft convention — `validate` does not hard-fail on it.
- **`skills/_attic/` convention** for retired skills: reference-only, off the resolver path (`activate`/`inject` resolve `skills/<name>/`, never `skills/_attic/<name>/`), not walked by `validate`, and listed in a separate "Attic" section of `list skills`. `skills/_attic/README.md` documents revival.
- `docs/skill-lifecycle.md` (the lifecycle, the field, the attic, retire/revive recipes); CONTRIBUTING "Skill lifecycle" subsection; `tests/test_skill_lifecycle.sh` asserting attic exclusion, soft status surfacing, and that `validate` ignores the attic.

### Changed
- All skill walkers — `status.sh` (`_list_skills`, `_list_dependents`), `collisions.sh`, `validate.sh`, `inject.sh` (`--tag`), and the skill-walking tests (`test_skill_scope_lint.sh`, `test_tags_vocabulary.sh`) — now skip every underscore-prefixed dir (`_*`), not just `_TEMPLATE`.
- `architecture-treemap` status normalized `probationary` → `experimental` (documented vocabulary).

### Removed
- `shinymultiome-uio-host` retired to `skills/_attic/` and dropped from `roles/base.yaml` — reference-only; revive per the attic README.

## [3.2.0] - 2026-06-23

Reproducible agentic science: centralize, inject, and enforce the owner's five
standing conventions (figure legibility, results placement, README adjacency,
planning decomposition, reproducibility) so a dropped-in agent follows them
without hand-steering. Three-layer model — centralized rule (CRAFT block),
capability (skill/helper/agent/command), guardrail (validate check + hook).

### Added
- **`SCIAGENT:CRAFT` managed block.** `block.sh` is parameterized by block id; `craft.yaml` is the single source of truth for the five-convention craft text + numeric floors; `lib/sciagent/craft.sh` renders it into AGENTS.md alongside `SCIAGENT:ROLES` on `activate`/`inject`/`eject`, SHA1 drift-detected and idempotently re-rendered. ROLES bytes unchanged.
- **Cross-language figure-style contract.** `lib/figure-style/figure_helpers.{R,py}` (parity-checked): `project_theme`/`set_paper_style`, `save_figure(variant="both")` (print+screen from one plot), `save_overview` (atomic figure+table+caption), `contrast_path`/`overview_path`, `style_series`, `purge_figures`, `write_caption`, `append_master_table`, `round_numeric_cols`, `direction_cue`. Symlinked into analysis repos by `activate` (relative, portable) with a per-project shim + fallback. `analysis_config.yaml:figures` carries dual-context floors (base ≥16 screen / ≥9 print). New `figure-style` skill + `figure-audit` agent.
- **Opt-in `sciagent validate --check` guardrails** (soft-warn default, `--strict` hard-fail): `figure-style`, `results-layout`, `captions`, `provenance`, `freshness` (CRAFT hash staleness + submodule-commit ancestry). `_scratch/`/`$TMPDIR` always exempt; default + activate-internal behavior unchanged.
- **Claude Code hooks** (scaffolded, project-scoped): PreToolUse no-ephemeral + figure-save nudge; Stop caption-sweep. `SCIAGENT_STRICT` toggles warn→block.
- **`sciagent update`** — re-pin the toolkit submodule + re-activate (re-render both blocks, re-link helpers) + report.
- **Planning suite:** gold-standard plan templates (`templates/plan/{00_INDEX,NN_slug}.md.template`); `reasoning-trace` skill; `science-architect` overlay role; orchestration commands `/pipeline-plan`, `/explore-and-plan`, `/add-figure-variant`, `/interpret-storm`.
- New tags `figure`, `provenance`, `planning`.

### Changed
- `scrna-pipeline-conventions` rewritten to the stage-based `03_results/<stage>/{figures,tables}/` layout (flat `plots/`/`checkpoints/` layout removed) and defers to the figure-style contract.
- `captions`, `doc-curator` (C3 path-qualified + C4 committed-script provenance), `handoff` (scripts/artifacts/decisions), `bio-interpreter`/`insight-explorer` (persist before returning) extended.
- Analysis `AGENTS.md.template` references the CRAFT block, figure-style shim, and `docs/_internal/plans/` namespace.

### Fixed
- `sciagent new project --type software` no longer creates an `analysis`-only `research/` namespace.

### Removed
- `docs/guidelines/visualization.md` `base_size=12`/`theme_publication` dead-end superseded by the figure-style contract (file kept as a redirect stub).

## [3.1.0] - 2026-06-02

### Changed
- Project type `software-tool` renamed to `software`; `--type software` now suggests activating `architect` (the software-design lane). No dedicated software role exists.

### Removed
- `roles/software-tool.yaml` removed; the `architect` role owns the software-design lane. No replacement software role — use `sciagent activate architect` in software projects.

### Fixed
- **`inject` resolves the requires-closure.** `sciagent inject <skill>` now walks the skill's transitive `requires:` graph and mounts any missing dependency (recorded with `via="requires:<root>"`), matching what `activate` does for role skills — injecting an orchestrator skill pulls its leaves. A skill already supplied by the active requires-closure is refused (`nothing to inject`) instead of creating a spurious manifest row.
- **`eject` no longer destroys shared dependencies.** Ejecting a closure root now prunes closure deps no longer needed by any other root, and direct eject of an auto-mounted dependency is refused (pointing you to eject the root instead).
- **`status` surfaces requires-inherited skills in all outputs.** Inherited skills now appear in text, `--effective`, `--json`, and `--source` — previously only the default text view. `--effective` now exits 0 on success (was 1).
- **`status` warns on stack drift.** When the manifest pins a role no longer in the catalog, the text view adds a `Notes:` line naming it and pointing at `sciagent deactivate`; `--json` gains a `stale_roles` array.

### Removed
- **`sciagent roster` command removed.** The verb and its module (`lib/sciagent/roster.sh`) have been deleted. Agent-metadata functionality (domain, description) is now served by:
  - `sciagent status` (text mode) — Sub-agents section now shows `domain[0]` and the first 60 chars of description for each activated agent.
  - `sciagent status --json` — each agent object now carries `"domain"` and `"description_brief"` fields.
  - `sciagent list roles` — roles listing now shows skill/agent/command counts per role.
  - `sciagent list role <name>` — new subcommand; prints role description, skills, agents, and commands in detail.
- Deleted 6 redundant roles that bloated the UI and duplicated coverage: `planning`, `annotator`, `scrna-atlas`, `multiome-analysis`, `multiome-grn`, `dc-dictionary`. The clean framework is now 4 roles: `base` (general scRNA foundation), `scatac-regulatory` (chromatin/ATAC overlay), `pathway-signature` (pathway/functional overlay), `architect` (software architecture, standalone).
- Updated `scatac-regulatory` and `pathway-signature` "Do NOT use when" comments to reference the surviving role names (`base`, `scatac-regulatory`) instead of the deleted ones.

### Added
- **`lib/sciagent/frontmatter.sh`** — new leaf module (no sciagent deps) with a pure bash YAML frontmatter parser. Exports `_fm_extract`, `_fm_scalar`, `_fm_list`, `_fm_nested_scalar`, `_fm_description`. Replaces the identical private functions that were baked into `roster.sh`.
- **`sciagent list role <name>`** — new subcommand; prints a detailed RPG-hero view of a role (description, skill list with count, agent list with count, command list with count). Returns exit 1 if the role does not exist.
- **`sciagent list roles`** — enhanced output now shows per-role skill/agent/command counts alongside the description.

## [3.0.0] - 2026-05-31

### Added
- `templates/GEMINI.md.template`: one-line `@AGENTS.md` shim — Gemini provider now has a canonical single-source shim matching the Claude pattern
- `templates/AGENTS.md.template`: appended `## Docs architecture` section with full taxonomy (`docs/stages/`, `docs/reference/`, `docs/_internal/{reasoning,research,plans,reports,handoffs}`), naming rules (YYYYMMDD_HHMMSS snake_case), archive convention (`_archive/YYYYMMDD/`), and boundary rules (no `.md` in `03_results/`)
- `sciagent gitignore [<path>]`: new verb that writes a `SCIAGENT:GITIGNORE` BEGIN/END managed block into `.gitignore`; idempotent; covers `docs/_internal/`, `.claude/`, `.agents/`, `.gemini/`, `.sciagent/`, `.mcp.json`, `.env`, `.env.*`
- `sciagent new project`: now scaffolds full `docs/` tree including `docs/_internal/{reasoning,research,plans,reports,handoffs}` with `.gitkeep` sentinels; writes `docs/README.md`; applies managed `.gitignore` block via `sciagent gitignore`
- `sciagent validate --project-dir <dir>`: six-check docs-layout linter — checks for `docs/`, `docs/_internal/`, hard-fails if `docs/_internal/` exists but is ungitignored, warns on `.md` files in `03_results/`, warns on coexisting archive-style dirs, warns on malformed root handoffs
- `templates/pre-commit.template`: pre-commit hook scaffold that runs `sciagent validate --project-dir` on every commit
- `tests/test_validate_docs_layout.sh`: 5 subtests covering all docs-layout linter checks

### Removed
- `mcp_servers/pal` venv and `deprecated/` harness installers (`install_claude.sh`, `install_codex.sh`, `install_gemini.sh`)

### Changed
- Reorganized `commands/` and `agents/` into role-family subfolders. The resolver in `lib/sciagent/symlinks.sh` (`resolve_canonical`) walks these subtrees recursively so symlinks under `.claude/` and `.agents/` remain flat — consumers see no shape change. Basename uniqueness is enforced by `tests/test_no_duplicate_basenames.sh`.
  - `agents/architect/` — 12 architect-pipeline agents
  - `agents/analysis-base/` — 7 base-role helper agents
  - `commands/architect/` — 14 architect commands; `commands/commit.md` stays top-level
- Renamed role `pathway-signature-agent` → `pathway-signature` to drop the misleading `-agent` suffix (roles live in `roles/`, LLM agents live in `agents/`). Fixed phantom skill references in this role's skills list — replaced 7 deleted skills with the consolidated triad `bulk-rnaseq-gsea`, `bulk-rnaseq-activity-inference`, `bulk-rnaseq-pathway-explorer`.
- Added `muon-multimodal-analysis` (originally added under its old name) to `roles/multiome-grn.yaml` to fill the Python-side 10x ATAC preprocessing / MuData / differential accessibility gap. Subsequently reverted; see below.
- Canonicalised all role-file docstrings (Purpose / Use when / Do NOT use when / optional Composes-with / optional Pipeline) for consistency and brevity.

### Added
- New role `multiome-analysis` (`roles/multiome-analysis.yaml`) for paired RNA+ATAC ingestion, preprocessing, integration, and differential accessibility — without GRN inference. Composes with `multiome-grn` for the downstream regulatory step.

### Changed
- Renamed skill `python-multimodal-10x` → `muon-multimodal-analysis`. The previous name conflated language with purpose; the new name aligns with the audit's namespace convention (defining library = muon/MuData). Directory, frontmatter `name:`, and all cross-references updated; history preserved via `git mv`.

### Reverted
- Dropped `muon-multimodal-analysis` (formerly `python-multimodal-10x`) from `roles/multiome-grn.yaml`. The skill belongs in the new `multiome-analysis` role; `multiome-grn` reverts to its original GRN-inference-only scope. Compose with `multiome-analysis` for preprocessing.

### Deferred (ADR-0002 placeholder)
- Skill-content dissection for the monolithic multimodal skills (`muon-multimodal-analysis`, `seurat-multimodal-analysis`, `multimodal-anndata-mudata`, `cellranger-arc-multiome`) is deferred to ADR-0002. The current role topology bundles these as-is; future work may split each skill into focused sub-skills (e.g. ATAC preprocessing vs. WNN integration vs. peak-gene linkage) without further role changes.

### Removed
- `roles/min.yaml` — unused minimal role.
- `docker/` — CI test scaffolding that was no longer wired into the test suite.

### Migration
After upgrade, projects that have an active sciagent stack must re-activate to refresh symlinks against the new source layout:

```
sciagent deactivate
sciagent activate <base> [overlay]
```

### Changed (prior)
- Removed MCP infrastructure (ToolUniverse, Serena, PAL, Sequential Thinking, Context7), profile switcher, and harness installers (Claude Code, Gemini CLI, Codex CLI). Toolkit now covers roles, agents, and skills only.
- Deleted obsolete docs: `docs/MCP-CONTEXT-MANAGEMENT.md`, `docs/INSTALLATION.md`, `docs/CONFIGURATION.md`, `docs/FAQ.md`, `docs/QUICKSTART.md`, `docs/ARCHITECTURE.md`, `docs/ARCHITECTURE_REVIEW.md`, `docs/CLI_IMPROVEMENT_PLAN.md`, `docs/ISSUES.md`.
- Stripped MCP/harness/profile sections from `README.md`, `CLAUDE.md`, `AGENTS.md`, `CONTRIBUTING.md`, `agents/README.md`, `templates/vendor/CLAUDE.md.template`, `templates/vendor/AGENTS.md.template`, `commands/verify.md`, `docs/workflows/architect/`.

### Fixed (historical — MCP addon tool permissions)
- `manage-addon.sh` now reads `tool_permissions` from addon templates and merges them into `settings.local.json` `permissions.allow`. Previously, addon tools were denied in subagents (Task tool) because subagents can't prompt interactively for permission.
  - `manage-addon.sh`: `update_settings_local()` strips stale `mcp__*` entries and rebuilds from enabled addons
  - `switch-mcp-profile.sh`: Settings generation includes addon tool permissions (survives profile switches)
  - `notebook-tools.addon.json`: Added `tool_permissions` array with all 11 tool names
  - `jupyter.addon.json`: Added empty `tool_permissions` array for future use

## [2.0.0] - 2025-12-16

### Added
- **Role System**: New declarative role-based configuration for agents and skills
  - `roles/base.yaml` - Default bioinformatics analysis role
  - `scripts/activate-role.sh` - Role activator script (symlinks agents/skills to `.claude/`)
  - `skills/` directory for custom skills
- **Template System**: Centralized AI context templates in `templates/vendor/`
  - `CLAUDE.md.template` - Claude Code project instructions
  - `GEMINI.md.template` - Gemini CLI project instructions
  - `AGENTS.md.template` - Universal AI rules for all agents
  - `context.md.template` - Scientific project context
  - `analysis_config.yaml.template` - Analysis parameters
- **Enhanced Profile System** (`switch-mcp-profile.sh`):
  - API key substitution (`${GEMINI_API_KEY}`, `${OPENAI_API_KEY}`, `${CONTEXT7_API_KEY}`)
  - `validate_profile()` function for dependency checking
  - Profile validation before switching
- **setup-ai.sh Enhancements**:
  - Template installation with placeholder substitution
  - Automatic role activation (`activate-role.sh base`)
  - Creates `02_analysis/config/analysis_config.yaml`
- **Architecture Tests** (`docker/test/Dockerfile.architecture-test`):
  - Validates role system, templates, profile switching
  - 8 sub-tests covering all modularization changes
- **Gemini CLI Test** (`docker/test/Dockerfile.gemini-test`)
- **Full test suite** now includes 11 tests (up from 10)

### Changed
- **scbio-docker Integration**: Directory renamed from `01_scripts/` to `01_modules/`
- **Template References**: All templates now reference `01_modules` (not `01_scripts`)
- **research-full.mcp.json**: Added explicit `--include-tools` with 14 curated tools to prevent context overflow

### Fixed
- **MCP Configuration**: Fixed incorrect arguments in `full.mcp.json` and `research-full.mcp.json` templates that caused `switch_mcp full` to fail.
  - Removed incorrect `uv` arguments (`--directory`, `run`) that were passed to the direct binary executable.
  - This ensures `tooluniverse` works correctly in `full` and `research-full` profiles, consistent with `research-lite`.
- **Dev Container npm EACCES**: Fixed npm global install failures in dev containers where nvm is installed to a read-only location (e.g., `/opt/nvm/` owned by root).
  - Added `ensure_npm_writable_prefix()` function to detect read-only npm prefix and auto-configure user-local prefix (`~/.npm-global`)
  - Added `ensure_npm_global_path()` function to ensure previously installed npm binaries are in PATH
  - Updated `install_gemini.sh` and `install_codex.sh` to call these functions before npm install
  - PATH updates are now persisted to `~/.bashrc` automatically
- **Prefix Conflict**: Fixed an issue where `nvm` usage was conflicting with a hardcoded `npm` prefix setting in `common.sh`, causing warnings and potential environment issues.
  - Removed `configure_npm_prefix` call from `ensure_nvm` in `scripts/common.sh`.
  - Removed `prefix` setting from `~/.npmrc`.
  - Removed `~/.npm-global/bin` from `~/.bashrc` to prevent shadowing of `nvm` managed binaries.
- **Claude Code Settings**: Fixed invalid permission glob patterns in `.claude/settings.local.json` generation.
  - Updated `scripts/switch-mcp-profile.sh` to use the correct `Bash(for f in :*)` wildcard syntax instead of invalid `Bash(for f in :*.ext)` patterns.
  - This resolves "Invalid Settings" warnings in `claude doctor` and ensures proper file globbing permissions.
- **MCP Profile Switching**: Fixed `switch-mcp-profile.sh` to correctly generate the `permissions` block in `settings.local.json`, ensuring MCP servers load correctly in Claude Code.

### Changed
- **Dependencies**: Updated `scripts/common.sh` to respect `nvm` environment management and avoid forcing global `npm` prefixes when `nvm` is active.
