# Architecture — philosophy of context propagation, made provider-agnostic

This doc has two halves. **Part A** states the toolkit's existing philosophy (unchanged — we are
extending it, not replacing it). **Part B** is the provider-agnostic model that realizes that
philosophy across five coding CLIs.

---

## Part A — The philosophy (what SciAgent context IS and how it propagates)

### A.1 The load-bearing thesis: the scaffold IS the interface

> *"In an agent-operated repository, the layout, the docs, and the context files are not documentation
> **of** the work — they are the **interface** through which the work is done."*
> — `2026-06-04_mindpalace-craft-as-architecture.md:18`

> *"The scaffold IS the interface. A convention the owner re-types is a convention that is not in the
> scaffold, not single-sourced, or not enforced."*
> — `2026-06-23-repoducible-agentic-science/00_ARCHITECTURE.md:16`

An agent behaves well not because it is told to in the moment, but because **the repository it woke up
in already encodes that disposition.** Reproducibility, figure style, provenance discipline,
planning-before-coding — these are properties of the *environment*, injected once and enforced, not
prompt hygiene the human re-types each session.

### A.2 The three-layer enforcement model (the guardrail mechanism)

Every convention is expressed three times, once per layer (`00_INDEX.md:95`,
`00_ARCHITECTURE.md:38`):

| Layer | Vehicle | Example |
|-------|---------|---------|
| **(a) Stated once** | An always-on **managed block** in `AGENTS.md` (`SCIAGENT:CRAFT`, `SCIAGENT:ROLES`) | "Every figure is saved via `save_figure`" |
| **(b) Made doable** | A **skill / helper lib / command** | `lib/figure-style/figure_helpers.{py,R}`, `skills/figure-style/` |
| **(c) Made unskippable** | A **`validate --check` + hook** | `validate --check figure-style`; PreToolUse/Stop hooks |

> *"hooks/validate > goodwill. Every convention is stated once (CRAFT), made doable (skill/helper/command),
> made unskippable (validate check + hook)."*

Crucially: **ethos lives in an always-on managed block, not a skill** — *"a skill the model has to
choose to consult is a skill it will skip under pressure"* (`2026-06-04…:68`). This is the reason the
provider-agnostic design must guarantee layer (a) reaches **every** harness's context file, and layer
(c) reaches **every** harness's hook surface — not just Claude's.

### A.3 The propagation primitives that already exist (and work)

- **Managed blocks** — single-sourced in `craft.yaml` / role YAML, rendered by `activate` into
  `AGENTS.md` between HTML-comment markers with a SHA1 drift guard (`block.sh:4`, `craft.sh`). *No one
  hand-edits them; `activate` re-renders; they cannot drift.* This is the proven mechanism for
  "centralized-and-propagating text."
- **Symlink trees** — skills/agents/commands symlinked from the canonical toolkit into the project's
  harness dirs (`symlinks.sh:260`). Relative, pin-respecting links when the toolkit is in-tree.
- **Shared code** — `lib/figure-style/` symlinked into `02_analysis/helpers/` so a style fix
  propagates on re-activate.
- **Manifest + ownership state** — `.sciagent/manifest.json`, `.sciagent/claude_settings.state` —
  every injected artifact is tracked so teardown removes only what the toolkit owns, and merges never
  clobber user values (reverse `jq` merge, existing-wins).

### A.4 What "guardrail to my coding style" concretely means

These files *are* the coding style; provider-agnosticism means all of them must reach every harness:

- `craft.yaml` → the six standing CRAFT conventions (Figures, Results, README, Planning,
  Reproducibility, Memory & traceability).
- `docs/guidelines/code_style.md` (526 lines) → the depth the `code-reviewer`/`doc-curator` agents cite.
- `skills/{figure-style,reasoning-trace,scrna-pipeline-conventions,architecture-first-dev}/`.
- `sciagent validate --check {figure-style,results-layout,captions,provenance,freshness}`.
- The Claude Code hooks (no-ephemeral PreToolUse; caption-sweep Stop).

---

## Part B — The provider-agnostic model

### B.1 The key structural fact: source is neutral, materialization is not

From `research/01_context_propagation_map.md`:

| Already provider-neutral (the SOURCE layer) | Claude-Code-hardcoded (the MATERIALIZATION layer) |
|---|---|
| `AGENTS.md` as source of truth | `.claude/` tree location |
| `roles/*.yaml` (concept bundles) | `.claude/settings.json` schema (editorMode, effortLevel, hooks) |
| `skills/*/SKILL.md` (agentskills.io standard) | `.claude/settings.local.json` `outputStyle` key |
| `system-prompts/*.md` (deliberately provider-agnostic name) | `.claude/output-styles/<name>.md` |
| `craft.yaml`, `tags.yaml`, managed-block markers | `.claude/agents/*.md`, `.claude/commands/*.md` (labeled "Claude-only") |
| `.agents/` dual-track mirror (already written today) | `.claude/statusline.sh`, `.claude/hooks/*.sh` |

**Design consequence:** we do **not** rewrite the source layer. We replace the single hardcoded
materialization backend with a **provider dispatch table**. Everything above the table stays.

### B.2 The portable substrate: `AGENTS.md` + `.agents/skills/`

From `research/02_cli_capability_matrix.md` — AGENTS.md is read **natively by 4 of 5** target CLIs:

| CLI | Reads `AGENTS.md`? | Skills convention | Shim needed |
|-----|-------------------|-------------------|-------------|
| codex | ✅ canonical | `~/.codex/skills/` (agentskills.io) | none |
| pi | ✅ (`AGENTS.md` **or** `CLAUDE.md`) | `.agents/skills/`, reads `~/.claude/skills` too | none |
| agy | ✅ (`AGENTS.md` **and** `GEMINI.md`) | `.agents/skills/<name>/SKILL.md` | none |
| opencode | ✅ | `.opencode/skills/`, `instructions[]` | none |
| claude code | ❌ (reads `CLAUDE.md`) | `.claude/skills/` | `CLAUDE.md` = `@AGENTS.md` (already done) |

**Therefore the canonical injection target flips from `.claude/` to `AGENTS.md`.** Claude Code becomes
the *special case* (needs a one-line import shim) instead of the default. This is the single most
important architectural move in the plan, and it costs almost nothing — the umbrella already ships
`CLAUDE.md` = `@AGENTS.md`.

### B.3 The dispatch table (shape B, finally built)

`docs/proposals/ai-research/11` recommended "harness-agnostic core (bash+Python) + thin per-harness
adapters (~150 LOC each)." Concretely, replace the direct calls in `activate.sh` /
`claude_settings.sh` with a table keyed by harness:

```
lib/sciagent/harness/
  common.sh        # shared: AGENTS.md block render, .agents/skills symlink, CRAFT
  claude.sh        # .claude/ tree, settings.json, output-style, statusline, hooks   (extract existing)
  codex.sh         # ~/.codex/config.toml keys, AGENTS.md (native), skills dir
  agy.sh           # ~/.gemini/antigravity-cli/, AGENTS.md/GEMINI.md (native), .agents/skills
  opencode.sh      # opencode.json, agents/*.md (native subagents!), AGENTS.md
  pi.sh            # .pi/settings.json, AGENTS.md (native), .agents/skills, → pi extension (tier 3)
```

Each adapter implements the same contract (`harness_ensure_context`, `harness_ensure_skills`,
`harness_ensure_settings`, `harness_ensure_hooks`, `harness_teardown`) and reports what it can/can't
do. `activate`/`provision` loop over the **detected or requested** harnesses. Detection extends the
existing `status.sh` harness probe (which already detects `.pi/`).

### B.4 Two scopes: project-level vs user-level

Today the toolkit only writes **project-level** settings, and only on `activate`. The provider-agnostic
model needs both scopes, because a fresh container install (Bug 2) and cross-provider defaults live at
**user level**:

| Scope | What lands | Trigger | New in this plan |
|-------|-----------|---------|------------------|
| **Project** | role skills/agents/commands, output-style, `.claude/settings.json`, CRAFT/ROLES blocks in `AGENTS.md` | `sciagent activate` | dispatch table (tier 2) |
| **User** | power-user defaults (`~/.claude/settings.json`), global `AGENTS.md`-class context, global skills | `sciagent provision` (new verb), run at devcontainer `postCreateCommand` | **entirely new (tier 1)** |

`provision` is the new verb. It is idempotent, non-clobbering (same reverse-`jq` merge as project
settings), and per-harness via the same dispatch table. It is what a devcontainer runs once at create
time so every harness comes up wearing the user's defaults.

### B.5 Idempotency & ownership carry over unchanged

The manifest / state-tag / reverse-merge / SHA1-drift machinery (A.3) is provider-neutral already.
Each adapter writes its own manifest section and its own ownership state. The invariants hold across
all harnesses: **toolkit owns only what it wrote; user values always win on merge; teardown never
deletes a non-symlink or a hand-edited block.**

### B.6 Where pi goes beyond adapters (tier 3)

For four of the five harnesses, an adapter (materialize context + skills + settings) is the whole
story. **pi is the exception**: it exposes a userland extension API (`pi.registerTool`,
`pi.registerCommand`, `pi.on(event)`) and the installed package ships a *working sub-agent example*
(`examples/extensions/subagent/`). So for pi we can do what no adapter can: route a SciAgent **role**
into **parallel or chained isolated sub-agents**, auto-selected from the role menu or named explicitly.
That is Tier 3, and it is the part of this plan that turns "inject context" into "orchestrate the way I
would." See `40_pi_extension_tier3.md`.
