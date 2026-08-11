# CLAUDE.md

This file provides guidance to Claude Code when working with the SciAgent-toolkit codebase.

---

## Repository Overview

**SciAgent-toolkit** is a per-project context manager for AI coding harnesses. `sciagent activate <role>` mounts the *entire* catalog of skills, sub-agents, and slash commands into a project's `.claude/` and `.agents/` directories — a role no longer gates which of these get mounted (that filtering was removed). A role now decides two things: **provenance** (which role a mounted name is attributed to in `sciagent status`/the managed block, and which wins on a name collision between the two stack slots) and, optionally, which Claude output-style gets applied. See "Role System" below and `docs/architecture.md` for the full design spec.

**Integration:** Used as a submodule at `01_modules/SciAgent-toolkit/` in analysis projects.

---

## Key Directories

| Path | Purpose |
|------|---------|
| `bin/sciagent` | CLI dispatcher |
| `lib/sciagent/` | Internal bash modules — `activate.sh`, `deactivate.sh`, `status.sh`, `new.sh`, `block.sh`, `craft.sh`, `craft_verb.sh`, `roles.sh`, `symlinks.sh`, `stack.sh` (catalog walk — mounts everything, roles supply provenance only), `frontmatter.sh`, `validate.sh` (frontmatter-shape + name-collision checks), `lint.sh` (opt-in PROJECT guardrail checks, the `lint` verb), `collisions.sh`, `claude_settings.sh`, `harness.sh`, `provision.sh`, `update.sh`, `gitignore.sh` |
| `roles/` | Role definitions (YAML) — provenance labels only, not a mount filter |
| `agents/` | Canonical sub-agent definitions (`.md` files) |
| `skills/` | Canonical skill definitions (`<name>/SKILL.md` format) |
| `commands/` | Canonical slash command definitions (`.md` files) |
| `system-prompts/` | System-prompt / output-style definitions (`.md` files; Claude Code consumer mounts into `.claude/output-styles/`) |
| `templates/` | Project scaffolding templates (`AGENTS.md.template`, `CLAUDE.md.template`; AI docs go in `docs/_internal/`) |
| `tests/` | Bash test suite (`run-all.sh`) |

---

## Common Commands

```bash
# Activate base role (symlinks the whole catalog of agents/skills/commands into
# .claude/ and .agents/; the role only picks provenance labels + output-style)
bin/sciagent activate base

# Check skill frontmatter shape + cross-namespace name collisions
bin/sciagent validate

# Run the opt-in PROJECT guardrail checks against an analysis repo
bin/sciagent lint --project-dir <analysis-repo-dir>

# Run the test suite
bash tests/run-all.sh

# Scaffold a new role
bin/sciagent new role my-role

# Check active stack
bin/sciagent status
```

---

## Role System

**Roles no longer gate what gets mounted.** `sciagent activate <base> [overlay]` always
mounts the entire catalog under `skills/`, `agents/`, `commands/` — a role's
`skills:`/`agents:`/`commands:` lists only decide **provenance attribution** for the
names they list (which role gets credited in `sciagent status`/the managed block, and
which one wins on a name collision between `base` and `overlay`, last-wins). Anything
in the catalog not named by either role's lists still mounts, attributed to `catalog`.
See `docs/architecture.md` §5 and `lib/sciagent/stack.sh` for the mechanism.

Role YAML schema:

```yaml
name: base
description: Default bioinformatics analysis role

skills:
  - skill-name
agents:
  - agent-name
commands:
  - command-name
```

There is no `output_style:` field in the schema anymore — a role YAML may still carry
a stale one from before this changed, but it is never read. Output-style selection is
a runtime choice: `sciagent activate --output-style <name>`, else `craft.yaml`'s
`output_style:` key, else none.

The stack is still capped at 2 (`base` + optional `overlay`) — that cap is a UX/
inspectability constraint on the provenance model, not a mount-capacity limit, since
mounting is unconditional either way. Last-wins on provenance for name collisions
between the two slots.

### Adding a New Agent

1. Create `agents/<name>.md` with YAML frontmatter (`name`, `description`, `model`, `color`)
2. Optionally add it to a `roles/*.yaml` if you want it attributed to that role in `sciagent status`/the managed block — it mounts either way (attributed to `catalog` if you skip this)
3. Test: `bin/sciagent activate base && echo OK`
4. Update `agents/README.md`
5. If the chosen name collides across namespaces (also exists as a skill/command/role), see `CONTRIBUTING.md` for the allowlist procedure.

**Optional frontmatter** (`sciagent status` and `sciagent status --json` read these; absence renders a blank column, never an error):

```yaml
domain:                # 1–3 taxonomy tags, e.g. session-management, documentation
  - session-management
outputs:               # only for agents that write files; omit for read-only agents
  default_path: docs/_internal/sessions/   # canonical default per the naming convention
  kind: session-handoff                     # matches the AGENTS.md routing-table kind
  path_source: AGENTS.md                    # constant; signals the Step-0 protocol applies
```

Any agent that writes a dated artifact must implement the Step-0 output-path protocol in its body — see `docs/architecture.md § Agent output path resolution (Step-0)`.

### Adding a New Role

1. `bin/sciagent new role <name>` — scaffolds `roles/<name>.yaml`
2. Edit the YAML to list the skills/agents/commands you want attributed to this role (provenance only — it does not change what's mounted)
3. Test: `bin/sciagent activate <name>`

### Adding a New Skill

1. `bin/sciagent new skill <name>` — copies `skills/_TEMPLATE/` to `skills/<name>/`
2. Edit `skills/<name>/SKILL.md`
3. It mounts automatically on the next `activate` — no role edit needed. Optionally add it to a `roles/*.yaml` if you want it attributed to that role rather than to `catalog`.
4. If the chosen name collides across namespaces (also exists as an agent/command/role), see `CONTRIBUTING.md` for the allowlist procedure.

---

## After Modifying Agents, Skills, or Roles

Test activation from the toolkit directory:

```bash
bin/sciagent activate base && echo "OK" || echo "FAILED"
ls -la .claude/agents/ .claude/skills/
```

If this toolkit is a submodule, re-activate from the parent:

```bash
./01_modules/SciAgent-toolkit/bin/sciagent activate base
```

---

## Architecture Reference

See [docs/architecture.md](docs/architecture.md) for the full design spec.

For skills that ship version-locked, tested executable code (not just prose), see [docs/packaged-skills.md](docs/packaged-skills.md) — the packaged tier, with `skills/mllmcelltype-consensus-annotation/` as the reference implementation.
