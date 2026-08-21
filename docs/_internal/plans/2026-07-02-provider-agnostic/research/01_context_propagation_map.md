# Context-propagation map (source layer vs materialization layer)

How `sciagent activate` injects context today, and which parts are Claude-hardcoded. Anchors are
file:line in the toolkit.

## What `activate <base> [overlay]` materializes

1. **Provider-harness scaffold, always first** (`activate.sh:97–98`): `claude_settings_ensure_statusline`
   (copy `statusline.sh.template`→`.claude/statusline.sh` if absent, chmod +x) and
   `claude_settings_ensure_project_defaults` (copy or reverse-`jq`-backfill `settings.json`).
2. **Dual-track symlink trees** (`activate.sh:226–242`, `symlinks.sh:260–299`) — `symlink_create_dual`
   writes BOTH `.claude/` and `.agents/` links to the canonical toolkit source (relative if in-tree):
   skills → `{.claude,.agents}/skills/<n>`; agents → `.../agents/<n>.md`; commands → `.../commands/<n>.md`;
   output-style → `.claude/output-styles/<n>.md` (**Claude-only, no `.agents/` mirror**); figure-style
   helper → `02_analysis/helpers/figure-style`.
3. **Output-style activation** (`activate.sh:244–251`): `claude_settings_apply` writes
   `outputStyle:<n>` into `.claude/settings.local.json` (symlink alone doesn't select it).
4. **Managed blocks in `AGENTS.md`** (`activate.sh:264–276`): `SCIAGENT:ROLES` (`block.sh`, SHA1
   drift-guarded) + `SCIAGENT:CRAFT` (`craft.sh`).
5. **Manifest** `.sciagent/manifest.json` (`symlinks.sh:34,116–200`) — records stack, every symlink,
   injected entries, block hash.

## Source (neutral) vs materialization (Claude-hardcoded)

| Provider-neutral SOURCE | Claude-hardcoded MATERIALIZATION |
|---|---|
| `AGENTS.md` source of truth | `.claude/` tree location |
| `roles/*.yaml` concept bundles | `.claude/settings.json` schema (editorMode/effortLevel/hooks) |
| `skills/*/SKILL.md` (agentskills.io) | `.claude/settings.local.json` `outputStyle` |
| `system-prompts/*.md` ("provider-agnostic naming", `AGENTS.md:47`) | `.claude/output-styles/*.md` |
| `craft.yaml`, `tags.yaml`, block markers | `.claude/agents/*.md`, `.claude/commands/*.md` ("Claude-only", `stack.sh:180,189`) |
| `.agents/` dual-track mirror (already written) | `.claude/statusline.sh`, `.claude/hooks/*.sh` |

## Role → {skills, agents, commands, style}

Role = `roles/<name>.yaml` (`roles.sh:24`): scalars `name/description/output_style`, arrays
`skills[]/agents[]/commands[]`. (a) skills → `skills/<n>/` (must have `SKILL.md`), transitive `requires`
closure. (b) agents → `agents/**/<n>.md`. (c) commands → `commands/**/<n>.md`. (d) output_style →
resolved by frontmatter `name:` in `system-prompts/*.md` (`roles.sh:115–160`), selected via
`settings.local.json`. Stack depth ≤2, last-wins (`activate.sh:66–69`).

## Idempotency / ownership

Manifest tracks every artifact; teardown removes only owned symlinks, refuses non-symlinks
(`symlinks.sh:304–329`). Managed blocks SHA1-drift-guarded (`block.sh:4`). settings non-clobber via
reverse `jq -s '.[0]*.[1]'` (existing wins, `claude_settings.sh:216–233`). outputStyle ownership state
in `.sciagent/claude_settings.state` (`created`/`existed-no-style`/`existed-with-style`). Toolkit-locality
guard: mutating verbs refuse an external toolkit against a project shipping its own
(`bin/sciagent:126–166`).

## Pre-existing multi-provider seams

- **`.agents/` dual-track** — every skill/agent/command already mirrored (`symlinks.sh:260–299`,
  `README.md:210–212`). pi/agy consume `.agents/skills/` natively.
- **`.gitignore` block** already ignores `.claude/`, `.agents/`, **`.gemini/`** (`gitignore.sh:13–15`).
- **`GEMINI.md.template`** exists but is **orphan** (no verb materializes it; `new.sh:169–170` only walks
  `_common` + type dirs).
- **`status.sh` detects `.pi/`** (`:300–302,631,647`) — no codex/gemini/opencode detection yet.
- **`delegate-cli` skill** — headless codex/agy invocation knowledge; wired into `science-architect.yaml:43`.

**Bottom line:** source is neutral; only the materialization backend is Claude-hardcoded. Replace it
with a dispatch table; keep everything above it.
