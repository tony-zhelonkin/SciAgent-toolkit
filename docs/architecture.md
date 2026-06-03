# sciagent — architecture

This is the canonical design spec for sciagent: a harness-agnostic, per-project context manager that injects curated bundles of skills, sub-agents, and slash commands into AI assistant sessions via role activation.

---

## 1. Goal

One job: **manage AI-harness context per project**. The user wears different hats (bioinformatician, software engineer, architect). Each hat = a *role*. Activating a role injects a curated bundle of skills, sub-agents, and slash commands into the project so the active harness session picks them up natively.

Non-goals: installing harnesses, managing MCPs, managing API keys, multi-provider posture, anything global to the user's home directory.

## 2. Mental model: roles as RPG combo classes

A role is a bundle: `skills + sub-agents + slash-commands + output-style`. A project has at most **two** active roles, stacked in order: `base` (the foundation) and optionally an `overlay` (the specialization). Last-wins on name collisions; the shadowed entry is visible in `sciagent status`.

Two roles, not N, because Claude Code's own three-tier resolution already creates "where did this come from?" debugging pain in practice. Two is enough for "wizard/knight" combos and stays inspectable.

## 3. Directory layout

### In the toolkit (canonical source)

```
sciagent-toolkit/
├── bin/sciagent              # CLI dispatcher (bash)
├── lib/sciagent/             # internal bash modules sourced by bin/sciagent
│   ├── activate.sh
│   ├── deactivate.sh
│   ├── inject.sh
│   ├── status.sh
│   ├── new.sh
│   ├── block.sh              # AGENTS.md managed-block read/write/verify
│   ├── roles.sh              # role YAML parser
│   └── symlinks.sh           # dual-track symlink helpers
├── skills/<name>/SKILL.md    # canonical skills (Anthropic SKILL.md format)
├── agents/<name>.md          # canonical sub-agents (Claude format)
├── commands/<name>.md        # canonical slash commands (Claude format)
├── output-styles/<name>.md   # canonical output styles
├── roles/<name>.yaml         # role definitions
└── templates/                # project scaffolding (new project bootstrap)
    ├── AGENTS.md.template    # points at docs/_internal/scientific-context.md
    └── CLAUDE.md.template    # 1-line shim: @AGENTS.md
```

### In an activated project (what `sciagent activate` writes)

```
project/
├── AGENTS.md                              # user-owned + managed block
├── CLAUDE.md                              # 1-line shim: @AGENTS.md (+ overrides)
├── .claude/                               # Claude harness native
│   ├── skills/<name> ──▶ toolkit/skills/<name>
│   ├── agents/<name>.md ──▶ toolkit/agents/<name>.md
│   ├── commands/<name>.md ──▶ toolkit/commands/<name>.md
│   └── output-styles/<name>.md ──▶ toolkit/output-styles/<name>.md
└── .agents/                               # harness-agnostic mirror
    ├── skills/<name> ──▶ toolkit/skills/<name>
    ├── agents/<name>.md ──▶ toolkit/agents/<name>.md
    └── commands/<name>.md ──▶ toolkit/commands/<name>.md
```

`.agents/` is a deliberate convention not yet adopted by Claude Code or Pi natively. Pi reads `.agents/skills/` walking up the tree — that already works. For sub-agents and commands, a Pi extension reading `.agents/agents/` and `.agents/commands/` is the user's job (out of scope for this toolkit).

Project-local only. No user-global (`~/.claude/`, `~/.agents/`) installation. The whole point is per-project context.

## 4. AGENTS.md managed block

> **Note on examples.** The role names (`reviewer`), skill names (`obsidian-vignette`, `simplify`, `cellrank-trajectory`), and command names (`/review`) used in code blocks throughout sections 4–6 are illustrative — substitute real names from `roles/`, `skills/`, and `commands/` when reading these as recipes. The schema and stack semantics are the load-bearing parts; the specific labels are not.

### Marker convention

```markdown
<!-- BEGIN SCIAGENT:ROLES v1 hash=<sha1-of-block-body> -->
# Active roles

Stack (in order, last-wins on name collisions):
1. **base** — `roles/base.yaml` — Default bioinformatics analysis role
2. **reviewer** — `roles/reviewer.yaml` — Code review overlay  *(overlay)*

## Skills (effective)
- `obsidian-vignette` — (base) capture session as Obsidian note
- `simplify` — (reviewer) [shadows base]
- ...

## Sub-agents (effective, Claude-only)
- `code-reviewer` — (reviewer)
- `docs-librarian` — (base)
- ...

## Slash commands (effective, Claude-only)
- `/commit` — (base)
- `/review` — (reviewer) [shadows base]
- ...

<!-- END SCIAGENT:ROLES -->
```

### Robustness rules

- **Markers**: HTML comments — invisible in rendered markdown, distinct namespace (`SCIAGENT:ROLES`), versioned (`v1`).
- **Drift detection**: header includes `hash=<sha1>` of the block body. On every `activate` / `inject` / `deactivate`:
  1. Read file, locate markers.
  2. If both markers found: recompute hash of current body. If hash mismatches stored, user has edited inside the block — print a unified diff and require `--force` (or `sciagent role stash` to preserve edits as a side file).
  3. If only one marker found: abort. Corrupted state. User must fix or pass `--force-reset` to nuke and rewrite.
  4. If neither marker found: append fresh block at EOF, preceded by one blank line. Never silently rewrite existing content.
- **Position-shift safe**: search file for markers every run. Never store byte offsets.
- **Single block, not per-role**: the entire stack is rendered into one block. Multi-block patterns invite drift between blocks.
- **CLAUDE.md**: a 1-line `@AGENTS.md` import is the default. Claude's native `@file` import means CLAUDE.md inherits AGENTS.md automatically without duplication or generation.

## 5. CLI surface

Binary: `bin/sciagent`. Alias: `si`. Both installable via symlink into PATH.

### Verbs

```
sciagent activate <base> [overlay]                # activate role(s); replaces current stack
sciagent deactivate [<role>]                      # deactivate stack (or single role)
sciagent inject [--skill|--agent|--command] <name>   # add one entry on top of current stack
sciagent inject --tag <tag>                       # add all skills carrying <tag>
sciagent eject  [--skill|--agent|--command] <name>   # remove one injected entry
sciagent validate [--quiet]                       # check toolkit integrity
sciagent status [--json] [--effective] [--source <name>]
sciagent roster [--json]                          # list active agents from .claude/agents/
sciagent list [roles|skills|agents|commands]
sciagent new role|skill|agent <name>              # scaffold from templates
```

### Activation semantics

- `sciagent activate base` — solo role (stack: `[base]`). If different stack was active, auto-deactivate first.
- `sciagent activate base reviewer` — two roles, ordered. `reviewer` overlays `base`. Last-wins on collisions.
- Three+ positional args → error: "Maximum stack depth is 2 (base + overlay)."
- Idempotent w.r.t. stack roles: `activate base reviewer` while the same stack is active re-verifies symlinks and re-writes the block if its hash drifted.
- Re-activating with a different stack tears down existing symlinks first (auto-deactivate then activate). One-step UX.
- **Clean-slate w.r.t. injected entries.** Re-activation always tears down the previous manifest, which discards any entries added via `sciagent inject` — even when the stack-roles dimension would otherwise be a no-op. Carrying injected state across `activate` invocations is out of scope (revisit if the depth-2 model changes). The CLI surfaces the loss on STDERR before teardown:

  ```
  sciagent: warning — activate is a clean-slate operation; dropping injected entries:
    - skill s_extra
    - agent ag_helper
    to preserve, run 'sciagent deactivate' first and re-inject after.
  ```

### Inject semantics

- `sciagent inject <name>` auto-detects the kind by scanning the three canonical directories (`skills/<name>/SKILL.md`, `agents/<name>.md`, `commands/<name>.md`).
- Unambiguous match → mounted into `_injected` with manifest `kind` set accordingly. Symlinks land in the kind-appropriate directories (skills in `.claude/skills/` + `.agents/skills/`; agents in `.claude/agents/`; commands in `.claude/commands/`).
- Ambiguous match (same name exists in 2+ canonical directories) → hard-fail; escape hatch is the explicit flag `--skill <name>`, `--agent <name>`, or `--command <name>`.
- When an explicit-flag inject succeeds but a companion entry of another kind also exists under that name, a stderr note advertises it. Companion entries are never auto-mounted.
- Unknown name (no canonical file under any kind) → hard-fail.
- `sciagent inject --tag <tag>` is the bulk form: adds all skills whose `SKILL.md` carries `<tag>`. Tag form is skill-only.
- If the active stack is `[base]`, inject creates an implicit anonymous overlay `_injected`. Stack becomes `[base, _injected]`. Inject never creates a third tier.
- `sciagent deactivate _injected` removes the injected-only overlay. `sciagent deactivate <role>` removes a named overlay AND any entries injected into it.

### Eject semantics

- `sciagent eject <name>` is symmetric to inject and removes one entry from the active overlay.
- Kind discriminator is the manifest `kind` field (not a re-scan of canonical directories), so the entry that gets removed is the one that was actually injected.
- Auto-detection on bare name: if `<name>` was injected under exactly one kind, eject removes it. If `<name>` was injected under 2+ kinds, eject hard-fails and requires `--skill` / `--agent` / `--command`.
- Eject refuses to remove an entry that came in via a stack-mounted role (those are owned by `activate`/`deactivate`).
- Unknown / not-injected name → hard-fail; manifest and symlinks are untouched.

### Validate

`sciagent validate` runs four checks against the toolkit content:

1. **Requires-graph** — every `requires:` reference resolves; the graph is acyclic.
2. **Tag vocabulary** — every tag used by a skill is declared in the tag registry.
3. **Optional skills-ref** — every `optional_skills:` reference resolves.
4. **Cross-namespace name collisions** — names appearing under more than one of `skills/`, `agents/`, `commands/`, `roles/`.

The first three are **hard-fail** (non-zero exit, error on stderr). The fourth is a **soft-warn** (warning on stderr, "all checks passed" on stdout, exit 0): mounting both is supported and sometimes deliberate. `--quiet` suppresses all output and emits only the exit code.

`validate` is allowlist-blind by design — it reports every collision. The allowlist-aware view lives in `sciagent status`'s Notes section, which annotates intentional family overlaps.

### Status output

Default (TTY, plaintext, no color magic):

```
Stack:
  1. base       roles/base.yaml      Default bioinformatics role
  2. reviewer   roles/reviewer.yaml  Code review overlay     [overlay]

Skills (12 effective):
  obsidian-vignette        base
  simplify                 reviewer    (shadows base)
  ...

Sub-agents (4 effective, Claude-only):
  code-reviewer            reviewer
  docs-librarian           base
  ...

Slash commands (5 effective, Claude-only):
  /commit                  base
  /review                  reviewer    (shadows base)
  ...

Managed block: AGENTS.md  lines 14–58   hash OK
Symlinks:      .claude/* OK   .agents/* OK
Harness:       Claude Code detected (.claude/ present)   Pi: not detected (.pi/ absent)
```

`--json` emits structured stack + effective table + shadow list + drift state.
`--effective` emits just the merged name list (for piping).
`--source <name>` resolves a single name to its providing role.

## 6. Role YAML schema

```yaml
# roles/base.yaml
name: base
description: Default bioinformatics analysis role
skills:
  - obsidian-vignette
  - cellrank-trajectory
agents:                    # Claude-specific; ignored by Pi until extension exists
  - code-reviewer
  - docs-librarian
commands:                  # Claude-specific
  - commit
  - review
output_style: architect-mentor   # optional, Claude-specific
```

No `mcp_profile`. No installer/harness fields. Just content bundles.

A role MAY declare empty arrays. A minimal role is just `name` + `description`.

## 7. Symlink topology

Every `activate` walks the union of declared `skills + agents + commands + output_style` across the stack (last-wins) and creates symlinks in BOTH `.claude/<category>/` AND `.agents/<category>/` (skills, agents, commands — *not* output-styles, which is Claude-only).

- Symlinks are absolute paths back into the toolkit. This survives toolkit being a submodule.
- Existing symlinks are removed before new ones are created (clean slate per activate).
- A name appearing in both base and overlay results in one symlink pointing at the overlay's canonical file (last-wins). The shadowed entry is recorded in the managed block, not in symlinks.

## 8. Deactivation

`sciagent deactivate` (no args) — full teardown:
- Remove the managed block from AGENTS.md (preserve everything outside markers byte-for-byte).
- Remove `.claude/skills/*`, `.claude/agents/*`, `.claude/commands/*`, `.claude/output-styles/*` symlinks created by sciagent only (track via a sidecar `.sciagent/manifest.json` so we never delete symlinks we didn't create).
- Remove `.agents/skills/*`, `.agents/agents/*`, `.agents/commands/*` analogously.
- Remove `.sciagent/manifest.json`.

`sciagent deactivate <name>` — partial: remove just that role from the stack. If it's the base, the overlay also goes (overlay without base is meaningless). If it's the overlay, base remains.

## 9. State files

```
.sciagent/
└── manifest.json       # canonical state: stack, created symlinks, block-hash, schema version
```

- AGENTS.md managed block is the **human-readable** state.
- `.sciagent/manifest.json` is the **machine-readable** state — used for safe teardown (we only remove symlinks we own) and drift detection.
- On conflict between block and manifest: print warning, manifest wins for teardown purposes, block is rewritten.

### Injected-entry row schema

Each injected entry is one row with shape `{overlay, skill, via, kind}` where `kind ∈ {skill, agent, command}`. Wire format is pipe-delimited (`overlay|skill|via|kind`). A row missing `kind` is read as `skill` (forward-compat with pre-extension manifests). The row's name field is `skill` for all kinds — a stable schema quirk retained for backward compatibility; `kind` is the authoritative discriminator.

## 10. Idempotency

Every operation is safe to re-run:
- `activate base` twice → no-op (verifies + repairs).
- `inject X` twice → first creates symlink, second is no-op.
- `deactivate` twice → first removes, second is silent no-op (manifest absent).

## 11. New-project scaffolding (`sciagent new`)

Replaces `setup-ai.sh`:

```
sciagent new project [<dir>] [--type analysis|software-tool]
                               # bootstrap a typed project tree, do NOT activate
sciagent new role <name>       # scaffold roles/<name>.yaml from template
sciagent new skill <name>      # copy skills/_TEMPLATE/ to skills/<name>/
sciagent new agent <name>      # scaffold agents/<name>.md from template
```

`sciagent new project` materializes the directory tree for `--type` (default `analysis`),
renders `templates/project/_common/` first then `templates/project/<type>/` over it, seeds a
shared `.gitignore`, and drops `.gitkeep` into empty dirs. It does NOT activate any role; the
user picks one with `sciagent activate` (see the type→default-role hint in §13).

## 12. Packaged skills

Most skills are flat docs-only SKILL.md files an agent reads and interprets. A **packaged skill** is the other tier: a thin SKILL.md interface over a deep, version-locked, tested module that an agent or human outsources execution to via one CLI — deterministic instead of re-derived each run. Marked with `packaged: true` in frontmatter (orthogonal to `scope:` / `tier:`, which describe the doc, not the packaging). The reference implementation is `skills/mllmcelltype-consensus-annotation/`; the code is the spec, so there is no separate packaged template.

See [docs/packaged-skills.md](packaged-skills.md) for the full contract, distribution stance, and the copy-this checklist.

## 13. Project types

`sciagent new project --type <t>` materializes one of two first-class project shapes. The
template root is `templates/project/{_common,analysis,software-tool}/`: `_new_project()`
renders `_common/` first, then overlays `<type>/`.

### First-class types

| Type | One-liner | Canonical top-level layout | Default role |
|------|-----------|----------------------------|--------------|
| `analysis` | A scientific analysis project. | `00_data/ 01_modules/ 02_analysis/ 03_results/<NN_phase>/ docs/ docs/_internal/` | `base` |
| `software-tool` | A standalone, packageable software library or CLI. | `src/ tests/ docs/ examples/ docs/_internal/` + `tool_config.yaml` / `README.md` | `software-tool` |

These are the only two `--type` values. Default role is a hint emitted in the "Next steps"
output; it is not auto-activated.

### Umbrella — layout variant of `analysis` (no separate `--type`)

An **umbrella** is an `analysis` project whose integration scope is cross-project:

- Root-level submodules are other *analysis* projects, gitlinked at the repo root — **not**
  inside `01_modules/`. (`01_modules/` holds only software-tool toolkits serving the umbrella.)
- An `integration/` directory plays the `02_analysis/` role; its scripts consume only each
  child's `03_results/` published surface — never child `02_analysis/` or intermediate state.

There is no `--type umbrella`. Initialize with `--type analysis`, then add child submodules
and create `integration/` by hand. Umbrella projects use `base` as their role (same as analysis).
DC-nexus is the canonical instance.

### Deferred / catalogued-only types

`pipeline` (Nextflow/Snakemake), `paper` (manuscript + figures), and `data-package` (curated
dataset + loader) are real categories with no present-day instance demanding a scaffold. The
vocabulary is reserved here; no `templates/project/<type>/` directory exists for them. They are
the extension points for a future `--type` value: add `templates/project/<type>/` plus a
`_dirs_for_type` branch.

### Nested-toolkit activation (depth-2 cap, independent per-root)

The stack depth cap of 2 (§5) is a load-bearing invariant. A child `software-tool` under a
parent's `01_modules/<tool>/` does **not** become a third stack tier on top of the parent's
analysis stack. Instead:

- **Activation root = CWD.** A child toolkit is activated by `cd`-ing into it and running
  `sciagent activate <role>` there. It gets its own `.claude/`, `.agents/`, and AGENTS.md
  managed block rooted at the child directory — its own depth-≤2 stack, independent of the parent.
- Parent analysis and child tool are *different working contexts*, not nested ones. When working
  on the child tool you want the `software-tool` role, not "analysis base + tool overlay".
- `sciagent new project --type software-tool` run inside an existing project's `01_modules/`
  emits an informational note pointing at `cd <dir> && sciagent activate software-tool`.

Cross-root awareness (auto-switching context on `cd`) is a future shell-hook concern, out of
scope for sciagent-core.

### `docs/_internal/` namespace by type

`reasoning/` is universal — agents write decision traces there regardless of type. The rest of
the namespace is type-conditional:

| Subdir | Holds | Present in |
|--------|-------|-----------|
| `reasoning/` | decision traces, why-not logs | analysis + software-tool |
| `sessions/` | session handoffs | analysis |
| `scratch/` | throwaway notes | analysis |
| `design/` | design records, API drafts, ADRs for the tool | software-tool |
| `benchmarks/` | benchmark results, profiling logs | software-tool |

The handoff agent targets `sessions/` for analysis and `design/` for software-tool.

Naming conventions for files within this namespace live in the scaffold itself
(`docs/_internal/README.md`), not in `CLAUDE.md` or `AGENTS.md`. Every new project ships
that reference; agents read it at runtime to name dated artifacts.

### Agent output path resolution (Step-0)

Any agent that writes a dated artifact resolves its output directory at runtime rather than
hardcoding it. The protocol is three steps, executed by the LLM running the agent:

```
Step 0: Read AGENTS.md (project root). Find the `## Documentation namespace` section.
        Find the routing-table entry for your output kind. Use that directory.

Fallback: if AGENTS.md has no Documentation namespace section, or no routing entry
          matches your kind, use `outputs.default_path` from your own frontmatter.
          Proceed silently — do not stop or ask.

Never write to the project root. Never hardcode a path that includes a project name,
user name, or absolute filesystem location.
```

This is prose in the agent body, not a library call. It keeps agent definitions
project-agnostic (no baked-in paths) while letting each project declare its own routing in
the always-present AGENTS.md. The `outputs.default_path` frontmatter field is the canonical
default per the convention, so an agent works correctly even before a project configures its
managed block.

---

**Status:** Spec version 1. Implementation tested via `tests/run-all.sh`.
