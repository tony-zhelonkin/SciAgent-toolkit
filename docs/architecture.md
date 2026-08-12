# sciagent — architecture

This is the canonical design spec for sciagent: a harness-agnostic, per-project context manager that mounts the full catalog of skills, sub-agents, and slash commands into AI assistant sessions via role activation. Roles no longer curate *which* content is mounted (that gating was removed — see §5); a role now decides one thing only: provenance labels — who gets credited for a name in `sciagent status`/the managed block, and which of the two stack slots wins a name collision. Output-style is **not** role-scoped; it is selected at activation time (§6).

---

## 1. Goal

One job: **manage AI-harness context per project**. The user wears different hats (bioinformatician, software engineer, architect). Each hat = a *role*. Activating a role mounts the whole catalog of skills, sub-agents, and slash commands into the project so the active harness session picks them up natively; the role's own `skills:`/`agents:`/`commands:` lists only decide provenance attribution (§5).

Non-goals: installing harnesses, managing MCPs, managing API keys.

**Scope of "per project."** Every *mounting* verb — `activate`, `deactivate`, `update`, `craft`, `lint`, `status` — is strictly project-scoped and writes nothing outside the project directory. That is what makes the reproducibility claim meaningful: the project, plus the toolkit commit it pins, fully determines the mounted context.

One capability sits deliberately outside that scope, and this list used to deny it. `sciagent provision` seeds *user-global* baselines across detected harnesses (see `lib/sciagent/provision.sh`) — the shared `SCIAGENT:CONTEXT` block in each harness's global root, and power-user settings defaults. It exists because a devcontainer wants to run it once at create time. Per **ADR-D4** it is an **opt-in personal bootstrap**, not a tier of normal activation:

- no project verb calls it, and it never calls `activate` or `validate`;
- it is not part of what the toolkit *distributes* — packaging must not select or configure harnesses (**ADR-D3**);
- only the Claude adapter is implemented fully; the other harnesses' settings adapters are logged honestly and skipped, never faked.

Multi-provider posture is therefore a *partial* goal, not a non-goal, and it is worth being exact about how partial: source directories are harness-neutral, `AGENTS.md` + `.agents/skills` are the neutral binding surface, and detection plus global-context paths exist for five harnesses — but full project materialization exists for Claude Code only, with `.agents` mirrors alongside. See `docs/proposals/2026-08-11-offline-distribution/` for the adapter layer this is headed toward (**ADR-D5**).

## 2. Mental model: roles as RPG combo classes

A role is a *label* over `skills + sub-agents + slash-commands` — not a bundle that
gates them, and not a carrier for output-style (§6). A project has at most **two** active roles, stacked in order: `base` (the foundation) and optionally an `overlay` (the specialization). Last-wins on name collisions; the shadowed entry is visible in `sciagent status`.

Two roles, not N, because Claude Code's own three-tier resolution already creates "where did this come from?" debugging pain in practice. Two is enough for "wizard/knight" combos and stays inspectable.

## 3. Directory layout

### In the toolkit (canonical source)

```
sciagent-toolkit/
├── bin/sciagent              # CLI dispatcher (bash)
├── lib/sciagent/             # internal bash modules sourced by bin/sciagent
│   ├── activate.sh
│   ├── deactivate.sh
│   ├── status.sh
│   ├── new.sh
│   ├── block.sh              # AGENTS.md managed-block read/write/verify
│   ├── craft.sh / craft_verb.sh   # SCIAGENT:CRAFT block renderer / `craft` verb
│   ├── roles.sh              # role YAML parser (provenance only — see §5)
│   ├── stack.sh              # stack-walk + catalog-fallback (mounts everything)
│   ├── collisions.sh         # cross-namespace name-collision enumeration
│   ├── validate.sh           # `validate` verb: frontmatter shape + collisions
│   ├── lint.sh               # `lint` verb: opt-in PROJECT guardrail checks
│   ├── ownership.sh          # materialize-and-keep-current discipline for copied
│   │                         #   template bodies (hooks, statusline, helper shims)
│   └── symlinks.sh           # dual-track symlink helpers + the 02_analysis/helpers
│                             #   seam (contract-lib mounts and their shim modules)
├── skills/<name>/SKILL.md    # canonical skills (Anthropic SKILL.md format);
│                             #   allowed top-level keys: name, description, license,
│                             #   allowed-tools, compatibility — no `metadata:` block in use
├── skills/_attic/<name>/     # retired skills — reference-only, off the resolver path (see docs/skill-lifecycle.md)
├── agents/<name>.md          # canonical sub-agents (Claude format)
├── commands/<name>.md        # canonical slash commands (Claude format)
├── system-prompts/<name>.md  # canonical output styles
├── roles/<name>.yaml         # role definitions (provenance labels only)
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
│   └── output-styles/<name>.md ──▶ toolkit/system-prompts/<name>.md
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

- **Markers**: HTML comments — invisible in rendered markdown, distinct namespace (one id per block: `ROLES`, `CRAFT`, `CONTEXT`). The `v1` in the marker is **decorative — nothing parses it**, so there is no version-negotiated upgrade path: a reader expecting `v2` would fail to match a `v1` BEGIN line while still matching the unversioned END line, and report corruption. Treat the format as effectively unversioned until that is fixed.
- **Drift detection**: the header carries `hash=<sha1>` of the block body, checksummed against the body itself — never against what the renderer *would* produce. That is what distinguishes a merely **stale** block (old body, self-consistent hash → re-render proceeds silently) from a **hand-edited** one (body mutated, hash no longer matching → refuse).
  `block_hash_check` returns: `0` match, `1` no block, `2` one marker only (corrupt), `3` drift, `4` no id argument passed (a caller bug, deliberately distinct from `1`).
  1. If both markers found and the hash matches: re-render in place.
  2. If both markers found and the hash mismatches: **`craft` refuses and names `--force`.** `activate`'s ROLES write path does **not** implement this guard — it removes and rewrites unconditionally, so hand-edits inside a ROLES block are lost silently. `status` reports ROLES drift as a note; that is the only warning. **This asymmetry is a known gap.**
  3. If only one marker found: refuse and tell the user to repair the markers by hand. There is no `--force-reset`.
  4. If neither marker found: append a fresh block at EOF, preceded by one blank line. Never silently rewrite existing content, and never change the file's permission bits.
  5. If the file contains **two blocks with the same id**, the state is unrecoverable through the CLI: the body read spans both, so the hash can never match, `craft` reports "drifted (hand-edited)" forever, and `--force` rewrites only the first. Delete the duplicate by hand. Known gap; the message misattributes the cause.
- **Position-shift safe**: search file for markers every run. Never store byte offsets.
- **Single block, not per-role**: the entire stack is rendered into one block. Multi-block patterns invite drift between blocks.
- **CLAUDE.md**: a 1-line `@AGENTS.md` import is the default. Claude's native `@file` import means CLAUDE.md inherits AGENTS.md automatically without duplication or generation.

## 5. CLI surface

Propagation — how a change in the toolkit reaches an already-provisioned project,
and which verb performs which step — is specified in `docs/propagation.md`.

Binary: `bin/sciagent`. Alias: `si`. Install as a relative-path shell alias
(`alias si='./01_modules/SciAgent-toolkit/bin/sciagent'`) once more than one project on
the host vendors its own toolkit copy — a PATH symlink bakes one fixed checkout as the
resolution target and misfires across projects. See README § Install and §13
below.

### Verbs

```
sciagent activate <base> [overlay] [--output-style <name>]   # activate role(s); replaces current stack
sciagent deactivate [<role>]                      # deactivate stack (or single role)
sciagent validate [--quiet]                       # check skill frontmatter shape + name collisions
sciagent lint [--project-dir D] [--check <name>...] [--strict] [--quiet]  # opt-in PROJECT guardrail checks
sciagent status [--json] [--effective] [--source <name>]
sciagent list [roles|skills|agents|commands]      # catalog view; `list role <name>` for one role's detail
sciagent new project|role|skill|agent [args]      # scaffold from templates
sciagent craft [--project-dir D] [--force] [--quiet]   # render/refresh the SCIAGENT:CRAFT block
sciagent gitignore [<path>]                       # add/update the SCIAGENT:GITIGNORE block
sciagent update [--to <ref>] [--no-pin] [--quiet] # re-pin toolkit submodule + re-activate current stack
sciagent provision [--harness <csv|all>] [...]    # seed user-level context + settings per harness
```

There is no `inject`/`eject` verb. That mechanism — adding/removing one entry on
top of the active stack — was removed along with role-based gating: since
`activate` always mounts the whole catalog (§5, and `stack_walk` in
`lib/sciagent/stack.sh`), there is nothing left for `inject` to add that isn't
already mounted, and nothing for `eject` to remove without also hiding it from
every other role.

### Activation semantics

- `sciagent activate base` — solo role (stack: `[base]`). If different stack was active, auto-deactivate first.
- `sciagent activate base reviewer` — two roles, ordered. `reviewer` overlays `base`. Last-wins on collisions.
- Three+ positional args → error: "Maximum stack depth is 2 (base + overlay)."
- Idempotent w.r.t. stack roles: `activate base reviewer` while the same stack is active re-verifies symlinks and re-writes the block if its hash drifted.
- Re-activating with a different stack tears down existing symlinks first (auto-deactivate then activate). One-step UX.

### Validate

`sciagent validate` runs against the toolkit content (see `lib/sciagent/validate.sh`):

1. **Frontmatter shape** — every skill has a `name:` matching its directory and a `description:` of at most `SCIAGENT_DESC_MAX` chars (350 by default). The length is measured on the *folded* value, so a `>-`/`|` block scalar is counted in full rather than by its first physical line — the earlier first-line-only measurement made the cap silently unenforceable for any multi-line description.
2. **Optional skills-ref** — if installed, invoked per skill; surfaced as a warning, silently skipped when absent.
3. **Cross-namespace name collisions** — names appearing under two or more of `skills/`, `agents/`, `commands/`, `roles/`.

Check 1 is **hard-fail** (non-zero exit, error on stderr). Checks 2–3 are **soft-warn** (warning on stderr, "all checks passed" on stdout, exit 0): mounting both is supported and sometimes deliberate. `--quiet` suppresses the success summary.

`validate --check <name>...` is a deprecated compatibility path: it delegates to `sciagent lint` (the opt-in PROJECT guardrail layer split out in Phase 5a) and prints a deprecation note to stderr. The default `validate` path (no `--check`) never runs lint checks and is unaffected.

Two checks that used to run unconditionally here are **no longer on the default path**, because `activate` uses `validate --quiet` as its pre-flight and a consumer-project finding must never be able to block activation:

- **`docs-layout`** moved to `lint` (`sciagent lint --check docs-layout`), where it shares the warn-by-default / `--strict`-escalates hardness of every other project check, and is included in `--check all`. It carries the same absent-subject guard as its siblings: **no `docs/` tree at all is not a finding**, it is an absent subject, so the check returns silently. `docs-layout` audits the structure of an existing `docs/` tree; it does not mandate that one exist. (Only once `docs/` exists does a missing `docs/_internal/` warn — the project has adopted the convention but not finished scaffolding it.) Without that guard the check fired on any directory whatsoever, breaking the "a software or empty project produces zero findings" invariant every other check upholds.
- **`env-hygiene`** stays in `validate` but behind an explicit `--env-hygiene` flag. It inspects the invoking shell's environment rather than anything under `--project-dir`, so `lint`'s "against a consumer project" framing does not fit it.

The "all checks passed" line is emitted only after every check capable of failing has run. It previously printed before the docs-layout check, so `validate` could report success and then exit 1.

There is no requires-graph or tag-vocabulary check anymore — both mechanisms (`metadata.requires:`, tag declarations) were removed along with role-based gating. See `skills/README.md` § Taxonomy.

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

Every skill not attributed to a role in the active stack is still mounted, attributed to `catalog` — this surfaces in every output mode (the text Skills section, `--effective`, `--json`, and `--source`), not just the default text view. Same treatment for agents/commands (Phase 5d).

If the manifest pins a role that has since left the catalog (stack drift), status flags it: the text view adds a `Notes:` line naming the stale role and pointing at `sciagent deactivate`, and `--json` carries a `stale_roles` array.

`--json` emits structured stack + effective table + shadow list + drift state (including `stale_roles`).
`--effective` emits just the merged name list (for piping) and exits 0 on success.
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
```

No `mcp_profile`. No installer/harness fields. No `output_style:` — that field is
not part of the schema anymore (Phase 5d); a role YAML may still carry a stale
`output_style:` key from before the change, but it is never read. Output-style
selection is a runtime choice: `sciagent activate --output-style <name>`, else
`craft.yaml`'s `output_style:` key, else none (see `lib/sciagent/activate.sh`).

`skills:`/`agents:`/`commands:` no longer gate what gets mounted (§5) — every
name in the catalog is mounted regardless of role. These lists only decide
**provenance attribution**: which role a mounted name is credited to in
`sciagent status`/the managed block, and (via last-wins across the stack) which
role wins the attribution on a name collision between base and overlay.

A role MAY declare empty arrays. A minimal role is just `name` + `description`.

## 7. Symlink topology

Every `activate` walks the whole catalog of skills, agents, and commands (§5) — not just what the stack's roles declare — resolves last-wins provenance across the stack, and creates symlinks in BOTH `.claude/<category>/` AND `.agents/<category>/` (skills, agents, commands — *not* output-styles, which is Claude-only and selected separately, see §6).

- Symlinks are **relative** paths back into the toolkit when the toolkit lives inside the project tree (e.g. a submodule under `01_modules/SciAgent-toolkit`): each link is relativized to its own directory (via `realpath -ms --relative-to`), so the committed `.claude/*`/`.agents/*` links stay portable across host/container/machine and travel with the submodule pin. When the toolkit is an **external/global** checkout outside the project tree, the target falls back to an **absolute** path (a relative link would be a fragile `../../../…` chain that breaks when either tree moves). The helper-lib link (`02_analysis/helpers/figure-style`) follows the same relativization rule.
- A toolkit-locality guard refuses `activate` (the only mutating verb — see `MUTATING_VERB` in `bin/sciagent`) when the project ships its own `./01_modules/SciAgent-toolkit` but the resolved `$SCIAGENT_TOOLKIT` is a *different* (external) toolkit — this prevents silently symlinking against, and pinning to, the wrong copy. Override with `--allow-external-toolkit` or `SCIAGENT_ALLOW_EXTERNAL_TOOLKIT=1`. `deactivate` is not guarded (it only removes what its own manifest owns).
- Existing symlinks are removed before new ones are created (clean slate per activate).
- A name appearing in both base and overlay results in one symlink pointing at the overlay's canonical file (last-wins). The shadowed entry is recorded in the managed block, not in symlinks.

## 8. Deactivation

`sciagent deactivate` (no args) — full teardown:
- Remove the ROLES and CRAFT blocks from AGENTS.md (preserve everything outside markers byte-for-byte, including the file's permission bits).
- Remove `.claude/{skills,agents,commands,output-styles}/*` and `.agents/{skills,agents,commands}/*` mounts, sweeping the **union** of two sources: (a) entries that are symlinks resolving inside the *currently active* `$SCIAGENT_TOOLKIT`, and (b) every path recorded in a still-present manifest that is still a symlink. Neither source alone is sufficient: (a) alone strands everything when teardown runs against a *different* toolkit checkout than the one that mounted (the links resolve elsewhere, so nothing matches); (b) alone strands everything when the manifest is lost or stale. (b) is safe to trust unconditionally because manifest entries are never guessed — each was written by `symlink_create_dual`/`symlink_create_helper_lib` at mount time. A file or symlink of the user's own is untouched either way: it is not toolkit-resolving and was never recorded.
- Remove `02_analysis/helpers/{figure-style,interactive-style}` on the same rule.
- Remove `.sciagent/manifest.json`.

**Also removed, each under its own ownership record** (2026-08-11/12; this used to
be documented as a known gap, and was one until those records existed):
`.claude/statusline.sh`, `.claude/hooks/{no_ephemeral,caption_sweep}.sh`, the
`.claude/settings.json` keys the backfill added (per-key, and the user's original
bytes restored verbatim when the file is still exactly what we wrote), the
`02_analysis/helpers/{figure_style.py,figure_style.R,interactive_style.py}` shim
modules, and the `@AGENTS.md` line `activate` prepends to `CLAUDE.md`.

The rule for every one of them is the same, and it is one-shot: reverse it **iff**
its content is still exactly what sciagent wrote, otherwise leave it in place with
a warning on stderr and drop the record — so the decision is made once and a
second `deactivate` is a silent no-op rather than a permanent nag. Directories are
only ever `rmdir`'d, never `rm -r`'d: a non-empty directory means real content
remains. See `lib/sciagent/ownership.sh` for the shared machinery and
`docs/propagation.md` §5 for the discipline in full.

`sciagent deactivate <name>` — partial: remove just that role from the stack. If it's the base, the overlay also goes (overlay without base is meaningless). If it's the overlay, base remains.

## 9. State files

```
.sciagent/
├── manifest.json            # canonical state: stack, created symlinks, schema version
├── claude_settings.state    # how settings.local.json's outputStyle key was set
├── project_settings.state   # which settings.json keys the backfill added (+ .orig)
├── claude_md.state          # created vs. prepended, for the CLAUDE.md import shim
├── statusline.sha1 / .ceded # content snapshot / ceded marker for statusline.sh
├── hook_state/              # one .sha1 (or .ceded) per materialized hook body
└── helper_shim_state/       # one .sha1 (or .ceded) per 02_analysis/helpers shim
```

Everything below `manifest.json` is an **ownership record**: a content snapshot
taken when a materialized body was written or adopted, read exactly once at
teardown (§8). A `.ceded` marker is the opposite — it records a body sciagent
found and did *not* write, so it warns once and never touches the file again.

Schema **v2** (2026-08-11). v1 also carried a `block_hash` field that `activate`
wrote and nothing read; it was dropped. Both readers (`manifest_stack`,
`manifest_symlinks`) are key-targeted, so a v1 manifest still on disk reads
fine and is rewritten to v2 by the next `activate`. A future reader must treat
an unknown `version` as readable, never as fatal.

- AGENTS.md managed block is the **human-readable** state.
- `.sciagent/manifest.json` is the **machine-readable** state — it records the
  active stack (the only machine-readable answer to *which* roles are active) and
  the links `activate` created, which lets `status` detect a mount that has been
  *deleted*. A target-ownership scan cannot do that: it sees what is present, not
  what is missing.
- The manifest is **not** the teardown authority. Teardown resolves ownership from
  each link's target (§8), so the manifest being stale, hand-edited or absent
  changes nothing about what gets removed. There is consequently no block-vs-manifest
  conflict to arbitrate.

## 10. Idempotency

Every operation is safe to re-run:
- `activate base` twice → no-op (verifies + repairs).
- `deactivate` twice → first removes, second is silent no-op (manifest absent).

## 11. New-project scaffolding (`sciagent new`)

Replaces `setup-ai.sh`:

```
sciagent new project [<dir>] [--type analysis|software]
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

Most skills are flat docs-only SKILL.md files an agent reads and interprets. A **packaged skill** is the other tier: a thin SKILL.md interface over a deep, version-locked, tested module that an agent or human outsources execution to via one CLI — deterministic instead of re-derived each run. Recognisable by its shape — `pyproject.toml`, a lockfile, `src/`, `tests/` — and by a tier note in the first line of its body; there is no frontmatter marker. The reference implementation is `skills/mllmcelltype-consensus-annotation/`; the code is the spec, so there is no separate packaged template.

See [docs/packaged-skills.md](packaged-skills.md) for the full contract, distribution stance, and the copy-this checklist.

## 13. Project types

`sciagent new project --type <t>` materializes one of two first-class project shapes. The
template root is `templates/project/{_common,analysis,software}/`: `_new_project()`
renders `_common/` first, then overlays `<type>/`.

### First-class types

| Type | One-liner | Canonical top-level layout | Default role |
|------|-----------|----------------------------|--------------|
| `analysis` | A scientific analysis project. | `00_data/ 01_modules/ 02_analysis/ 03_results/<NN_phase>/ docs/ docs/_internal/` | `base` |
| `software` | A standalone, packageable software library or CLI. | `src/ tests/ docs/ examples/ docs/_internal/` + `tool_config.yaml` / `README.md` | `architect` |

These are the only two `--type` values. Default role is a hint emitted in the "Next steps"
output; it is not auto-activated.

### Umbrella — layout variant of `analysis` (no separate `--type`)

An **umbrella** is an `analysis` project whose integration scope is cross-project:

- Root-level submodules are other *analysis* projects, gitlinked at the repo root — **not**
  inside `01_modules/`. (`01_modules/` holds only software toolkits serving the umbrella.)
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

The stack depth cap of 2 (§5) is a load-bearing invariant. A child `software` project under a
parent's `01_modules/<tool>/` does **not** become a third stack tier on top of the parent's
analysis stack. Instead:

- **Activation root = CWD.** A child toolkit is activated by `cd`-ing into it and running
  `sciagent activate <role>` there. It gets its own `.claude/`, `.agents/`, and AGENTS.md
  managed block rooted at the child directory — its own depth-≤2 stack, independent of the parent.
- Parent analysis and child tool are *different working contexts*, not nested ones. When working
  on the child tool, use the `architect` role (the software-design lane); not "analysis base + tool overlay".
- `sciagent new project --type software` run inside an existing project's `01_modules/`
  emits an informational note pointing at `cd <dir> && sciagent activate architect`.

Cross-root awareness (auto-switching context on `cd`) is a future shell-hook concern, out of
scope for sciagent-core. In the meantime, a CWD-relative `si` alias (§5, README § Install)
gets the practical benefit without a shell hook: it targets whichever root's toolkit copy
you're currently in, safely, in an umbrella with several sibling copies mounted at once.

### `docs/_internal/` namespace by type

`reasoning/` is universal — agents write decision traces there regardless of type. The rest of
the namespace is type-conditional:

| Subdir | Holds | Present in |
|--------|-------|-----------|
| `reasoning/` | decision traces, why-not logs | analysis + software |
| `sessions/` | session handoffs | analysis |
| `scratch/` | throwaway notes | analysis |
| `design/` | design records, API drafts, ADRs for the tool | software |
| `benchmarks/` | benchmark results, profiling logs | software |

The handoff agent targets `sessions/` for analysis and `design/` for software.

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
