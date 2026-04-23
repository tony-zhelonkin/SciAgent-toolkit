---
parent: ./README.md
view: how-to
---

# Extending the architect role

How to add new reviewers, commands, or adapt the role for other projects.

---

## Add a new reviewer

### Example: add a `security` reviewer (threat modelling lens)

**Step 1** — write the agent file at `agents/security.md`. Use any existing reviewer (e.g., `agents/stat.md`) as a template. The YAML frontmatter must include:

```yaml
---
name: security
description: |
  Security reviewer. Lens: threat modelling, input validation, auth boundaries,
  secrets handling, supply-chain risk, least-privilege.
  <example>
  user: "/review auth-refactor --as security,divergent"
  assistant: "Dispatching security to review threat surface."
  </example>
tools: Read, Grep, Glob, Write, WebSearch
model: sonnet
color: brown
---
```

Body must state:
1. What to evaluate (lens-specific checklist)
2. What this agent does NOT do (scope bounds — don't overlap with `bioinf`/`stat`/etc.)
3. Input contract: "Read `docs/{feature}/map.md` first; do not re-grep"
4. Output format: exactly one file `docs/{feature}/review/security.md`

**Step 2** — add `- security` to `roles/architect.yaml` under `agents:`.

**Step 3** — update `commands/review.md`:
- Add `security` to the roster table
- Add `security` to the `all` set definition in Phase 2

**Step 4** — re-activate the role:
```bash
scripts/activate-role.sh architect --project-dir /path/to/project
```

The script clears old symlinks and creates new ones — the new agent is live.

---

## Add a new meta-command

The meta-layer (`/meta-map`, `/meta-design`, `/meta-plan`) is intentionally minimal — three commands sharing one agent (`meta-architect`). If a project needs more cross-feature orchestration (e.g., `/meta-implement` to drive portfolio-phase rollouts, or `/meta-review` for portfolio-level critique), add it the same way you'd add a per-feature command, with these differences:

**Step 1** — write the command at `commands/meta-<verb>.md`. Use `commands/meta-map.md` as a template. Conventions:
- Phase 0: parse args; default scope = "all features in `docs/_meta/{prior-artifact}.md`'s `features:` frontmatter"
- Phase 1: verify the upstream `_meta/*.md` artifact exists (if your command depends on one)
- Phase 2: check whether the target file exists; offer use / regenerate / iterate
- Phase 3: dispatch `meta-architect` with `subagent_type: "meta-architect"`, naming the target file in the prompt
- Phase 4+: gate with the user; never auto-progress to another meta-command

**Step 2** — if your command needs a different system prompt than the existing `meta-architect`, add a new agent (e.g., `meta-reviewer.md`) following the meta-architect template. Write only inside `docs/_meta/`; never inside per-feature directories.

**Step 3** — add `- meta-<verb>` (and any new agent) to `roles/architect.yaml`.

**Step 4** — re-activate the role.

**Why this is rare**: most cross-feature concerns map cleanly onto the existing three commands. Add a new meta-command only when you have a recurring portfolio-level activity that doesn't fit map / design / plan.

---

## Add a new pipeline command

### Example: add `/audit <feature>` to run a retrospective post-implementation

**Step 1** — decide the command's contract:
- Input: feature slug
- Output: `docs/{feature}/audit.md`
- Dispatches: new `auditor` agent (or reuses existing agents)
- Gate: user review

**Step 2** — write the command at `commands/audit.md`. Use `commands/map.md` as a template. Sections: Phase 0 (parse args), Phase 1 (verify prerequisites), Phase 2 (dispatch), Phase 3 (surface result + next step), Rules.

**Step 3** — if the command dispatches a new agent, write that agent's `.md` file.

**Step 4** — update `roles/architect.yaml`:
- Add `- audit` under `commands:`
- Add the new agent (if any) under `agents:`

**Step 5** — re-activate the role.

---

## Add a new role that reuses the architect role's pieces

### Example: `archdoc-agent` (architecture + documentation writer)

Architect-agent is focused on software changes. You might want a role for producing standalone architecture documentation from an existing codebase — no implementation, no plan.

**Step 1** — create `roles/archdoc-agent.yaml`:

```yaml
name: archdoc-agent
description: Architecture documentation writer — map-only workflow, no implementation
mcp_profile: pal-only

agents:
  - mapper             # reuse
  - bioinf             # reuse (if scientific codebase)
  - ml                 # reuse (if ML codebase)
  - doc-curator        # from base role (SciAgent-toolkit/agents/doc-curator.md)

skills:
  - architecture-first-dev   # reuse

commands:
  - map                # reuse
  - review             # reuse
  - synthesize         # reuse
  # no /design, /plan, /implement
  # no /meta-* — documentation roles work feature-by-feature; meta-layer
  # is only useful when ≥2 features are in flight together
```

**Step 2** — activate the new role:
```bash
scripts/activate-role.sh archdoc-agent --project-dir /path/to/docs-project
```

The activation script is idempotent — you can activate multiple roles in sequence in the same project, but only the last activation's agent/skill/command set will be linked (the symlinks get cleared each time). To compose roles, you'd need a multi-role activation script — not currently implemented.

---

## Adapt the architect role for a non-Claude provider

The role YAML and markdown artifacts are portable. What isn't:

- **Agent files** (`agents/*.md`) use Claude Code's subagent frontmatter (`model: sonnet/opus`, `tools:` list). Other providers don't recognize these fields — they'd be ignored but the body text still works as a system prompt.
- **Commands** are Claude-Code-specific (slash command mechanism). Other providers use different triggering (e.g., Cursor rules, Aider conventions).
- **The activation script** symlinks to `.claude/` — unaware of `.cursor/`, `.aider/`, etc.

### Portability checklist
- Agent bodies port as system prompts — yes
- Artifact conventions (map.md / review/ / design/ / plan/) port — yes, any provider can produce them
- Provider-specific dispatch (parallel subagents, slash commands) — re-implement per provider
- Role YAML as a manifest — portable if you re-write activation for the target provider

---

## Customize a reviewer for domain specialization

### Example: sharpen `bioinf` for single-cell work specifically

The generic `bioinf` reviewer covers bulk + sc + ATAC. If a project is 100% single-cell, you might want a sharper lens.

**Option A — edit in place.** Change `agents/bioinf.md`'s body to focus on sc conventions. The downside: this weakens the reviewer for other projects using the same toolkit.

**Option B — create a variant.** Copy `agents/bioinf.md` → `agents/bioinf-sc.md`, edit for sc focus, and create a new role `roles/architect-sc.yaml` that uses `bioinf-sc` instead of `bioinf`.

**Option C — prompt the reviewer at invocation time.** Pass project context in the feature description when running `/map` — the reviewer reads `docs/{feature}/map.md` and specialises based on what's there. The prompt adapts; the agent file doesn't change.

Preferred: option C for most cases, option B if a team wants a durable specialization, option A if the project is the only consumer of the toolkit.

---

## Migrating pre-harness feature folders

If a project had feature docs before this harness was activated, the layout may not match the harness's expected structure. Pathway-explorer is the canonical example:

```
docs/umap/                                    (pre-harness layout)
├── 01_neighborhood_preservation.md
├── 02_entity_type_mixing.md
├── 03_contrast_dependence.md
├── 04_recommendations.md
└── README.md

docs/normalisation/                           (pre-harness layout, variant)
└── review/
    ├── 1_statistician.md
    ├── 2_computer_scientist.md
    ├── synthesis.md
    └── audit.md
```

Harness convention (what fresh features produce):

```
docs/{feature}/
├── map.md
├── review/
│   ├── bioinf.md
│   ├── wetlab.md
│   └── ...
├── synthesis.md
├── design/
│   ├── README.md
│   ├── 01-architecture.md
│   └── ...
└── plan/
```

### How the harness handles pre-harness layouts

- `/map` writes `docs/{feature}/map.md` alongside whatever was already there. No conflict.
- `/design` reads `synthesis.md` / `audit.md` at the feature root (normalisation-style layout) OR at the review/ subdirectory — both paths work. Pre-existing top-level `NN_*.md` files (umap-style) are read as context but are NOT written to by any harness command.
- `/review` always writes to `docs/{feature}/review/{name}.md`. If the feature doesn't have a `review/` subdirectory yet, it creates one. This means a fresh `/review` run on an umap-style pre-harness folder creates `review/` alongside the legacy top-level files — a split archive.

### Recommended migration (optional)

If you want a unified archive before running `/review` on a pre-harness feature:

```bash
# For umap-style layouts (top-level NN_*.md)
cd docs/{feature}
mkdir -p review
git mv 01_*.md 02_*.md 03_*.md review/   # or whatever the prefix convention was
```

The harness does not automate this — it's a one-time human decision, and the target filenames inside `review/` may also want renaming (e.g., `01_neighborhood_preservation.md` → `review/ml.md` if that maps cleanly). Leave the migration until you actually re-review the feature; legacy files are harmless as read-only history otherwise.

### What not to do

- Don't rename legacy files to match reviewer short names (`bioinf`, `wetlab`, etc.) unless their content actually represents that reviewer's lens — a misnamed legacy file confuses future readers worse than a clearly-legacy filename.
- Don't move files mechanically into `review/` if they're recommendations/syntheses, not reviews. `04_recommendations.md` belongs at the feature root or in a `synthesis.md`, not inside `review/`.

---

## Debug checklist

Symptoms and fixes:

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `/map` says "Command not found" | Role not activated in this project | `scripts/activate-role.sh architect --project-dir .` |
| `/review` dispatches but no files appear | Reviewer agents are writing outside expected paths, or sandboxed | Check reviewer's output message for the actual path; check `.claude/agents/<name>.md` symlink points to a real file |
| `/synthesize` silently does nothing | <2 files in `docs/{feature}/review/` | That's correct behavior — synth is conditional |
| `architect` returns `NEEDS ITERATION` repeatedly | Design docs have persistent consistency gaps | Read `design/review.md` carefully; the Verdict section lists specific items. Fix those, re-run. |
| `/design` stops with "MADR-N conflicts with this feature's decision" | A new MADR landed since the last `/design` run, or the local design legitimately diverges | Either re-run `/meta-design` to revise the MADR, or accept the local override (gets recorded in `03-decisions.md § Meta overrides`) and re-run `/design`. |
| `/meta-design` produces no MADRs | `_meta/map.md`'s "Convergent concerns" section is empty | That's correct: the features in scope don't actually overlap. Skip the meta-layer for this wave. |
| `/meta-plan` flags a collision in every portfolio phase | Two features touch the same file in their phase-01s | Resolve by sequencing one feature's phase-01 to a later portfolio phase, or by splitting the file change into two phases. |
| Per-feature `/design` doesn't read `_meta/design.md` | `_meta/` directory is missing or in the wrong place | `_meta/` must live at `docs/_meta/`, sibling to per-feature dirs. Check the path. |
| Reviewers re-grep the codebase | Reviewer isn't reading `map.md` first | Check the reviewer's system prompt (body of `agents/<name>.md`) — the "Core principles" section should emphasize map.md. If customized, restore the clause. |
| Activation silently skips a command | `commands/<name>.md` doesn't exist | Check the toolkit's `commands/` directory; the activation script warns on missing files. Look for `[WARN] Command not found:` in the activation log. |

---

## Testing changes

There is no unit-test infrastructure for the role system yet. Verification is manual:

1. Activate the role into a test project
2. Check that `.claude/agents/`, `.claude/commands/`, `.claude/skills/` contain the expected symlinks
3. Run each slash command once with a trivial feature slug
4. Check that the expected artifact is produced at the expected path

For comprehensive testing, see `docs/TESTING.md` (SciAgent-toolkit's test doc).
