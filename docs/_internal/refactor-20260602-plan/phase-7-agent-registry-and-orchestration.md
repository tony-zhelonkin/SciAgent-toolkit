# Phase 7: Agent Output Contracts and Orchestration Discoverability

## Summary

This phase solves two connected problems. First: agents write to wrong or inconsistent
locations because they have either hardcoded stale paths (handoff → project root,
bio-interpreter → `research_notes.md` at root) or no declared path at all. Second: an
orchestrator arriving cold — whether that is the main Claude Code session at the start of
the day, or a future workflow agent — cannot enumerate the active agent roster and map a
task to the right agent without being explicitly told.

The design stays deliberately light. This is not a schema registry, not a compiled
manifest, not a dependency graph. The repeatable structure we want is the kind that emerges
from a simple convention reliably applied, not from bolts tightened until the architecture
becomes rigid. Two conventions accomplish both goals:

1. **The Step-0 protocol**: every agent that writes a dated artifact reads its output path
   from the project's AGENTS.md routing table at runtime rather than hardcoding it.
   AGENTS.md is already the always-present runtime communication channel between project
   and agent — this makes it carry one more piece of information it is perfectly suited
   for. The agent doesn't know anything about the project by name or structure; it only
   asks "where does the routing table say I should write this kind of file?" and falls
   back to the canonical default if no table is present.

2. **Lightweight frontmatter metadata**: two optional fields — `domain:` (taxonomy tags)
   and `outputs.default_path` (the canonical path per the convention) — make agent
   definitions machine-scannable without imposing a rigid schema. An orchestrator can
   `grep` frontmatter from `.claude/agents/*.md` and produce a roster without any
   additional tooling. A `sciagent roster` command materializes this into a readable
   table on demand.

Everything else this phase does is repair work: the naming conventions agreed in the
phase-3 discussion were never written down; the handoff agent body has never been
updated to follow them; the doc-curator body was written before the phase-3/4 checks
existed; bio-interpreter outputs need a path update. Phase 7 executes those repairs and
closes the gap between what the plan docs say and what the agent files actually do.

---

## Current State

### Handoff agent (`agents/analysis-base/handoff.md`)

- **Output path**: project root, filename `handoff_YYYYMMDD_HHMMSS.md`
- **Archive logic**: moves old files to `.handoff_archive/` — exactly the root clutter
  the docs namespace is designed to eliminate
- **Context files it reads**: `plan.md` and `tasks.md` at root — neither of which is
  part of the phase-3 scaffold anymore; `plan.md` does not exist and `tasks.md` is
  user-maintained
- **No Step-0**: hardcoded paths, no convention lookup

Phase-3 (section 3.4) specifies the target: `docs/_internal/ai-generated/sessions/YYYY-MM-DD_<slug>.md`,
no archive dir, dated files accumulate as continuity record. That spec was written;
the agent body was not updated.

### bio-interpreter (`agents/analysis-base/bio-interpreter.md`)

- **Output path**: `research_notes.md` — no directory, implied project root
- This contradicts the phase-3 routing table which puts lit-search output in
  `docs/_internal/ai-generated/research/YYYY-MM-DD_NN_<topic>.md`
- Also contradicts the description's own handoff note: "Mechanism explanation in
  `research_notes.md`" needs to be updated

### doc-curator (`agents/analysis-base/doc-curator.md`)

- Body was written before phase-3/4 conventions existed
- Missing three checks specified in phase-3 section 3.6 and phase-4 section 4.9:
  1. Flag AGENTS.md/CLAUDE.md > ~150 lines and context.md that is not a pointer
  2. Scan public-facing READMEs for `_internal/` path references (one-way rule violations)
  3. List files in `03_results/**` with no matching `## <filename>` heading in the
     sibling README.md (uncaptioned artifacts)
- The description is a single 600-character prose string with embedded `\n` escapes —
  hard to read in source, hard for Claude Code to parse for dispatch decisions

### Frontmatter across all agents

Current fields: `name`, `description`, `model`, `color`. No `domain:` tags, no `outputs:`
declarations. An orchestrator or human arriving cold at `.claude/agents/` must read
every full description to understand the taxonomy, and cannot know where any agent writes
without reading its full body.

### Naming conventions

The conventions agreed during phase-3 discussion were never written into any file:

| Artifact type | Agreed naming |
|---|---|
| Research / lit search | `YYYY-MM-DD_NN_<topic>.md` |
| Codebase / data exploration | `YYYY-MM-DD_NN_<what>.md` |
| Session handoff | `YYYY-MM-DD_<slug>.md` (NN only on same-day collision: `YYYY-MM-DD_NN_<slug>.md`) |
| Decision records | `<slug>.md` (dateless; date lives in git log) |
| Phase plans | `YYYY-MM-DD__<phase-name>/` (directory) |

The `NN` sequence number for within-day artifacts: `printf "%02d" $(( $(ls docs/_internal/ai-generated/<subdir>/ 2>/dev/null | grep "^$DATE" | wc -l) + 1 ))` — no wall-clock dependency, no risk of collision from parallel agents.

---

## Changes

### 7.1 Document naming conventions as an authoritative reference

**Files:** `templates/docs/_internal/README.md` (new section); `docs/architecture.md` (reference)

**What:** Write the naming convention table above into `templates/docs/_internal/README.md`
so that every new project ships with it and every agent can find it by reading `docs/_internal/README.md`.
The table should also specify:
- The sequence-number computation for within-day `NN` (the bash one-liner above)
- The fallback: if the directory doesn't exist yet, `NN = 01`
- For decisions: explicitly dateless — the slug must be stable across time, the git log
  carries the date

This section belongs in `_internal/README.md` because that file is the namespace contract.
It should also be referenced (not duplicated) from `AGENTS.md.template`'s Documentation
Namespace section with a single pointer line: "Naming conventions: see `docs/_internal/README.md`."

**Why:** Agent bodies that need to create a dated file currently have no authoritative
reference; they invent names. Writing the convention into the always-present README of the
`_internal/` directory gives agents one place to read that is stable, project-local, and
not inside a harness config file.

**How:**
1. Add a `## Naming conventions` section to `templates/docs/_internal/README.md` with the
   table above and the NN computation.
2. Add one pointer line to `templates/AGENTS.md.template` in the Documentation Namespace
   section: "File naming: `docs/_internal/README.md`."
3. Note in `docs/architecture.md` that naming conventions live in the scaffold, not in
   CLAUDE.md or AGENTS.md.

---

### 7.2 Define the Step-0 convention-reading protocol

**Files:** Agent bodies for `handoff.md`, `bio-interpreter.md`, `captions.md` — and the
spec text in `docs/architecture.md`

**What:** Establish a named, reusable protocol that any agent can adopt. The protocol is
three steps:

```
## Step 0: Resolve output path

1. Read `AGENTS.md` (project root). Look for the `## Documentation namespace` section.
2. Find the routing table entry that matches this agent's output kind.
3. Use that path as the target directory.

Fallback: if AGENTS.md has no Documentation Namespace section, or the routing
table has no entry for this kind, use the default path declared in this agent's
`outputs.default_path` frontmatter field. Proceed silently — do not stop or ask.

Never write to the project root. Never hardcode a path that includes a project
name, user name, or absolute filesystem location.
```

This is inserted at the top of the Workflow section of any agent that writes dated
artifacts. It is NOT a library call; it is literally prose in the agent body that the
LLM executing the agent reads and acts on.

**Why:** This is the architectural keystone of phase 7. It satisfies three constraints
simultaneously:

1. **Generalizable**: no project-specific knowledge baked into the agent definition in
   SciAgent-toolkit
2. **Convention-following**: reads the project's declared routing at runtime — the same
   table the human and every other agent read
3. **Friction-free**: the fallback to `outputs.default_path` means the agent works
   correctly even in a project that hasn't set up a managed block yet; no stopping, no
   asking

The alternative — hardcoding the canonical path (`docs/_internal/ai-generated/sessions/`)
directly in the agent body — would work for the default convention but silently break
when a project legitimately deviates. The alternative — requiring the caller to pass the
path on every invocation — puts knowledge in the wrong place (caller) and adds friction
to the most common case (single user, no deviation). Step-0 adds zero friction when
the convention is set, and graceful defaults when it isn't.

**How:**
1. Write the Step-0 block as a named section in `docs/architecture.md` so it can be
   referenced: "Any agent that writes dated artifacts must implement the Step-0
   convention-reading protocol. See `docs/architecture.md § Agent output contract.`"
2. Insert the Step-0 block into `handoff.md`, `bio-interpreter.md` (see 7.3, 7.4).
3. For `captions.md`: add a lighter variant — the caption target is always
   `03_results/{phase}/README.md` where `{phase}` is provided by the caller; Step-0
   is not needed, but the agent should assert that the target matches the phase-4
   schema (section 4.9 format) before writing.

---

### 7.3 Rewrite the handoff agent body

**Files:** `agents/analysis-base/handoff.md`

**What:** Full body rewrite (not just a path change). The current body has five structural
problems beyond the wrong path.

**Changes:**

1. **Filename format**: `YYYY-MM-DD_<slug>.md`. Add time suffix only on same-day collision
   (ADR-3.5 recommendation). Slug = a 2–4 word description of what the session did
   (`integration-working`, `qc-threshold-tuning`, `handoff-agent-rewrite`).

2. **Output path**: resolved via Step-0 protocol from AGENTS.md routing table.
   Default: `docs/_internal/ai-generated/sessions/`. No `.handoff_archive/` — dated
   files accumulate in place; `ls -t docs/_internal/ai-generated/sessions/` is the
   manifest; newest is current.

3. **Context files the agent reads** (Step 1 — Gather Context):
   - `context.md` — the pointer (30-line orientation)
   - `docs/_internal/scientific-context.md` — biology, hypotheses
   - `docs/_internal/ai-generated/sessions/<most-recent file>` — prior handoff for
     continuity (not `plan.md`, which doesn't exist)
   - `tasks.md` if it exists (user-maintained, optional)

4. **Pre-write artifact check** (Step 2, per phase-3 section 3.4 step 5):
   Before writing the session file, scan `03_results/` for any artifact files
   (`*.pdf`, `*.png`, `*.svg`, `*.csv`, `*.tsv`, `*.html`) with no matching
   `## <filename.ext>` heading in the sibling phase `README.md`. List missing captions
   under `## Uncaptioned artifacts` in the handoff. This is a flag, not a blocker —
   the handoff is written regardless; the list gives the next session a concrete
   starting task.

5. **Section template**: simplify from the current heavily-formatted template. The
   existing template is 80+ lines of mandatory headings that make the agent produce
   bureaucratic documents. Reduce to:
   - Quick Orientation: stage, last completed, next step (≤ 5 lines)
   - What happened this session: concrete (code, files, decisions)
   - Technical state: working checkpoints (paths + sizes), any blockers
   - Uncaptioned artifacts (Step-2 output)
   - Next session: 2–3 prioritized next steps with exact commands if applicable

   Target: 30–60 lines of markdown, readable in 3 minutes, not 15.

**Why:** The current handoff produces 80+ line bureaucratic documents because the template
front-loads mandatory sections that are often empty. Lean documents get read; comprehensive
documents get skipped. The `context.md` read was pointed at `plan.md` which is a leftover
from a different project template; reading the actual phase-3 context files makes the
handoff properly oriented.

**How:**
1. Replace the entire body after the YAML frontmatter separator.
2. Update `description:` in frontmatter to remove references to `.handoff_archive/`.
3. Add `outputs.default_path: docs/_internal/ai-generated/sessions/` to frontmatter.
4. Add `domain: [session-management]` to frontmatter.

---

### 7.4 Fix bio-interpreter output path

**Files:** `agents/analysis-base/bio-interpreter.md`

**What:**

1. Add Step-0 protocol to the body: resolve the research output path from AGENTS.md
   routing table ("Literature / web research" entry). Default:
   `docs/_internal/ai-generated/research/`.

2. Change output filename from `research_notes.md` to
   `YYYY-MM-DD_NN_<topic-slug>.md` per the naming convention (section 7.1). The topic
   slug derives from the genes/pathway/question the user gave (2–4 words, snake_case).

3. Update the description's handoff pattern note (currently says
   "Mechanism explanation in `research_notes.md`") to reflect the new path.

4. Add to frontmatter:
   ```yaml
   domain: [literature-research]
   outputs:
     default_path: docs/_internal/ai-generated/research/
     kind: research-note
   ```

**Why:** `research_notes.md` at the project root is exactly the root-clutter pattern
identified in the phase-3 analysis (the AdaW project had `research_notes.md`, `ideas.md`,
`docs_notes.md` all at root). The fix is simple; the only reason it wasn't done yet is
that the agent predates the convention.

**How:** Edit the "Synthesize Research into Structured Documentation" section of the body
to change the target path and filename. Keep all other behaviour.

---

### 7.5 Add lightweight capability metadata to frontmatter

**Files:** All agent `.md` files in `agents/`

**What:** Two optional fields added to YAML frontmatter. Optional — absence does not
break anything; absence only means the roster table (section 7.6) shows a blank for
that column.

```yaml
domain:                        # list of concise taxonomy tags (1–3 tags)
  - session-management
outputs:                       # where this agent writes, and how the path is resolved
  default_path: docs/_internal/ai-generated/sessions/
  kind: session-handoff        # matches the routing-table "You are writing…" column
  path_source: AGENTS.md       # constant; signals Step-0 applies
```

**Canonical domain tags** (open list — extend as new agents are added):

| Tag | Agents |
|---|---|
| `session-management` | handoff |
| `literature-research` | bio-interpreter |
| `data-exploration` | insight-explorer |
| `code-review` | code-reviewer |
| `documentation` | doc-curator, captions |
| `tool-lookup` | docs-librarian |
| `architecture` | architect, mapper, meta-architect, slicer |
| `review` | bioinf, wetlab, graphic, stat, divergent, ml |
| `synthesis` | synth, feature-reviser |
| `status` | status-reporter |
| `orchestration` | meta-architect |

For agents that don't write files (`status-reporter`, `synth` in pure-read mode,
`docs-librarian`): omit `outputs:` entirely.

**Why:** Two consumers:
1. `sciagent roster` (section 7.6): reads frontmatter to build the table and manifest.
   Without domain tags the table exists but the taxonomy column is blank — still useful,
   just less structured.
2. An orchestrator reading `.claude/agents/` at session start: can `grep "^domain:" *.md`
   to get a coarse map before reading full descriptions. Domain tags are faster to scan
   than prose descriptions when the roster is large.

**Why lightweight / why optional:** Mandatory schema validation would make adding a new
agent a two-step process (write the agent, satisfy the schema). We are not trying to build
a type system — we are trying to make the roster scannable. Optional fields added
incrementally as agents are touched for other reasons is the right pace.

**How:**
1. During the implementation sprint for this phase, add `domain:` and `outputs:` to the
   seven analysis-base agents (handoff, bio-interpreter, captions, doc-curator,
   docs-librarian, insight-explorer, code-reviewer) since they are being rewritten anyway.
2. Add to architect agents opportunistically but do not block the sprint on it.
3. Update `CLAUDE.md` in the toolkit (the developer guide) to document the frontmatter
   schema, including that `domain:` and `outputs:` are optional.

---

### 7.6 Add `sciagent roster` command

**Files:** `lib/sciagent/roster.sh` (new), `bin/sciagent` (add roster subcommand)

**What:** New CLI subcommand `sciagent roster` that reads the active agent set from
`.claude/agents/` and prints a human-readable table plus optionally generates a JSON manifest.

**Table output** (human-readable, always printed):
```
AGENT             DOMAIN                  MODEL    DESCRIPTION (first 72 chars)
handoff           session-management      sonnet   Create session handoff after completing …
bio-interpreter   literature-research     sonnet   Research biological mechanisms via litera…
captions          documentation           sonnet   Generate figure legend captions for 03_re…
…
```

**JSON output** (written to stdout with `--json` flag, or to `.claude/agent-manifest.json`
with `--write`):
```json
{
  "generated_date": "YYYY-MM-DD",
  "role_stack": ["base"],
  "agents": [
    {
      "name": "handoff",
      "model": "sonnet",
      "domain": ["session-management"],
      "description_brief": "Create session handoff after completing work",
      "outputs": {
        "kind": "session-handoff",
        "default_path": "docs/_internal/ai-generated/sessions/"
      }
    }
  ]
}
```

**Implementation:** Pure bash + awk, no external dependencies (`yq` not assumed).
Frontmatter extraction: `awk 'BEGIN{in_fm=0} /^---/{in_fm++; next} in_fm==1{print; next} in_fm>1{exit}' agent.md`
Field extraction: `grep "^name:" | cut -d: -f2- | tr -d ' "'`

**Why:** Three consumers:
1. **Human in terminal**: `sciagent roster` gives a one-line-per-agent orientation
   faster than `ls .claude/agents/` followed by reading 19 agent files.
2. **Orchestrator bootstrapping**: a subagent or workflow script can call
   `sciagent roster --json` to load the roster without reading individual files.
   This is especially useful when the active role stack isn't known ahead of time.
3. **Session start ritual**: `sciagent status` shows the role stack; `sciagent roster`
   shows what that stack provides. Together they give full orientation.

**Why NOT committed as a static file**: a static `.claude/agent-manifest.json` committed
to the repo would go stale whenever agents change, become a required update step for
contributors, and need a gitignore decision. Generated on demand with `--write` when an
orchestrator needs it is better: always current, no maintenance burden.

**How:**
1. Write `lib/sciagent/roster.sh`:
   - Accept `--json` and `--write` flags
   - Exit gracefully if `.claude/agents/` does not exist (with message: "no active role — run `sciagent activate <role>`")
   - Truncate description to first 72 chars for table column; use full description in JSON
2. Register in `bin/sciagent`:
   ```bash
   roster) source "$LIB_DIR/roster.sh"; cmd_roster "$@" ;;
   ```
3. Add to `--help` output and `README.md`.
4. Add test: `tests/test_roster.sh` — activates base role, calls roster, asserts handoff
   appears in output.

---

### 7.7 Overhaul doc-curator with phase-3/4 checks

**Files:** `agents/analysis-base/doc-curator.md`

**What:** The existing body is good for general documentation curation but predates the
phase-3/4 conventions. Add three explicit checks as a named phase in the workflow.

**New Phase: Convention compliance checks** (insert before the existing "Consolidation"
phase):

**Check C1 — File size discipline:**
Flag any of the following over their limits:
- `AGENTS.md` > ~150 lines: recommend decomposing project-specific content to `docs/_internal/`
- `CLAUDE.md` > 2 substantive lines (anything beyond `@AGENTS.md` + a comment): flag as drift
- `context.md` that is not a pointer (> 30 lines, or does not contain a link to `scientific-context.md`): recommend converting to pointer + creating `docs/_internal/scientific-context.md`

**Check C2 — One-way reference rule:**
Scan all public-facing files:
- `03_results/**/README.md`
- `docs/guides/**`
- `docs/*.md` (outside `_internal/`)

For any file containing a path or link with `_internal/` in it: report as a violation
with the file and line number. Suggest rewriting to state the outcome rather than cite
the internal reasoning doc.

**Check C3 — Uncaptioned artifacts:**
```bash
find 03_results -name "*.pdf" -o -name "*.png" -o -name "*.svg" \
     -o -name "*.csv" -o -name "*.tsv" -o -name "*.html" | \
while read f; do
  dir=$(dirname "$f")
  base=$(basename "$f")
  if ! grep -q "^## $base" "$dir/README.md" 2>/dev/null; then
    echo "UNCAPTIONED: $f"
  fi
done
```
Report all uncaptioned files. This is not an error that blocks the curator — it is a
finding for the human to address (or to assign to the `captions` agent).

**Description field rewrite:** Convert the current prose-with-`\n`-escapes description
to a YAML block scalar (`|`) so it renders correctly in the source file. Trim to the
3-part template: "Use this agent when [trigger]. Do NOT use for [exclusions]."

**Add to frontmatter:**
```yaml
domain: [documentation]
```

**Why:** These three checks were specified in the phase-3/4 plan docs as additions to
doc-curator.md. The agent body predates those specs. This phase closes the gap. The
checks are not new ideas — they are the enforcement arm of conventions already decided.

**How:**
1. Insert "## Phase N: Convention compliance" as the first phase of the body workflow.
2. Add the three check implementations as bash snippets the agent should run (or adapt
   for its tool calls — doc-curator has `Bash` access via its existing tool list).
3. Update description frontmatter.

---

## Open ADRs

### ADR-7.1: How strict should the Step-0 fallback be?

**Options:**
- A — Use frontmatter `outputs.default_path`, proceed silently (no user interrupt)
- B — Use the frontmatter default, but emit a visible notice: "AGENTS.md has no routing
  table; writing to default path [path]. Run `sciagent activate <role>` to configure."
- C — Halt and ask the user to configure AGENTS.md before proceeding

**Recommended:** A for new artifacts; B is acceptable if the agent notices AGENTS.md
exists but has no Documentation Namespace section (suggests incomplete setup rather than
a new project). C is rejected: stopping an agent mid-session to ask for configuration
is exactly the friction we are trying to eliminate.

**Blocking implementation:** no

### ADR-7.2: Committed manifest vs generated on demand

**Options:**
- A — `sciagent roster --write` generates `.claude/agent-manifest.json` on demand; not
  committed; used by orchestrators that call it fresh
- B — Manifest committed alongside activation artifacts (updated whenever role changes)
- C — No JSON manifest at all; roster table only

**Recommended:** A. A committed manifest goes stale (agents change, manifest doesn't),
adds a required update step to `sciagent activate`, and needs a gitignore decision in
every project. Generate on demand is always current. C is fine for now (phase 7 can
ship without the manifest if no concrete external orchestrator consumer exists yet);
the roster table still achieves the human-orientation goal.

**Blocking implementation:** no

### ADR-7.3: Where does insight-explorer write its output?

The insight-explorer description says "Data patterns + viz recommendations" — it is
unclear whether it writes a file or only produces conversational output. The routing table
in AGENTS.md.template has no entry for data-exploration findings.

**Options:**
- A — insight-explorer writes `docs/_internal/ai-generated/explorations/YYYY-MM-DD_NN_<what>.md`
  (fits the "codebase / data explorations" routing table entry)
- B — insight-explorer is conversational only; no file output, no Step-0
- C — Add a new routing table entry "Data exploration findings" if A

**Recommended:** A and C together. Data exploration findings have the same durability
need as codebase explorations — losing them to the session context is exactly the
observability loss the namespace solves. The routing table entry in the template should
be expanded from "codebase / data exploration" to "codebase / data / result exploration"
to cover this agent's output.

**Blocking implementation:** no — can be done in the phase-3 template pass alongside 7.4.

### ADR-7.4: Should `sciagent activate` auto-run `sciagent roster --write`?

**Options:**
- A — Activate always regenerates the manifest (always current after activate)
- B — Activate does not touch the manifest; roster is explicit only
- C — Activate regenerates only if `--manifest` flag is passed

**Recommended:** B for now. Adding manifest generation to every activate increases the
surface area of a critical command. Once a concrete orchestrator consumer exists and
manifest freshness matters, promote to A.

**Blocking implementation:** no

---

## Dependencies

- **Depends on Phase 3**: AGENTS.md.template must include the Documentation Namespace
  section with the routing table. The Step-0 protocol reads that table. If phase 3
  is not yet landed, the Step-0 fallback (frontmatter `default_path`) must carry all
  agents through without AGENTS.md routing.
- **Depends on Phase 4**: Caption convention in `03_results/{phase}/README.md` (section
  4.9) is referenced by the handoff agent's pre-write check. If phase 4 hasn't landed
  yet, the check is a no-op (no README.md stubs exist to compare against).
- **Enables**: any future workflow orchestration that dispatches agents by domain tag or
  reads the roster manifest. Enables `doc-curator` enforcement (C1–C3) to be meaningful
  once the docs namespace and caption convention exist in projects.

---

## Breaking Changes

- **handoff agent**: output path changes from project root + `.handoff_archive/` to
  `docs/_internal/ai-generated/sessions/`. Old root `handoff_*.md` files in existing
  projects are not auto-migrated; they remain at root and can be manually moved.
- **bio-interpreter**: writes `YYYY-MM-DD_NN_<topic>.md` to `docs/_internal/ai-generated/research/`
  instead of `research_notes.md` at root. Any existing `research_notes.md` files are not
  touched.
- **`sciagent` CLI**: gains a new `roster` subcommand; no existing subcommands change.

All changes are acceptable per project ground rules (single user, no downstream consumers
of the agent API, clean modularity preferred over backwards compatibility).

---

## Estimated Scope

| Item | Files | Approx. delta |
|---|---|---|
| Naming conventions in `_internal/README.md` (7.1) | 1 | +20 lines |
| Step-0 protocol in `docs/architecture.md` (7.2) | 1 | +25 lines |
| Handoff agent full rewrite (7.3) | 1 | ~ −80 / +60 lines (body replacement) |
| bio-interpreter path fix + Step-0 (7.4) | 1 | −5 / +20 lines |
| Frontmatter metadata additions (7.5) | 7 analysis-base agents + all architect agents | +4–6 lines per file |
| `lib/sciagent/roster.sh` new (7.6) | 1 new | ~80 lines |
| `bin/sciagent` roster dispatch (7.6) | 1 | +3 lines |
| `tests/test_roster.sh` (7.6) | 1 new | ~30 lines |
| doc-curator overhaul (7.7) | 1 | +60 lines |
| AGENTS.md.template pointer + naming ref (7.1) | 1 | +2 lines |

**Total:** ~12 files touched/created; net ≈ +280 / −90 lines. No changes to the
activation pipeline (`lib/sciagent/activate.sh`, `symlinks.sh`, `stack.sh`) — the
roster command reads the already-activated `.claude/agents/` directory without touching
the activation mechanism.
