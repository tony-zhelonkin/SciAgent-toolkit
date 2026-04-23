---
parent: ./README.md
view: lifecycle-proposal
status: PROPOSED
date: 2026-04-22
---

# Pipeline lifecycle — proposed refinements

This document is a design proposal, not an accepted contract. It names the gaps in the current architect harness lifecycle and proposes five refinements (R1–R5, plus an optional R6). It is filed here — alongside `02-behavior.md` and `05-meta-approach.md` — so the user can reason about scope before any of it is implemented.

Nothing in this doc changes the harness. Acceptance would mean:
1. The ADRs in `03-decisions.md` are extended with one ADR per accepted refinement.
2. The command prompts in `../../commands/` and agent prompts in `../../agents/` are edited per the "Files modified" list below.
3. `02-behavior.md` absorbs the new emissions into its stage descriptions.

---

## 1. Problem statement

### 1.1 Triggering observation

On 2026-04-22 the meta-architect routed scope-triage items from a round-2 `leiden-stability` synthesis:

- A `gene_namespace` block was added to MADR-008 in `docs/_meta/design.md`.
- A `cluster_labels` sidecar schema was added to the same MADR.
- A label-drift CI guard was filed as `D-026` in `docs/_meta/deferred.md`.
- A `### Revision log` section was appended inside MADR-008 noting the additive expansion.

The plumbing worked. Every item landed in the right file. But two questions from the user exposed a discoverability gap that the raw material does not close:

1. **"Where is this MADR-008 revision log?"** The command described what it did but did not cite the anchor. The user had to grep.
2. **"The whole life cycle is not defined … what is WIP, what is the backlog, where are complete tasks? The production pipeline is not defined."** `/status` exists and reads frontmatter + grep, but its output evaporates with the chat. `_meta/deferred.md` exists as the canonical backlog, but it is discoverable only if you already know to look. No file on disk answers "what state is the portfolio in right now."

### 1.2 What the harness already provides

Confirmed by reading `02-behavior.md`, `commands/status.md`, `agents/status-reporter.md`, and the live state under `docs/_meta/` + per-feature `docs/{feature}/`:

| Surface | File(s) | What it tracks | Durability |
|---|---|---|---|
| Centralised backlog | `_meta/deferred.md` | One monotonic table, 26 rows today. Columns: `ID | Source | Source doc | Description | Reason | Cost | Depends | Status`. Status values `DEFERRED | IN-CONSIDERATION | REJECTED | ACTIVATED`. Append-only. | Durable |
| Cross-feature decisions | `_meta/design.md` | 11 MADRs with frontmatter `status: PROPOSED \| APPROVED` + per-MADR inline Status. Contains "Deliberate meta-level rejections" for permanent (non-time-boxed) no-gos. | Durable |
| Per-feature decisions | `{feature}/design/03-decisions.md` | ADRs (ACCEPTED), Open Questions (`OQ-N`), Hard constraints, "Deliberate departures from synthesis/audit". Deferrals in the last section also appended to `_meta/deferred.md` with provenance. | Durable |
| Design readiness | `{feature}/design/*.md` frontmatter | `status: DRAFT \| APPROVED`, flipped by `/design` Phase 4 and `/meta-apply` Phase 6; rolled back by architect `NEEDS ITERATION`. | Durable |
| Architect verdict | `{feature}/design/review.md` | `Verdict: READY \| NEEDS ITERATION \| NEEDS DISCUSSION`. Overwritten on re-run. | Durable-ish |
| Iteration archives | `{feature}/review/.history/round-N/*.md`, `{feature}/synthesis/.history/round-N.md` | Prior-round reviews/syntheses, pinned by `commands/review.md` ADR-019. | Durable |
| Portfolio status | `/status` command via `status-reporter` Sonnet sub-agent | Per-feature stage + verdict + blocker + cross-feature conflicts. Chat-only output. | **Ephemeral** |

### 1.3 What is missing

1. **No durable portfolio dashboard.** `/status` emits a markdown table to chat; nothing on disk answers "what is the state of this portfolio right now." After the chat scrolls, the view is gone.
2. **No pinned location for "Revision log" inside a MADR or a design doc.** MADR-008 has one because `meta-architect` chose to add it during the 2026-04-22 edit. The other ten MADRs do not. The convention is ad-hoc — a reader cannot predict whether a given MADR has a revision log, or where in the MADR it lives.
3. **No SHIPPED / DONE state.** Completion is inferred from git history and absence of blockers. Six months from now, nobody will remember which features shipped vs which stalled at `plan/` without a durable marker per feature.
4. **No `_meta/README.md`.** The `_meta/` directory is the portfolio control surface, but there is no index that names the lifecycle, the ID conventions (`D-NNN`, `MADR-NNN`, `OQ-M-N`), and the purpose of each file within. New readers land in `_meta/design.md` without context.
5. **Commands do not cite anchors they wrote to.** The meta-architect's 2026-04-22 output described the three edits ("Added `gene_namespace` block", "Added `cluster_labels` block", "Filed D-026", "Revision-log note added") but did not print `file:line` or `§ section-heading` anchors for any of them. The user had to ask "where is this MADR-008 revision log?" and grep to find it.

**What this is not.** These gaps are discoverability + durability failures, not mechanism failures. The backlog exists and is centralised. The decisions exist and are durable. The lifecycle states (DRAFT → APPROVED → READY) exist on design docs. The refinements below do not introduce new mechanisms; they make the existing ones visible, predictable, and queryable without a grep.

---

## 2. Lifecycle states (proposed canonical model)

This is the model the refinements below encode. A feature moves through seven states, each defined by a durable artifact or frontmatter flag, not by a human judgement:

```
┌───────┐   /map      ┌────────┐   /review    ┌──────────┐   /synthesize   ┌──────────────┐
│  NEW  │────────────▶│ MAPPED │─────────────▶│ REVIEWED │────────────────▶│ SYNTHESIZED  │
└───────┘             └────────┘              └──────────┘                 └──────────────┘
                                                                                  │
                                                                                  │ /design
                                                                                  ▼
┌──────────┐  /implement completes  ┌──────────────┐  /plan   ┌─────────┐ architect gate  ┌──────────┐
│ SHIPPED  │◀──────all phases done──│ IMPLEMENTING │◀─────────│ PLANNED │◀────READY───────│ DESIGNED │
└──────────┘                        └──────────────┘          └─────────┘                 └──────────┘
```

Detection rule for each state (what makes it observable without asking):

| State | Artifact that makes it true |
|---|---|
| NEW | feature directory does not exist |
| MAPPED | `{feature}/map.md` exists |
| REVIEWED | `{feature}/review/*.md` has ≥1 file |
| SYNTHESIZED | `{feature}/synthesis.md` exists |
| DESIGNED | `{feature}/design/README.md` frontmatter `status: APPROVED` AND `{feature}/design/review.md` verdict is READY |
| PLANNED | `{feature}/plan/phase-NN.md` files exist with all phases enumerated |
| IMPLEMENTING | ≥1 `plan/phase-NN.md` carries `completed:` in frontmatter but not all |
| SHIPPED | All `plan/phase-NN.md` carry `completed:` AND `{feature}/design/README.md` frontmatter `status: SHIPPED` |

The rules are mechanical — no state requires human interpretation. This is the invariant the refinements enforce.

---

## 3. Proposed refinements

Five refinements (R1–R5) form a coherent bundle. R6 is optional polish.

### R1 — Pin "Revision log" as a required section inside every MADR

**Rule.** Every MADR in `_meta/design.md` ends with a `### Revision log` section. Empty on creation; every subsequent edit appends one row. Same rule mirrors into per-feature `design/03-decisions.md` for every ADR that gets edited post-acceptance.

**Shape.**
```markdown
### Revision log
| Date | Source | Change |
|------|--------|--------|
| 2026-04-22 | leiden-stability round-2 synthesis | Added `gene_namespace` + `cluster_labels` blocks (additive; status unchanged). |
```

**Who enforces.** The `meta-architect` agent prompt adds: "Every MADR you emit or edit ends with a `### Revision log` section. If it doesn't exist, create it empty. If you are editing a MADR, append a row." The `architect` consistency gate checks this on every `/design` and `/meta-design` verdict.

**Backfill.** On next `/meta-design --iterate`, existing MADRs without a revision log get an empty section added. No historic reconstruction attempted.

**Why this and not a top-of-file CHANGELOG.** A CHANGELOG concentrates all change history in one place far from the thing that changed. A per-MADR revision log sits inside the object it describes — the reader always knows where to look. This mirrors how ADRs work in larger systems: the decision record *includes* its own history.

**Cost.** Small — one block added to the `meta-architect` prompt, one check added to the `architect` gate.

### R2 — Durable portfolio dashboard: `docs/_meta/status.md`

**Rule.** `/status` writes its output to `docs/_meta/status.md`. The current chat-only output becomes a terminal render of the same content. The file is the canonical view; the chat is a convenience.

**Sections (pinned).**
```markdown
---
generated: 2026-04-22T14:03:00Z
by: /status --write
---

# Portfolio status

## Pipeline
| Feature | State | Design status | Verdict | Open blocker | Last touched |
|---------|-------|--------------|---------|--------------|--------------|
| biological-workflow | DESIGNED  | APPROVED | READY            | —              | 2026-04-18 |
| leiden-stability    | SYNTHESIZED | —      | —                | 8 P0 ADRs pending /design | 2026-04-22 |
| normalisation       | PLANNED   | APPROVED | READY            | —              | 2026-04-15 |
| theme-bundles       | DESIGNED  | APPROVED | READY            | —              | 2026-04-10 |
| umap                | IMPLEMENTING | APPROVED | READY         | phase-04 in progress | 2026-04-21 |

## Shipped
(none yet)

## Deferred backlog (snapshot)
5 most recent; 22 total DEFERRED, 3 ACTIVATED, 1 REJECTED.
| ID | Source | Description | Status |
|----|--------|-------------|--------|
| D-026 | leiden-stability | Label-drift CI guard | DEFERRED — owner unassigned |
| D-025 | normalisation | Golden-file HTML snapshot test | DEFERRED |
| ... | | | |
Full table: `docs/_meta/deferred.md`.

## Cross-feature blockers
| ID | Text | Blocks | Status |
|----|------|--------|--------|
| OQ-M-6 | Namespace authority for Ochiai | leiden-stability /design, umap Round 2 | open |
| ... | | | |
Full list: `docs/_meta/map.md § Open cross-feature questions` and `docs/_meta/design.md § Open inter-feature questions`.
```

**Who writes.** The `status-reporter` agent is extended to accept a `--target` arg (default `docs/_meta/status.md`). It still reads only frontmatter + grep — never opens full file bodies.

**Refresh semantics.** Regenerating overwrites the file atomically (write to `.status.md.tmp`, rename). `generated:` timestamp in frontmatter is the freshness signal.

**This file is a view, not a source.** It is generated, not edited by hand. `deferred.md` and `design.md` and per-feature frontmatter remain canonical; `status.md` is the projection. If someone edits `status.md` directly, the next `/status` overwrites it. The frontmatter line `by: /status --write` makes this visible.

**Why not just leave it to chat.** Because the user asked "where do I look to see what is going on?" and the honest answer today is "run `/status` again and hope it catches everything." A single file on disk that always answers the question is worth one sub-agent invocation per refresh.

### R3 — Explicit SHIPPED state

**Rule (two parts).**

*Part 3a: `/implement` marks phase completion.* On successful phase execution (all checklist items pass), `/implement` appends `completed: YYYY-MM-DD` to the frontmatter of `{feature}/plan/phase-NN.md`. If phase execution halts (uncertainty, test failure), no completion marker is written.

*Part 3b: Feature flips to SHIPPED automatically.* When all `phase-NN.md` files in `{feature}/plan/` carry `completed:`, `status-reporter` flips `{feature}/design/README.md` frontmatter from `status: APPROVED` to `status: SHIPPED`, and the feature moves to the **Shipped** section of `status.md` next refresh.

**Discussion.** The auto-flip could alternatively be manual (a new `/ship` command, R6 below). Auto is preferred because: (a) it removes a step the user would forget, (b) the detection rule is mechanical — all phases carry `completed:` is a deterministic check, (c) manual rollback is easy (edit the frontmatter).

**What SHIPPED means.** "The plan, as approved, has been executed end-to-end." It does NOT mean "bug-free" or "production-ready" or "no future work needed." Future work on a shipped feature routes through `/design` again, which cannot re-open a SHIPPED doc — it creates a new feature slug or a new design round with its own artifacts. This is the invariant that makes SHIPPED durable.

**Cost.** Small — one frontmatter-append line in `/implement` on phase success, one check in `status-reporter`.

### R4 — Commands must cite concrete anchors they wrote to

**Rule.** Every command that appends content to a shared file prints an anchor line in its chat output with `file:line` or `file § section` precision. No descriptions — anchors.

**Shape.**
```
Wrote:
  docs/_meta/deferred.md § D-026 (line 39)
  docs/_meta/design.md § MADR-008 § Revision log (line 403)
  docs/_meta/design.md § MADR-008 § Consequences per feature (line 385)
```

**Who enforces.** The prompt addition lives in:
- `commands/meta-design.md` — meta-architect edits to MADRs.
- `commands/meta-apply.md` — feature-reviser edits to per-feature `03-decisions.md` + deferrals.
- `commands/design.md` — design agent edits to per-feature design docs + deferrals.

Agents (`meta-architect`, `feature-reviser`) inherit the rule in their system prompts.

**Why anchors and not descriptions.** The 2026-04-22 transcript shows a thorough description ("MADR-008 revision log notes the 2026-04-22 additive expansion — no status change; still ACCEPTED") followed by the user asking "where is this MADR-008 revision log?" The description answered *what* happened; the anchor would have answered *where*. Both are useful; both should be present; today only one is.

**Implementation note.** Agents already know line numbers after they write (they just called Edit/Write). This is surfacing existing knowledge, not computing new knowledge.

**Cost.** Trivial — one rule in four prompts.

### R5 — `docs/_meta/README.md` as the lifecycle index

**Rule.** On first `/meta-map` (when `_meta/` is first created), a `_meta/README.md` is emitted from a toolkit template. Thereafter it is read-only from the harness's perspective — users edit it if they want to add project-specific conventions.

**Contents (template).**
```markdown
# Portfolio control surface — _meta/

This directory holds the cross-feature artifacts that govern work across
≥2 features. Per-feature artifacts live in `docs/{feature}/`.

## Lifecycle

[Mermaid diagram of the seven states from §2.]

## Files in this directory

| File | Role | Written by | Canonical? |
|------|------|------------|------------|
| map.md | Portfolio inventory: shared touch-points, convergent concerns, cross-feature open questions | /meta-map | Yes |
| design.md | MADRs (meta-ADRs), inter-feature decisions inherited by every /design | /meta-design | Yes |
| deferred.md | Append-only backlog of time-boxed deferrals across all features | /design, /meta-design, /meta-apply | Yes |
| plan.md | Portfolio-phase sequencing, collision check, parking items | /meta-plan | Yes |
| status.md | Dashboard view — regenerated, not edited | /status --write | View (not source) |
| README.md | This file | /meta-map (first run, template) | Manual thereafter |

## ID conventions

- `MADR-NNN` — meta-ADR in `design.md`
- `D-NNN` — deferred item in `deferred.md`
- `OQ-M-N` — open cross-feature question in `map.md` or `design.md`
- `ADR-NNN` — per-feature ADR in `{feature}/design/03-decisions.md`
- `OQ-N` — per-feature open question in `{feature}/design/03-decisions.md`

## Where to look

- "What is the state of things right now?" → `status.md`
- "What is in the backlog?" → `deferred.md`
- "What have we decided across features?" → `design.md`
- "What depends on what?" → `plan.md` (dependency graph)
- "What does each feature touch?" → `map.md`
```

**Why a README and not a single "control panel" command.** A command evaporates; a README persists and is the first thing a new reader (or returning user) sees when they `ls docs/_meta/`. The command (`/status`) regenerates the dashboard view; the README orients the reader to the directory that contains both.

**Cost.** Template file + three lines in the `/meta-map` command to emit it on first run.

### R6 — `/ship <feature>` (optional, deferrable)

A thin command that:
1. Verifies all `plan/phase-NN.md` carry `completed:`.
2. Runs `/verify` and requires CLEAN.
3. Requires zero open P0 items in `03-decisions.md`.
4. Flips `design/README.md` frontmatter to `status: SHIPPED`.
5. Writes an entry to `_meta/shipped.md` (new file) with feature name, ship date, and a one-line summary pulled from `design/README.md § Approach summary`.

R3 already does the flip automatically. R6 is for users who want an explicit confirmation step + a durable shipped ledger separate from `status.md`. **Defer until R1–R5 are proven.**

---

## 4. Files to modify (if adopted)

### Command prompts
- `commands/meta-design.md` — R1 (MADR template adds Revision log section), R4 (output anchors).
- `commands/meta-apply.md` — R4 (feature-reviser cites anchors for every edit).
- `commands/design.md` — R1 (per-ADR revision logs), R4 (output anchors for deferrals + design doc edits).
- `commands/status.md` — R2 (add `--write`, default target `docs/_meta/status.md`), R3 (recognise SHIPPED state).
- `commands/implement.md` — R3 (append `completed:` to phase frontmatter on success).
- `commands/meta-map.md` — R5 (emit `_meta/README.md` from template on first run).

### Agent prompts
- `agents/meta-architect.md` — R1 (every MADR emits with Revision log), R4 (cite anchors).
- `agents/feature-reviser.md` — R1 (ADR revision logs), R4 (cite anchors).
- `agents/status-reporter.md` — R2 (write target, markdown schema), R3 (SHIPPED detection + frontmatter flip).
- `agents/architect.md` — R1 (consistency gate checks for Revision log section in every MADR and every edited ADR).

### Templates (new files)
- `templates/_meta-README.md.template` — R5.

### Docs
- `docs/workflows/architect/02-behavior.md` — add emissions per R1/R2/R3 to every relevant stage description.
- `docs/workflows/architect/03-decisions.md` — one ADR per accepted refinement, recording the decision + trade-offs.

### No edits to user project `docs/_meta/` or per-feature docs. Those are written by the updated commands on next run.

---

## 5. Migration plan

1. **Accept R1–R5 as ADRs** in `03-decisions.md`.
2. **Edit command + agent prompts** per §4. Run the updated `/meta-design --iterate` on an existing `_meta/design.md` — every MADR grows an empty Revision log. Existing MADR-008 keeps its.
3. **Run `/status --write`** — `docs/_meta/status.md` appears. Confirm the detection rules in §2 map to the live state of the five in-flight features.
4. **Run `/meta-map`** — if `_meta/README.md` is absent, template is emitted. If present, leave untouched.
5. **Test R3 with a manual frontmatter edit**: add `completed: 2026-04-22` to one `normalisation/plan/phase-NN.md`, re-run `/status --write`. Confirm the phase registers. Add it to the remaining phases; confirm `normalisation` moves to SHIPPED in `status.md` and its `design/README.md` frontmatter flips.
6. **Test R4 by simulating a `/meta-apply`** — confirm the output contains anchor lines, not just descriptions.

---

## 6. Open questions

### OQ-L-1 — Auto-flip SHIPPED or manual (R3 vs R6)?
Auto-flip is mechanical and removes a step; manual is explicit and creates a confirmation gate. Recommendation: auto-flip (R3), reserve `/ship` (R6) for later if users report wanting the gate.

### OQ-L-2 — `status.md` refresh cadence?
Regenerating on every `/status` invocation is the minimum. Should it also be regenerated automatically at the end of `/design`, `/meta-apply`, `/implement`? Pro: always fresh. Con: adds latency and sub-agent cost to common commands. Recommendation: on-demand only for now; add auto-refresh hooks only if users report staleness.

### OQ-L-3 — Where do permanent rejections surface in `status.md`?
Today `_meta/design.md § Deliberate meta-level rejections` and per-feature `03-decisions.md § Deliberate departures (permanent)` are durable but not dashboarded. `status.md` could include a "Declined (no time horizon)" section. Recommendation: add in a follow-up — keep the first dashboard version small.

### OQ-L-4 — What about OQs that never close?
Open Questions (`OQ-N`, `OQ-M-N`) are durable but have no intrinsic lifecycle. An OQ open for six months is either resolved-but-not-closed, or deprioritised. Recommendation: `status.md § Cross-feature blockers` shows age; ageing OQs become visible.

---

## 7. Scope decision (for the user)

Three coherent slices, each independently shippable:

1. **All-in (R1 + R2 + R3 + R4 + R5).** Full lifecycle contract. Roughly one working day of toolkit edits + iteration. Recommended.
2. **Discoverability first (R1 + R4 + R5).** Solves "where did it go?" today. Defers durable dashboard and SHIPPED state.
3. **Dashboard first (R2 + R3).** Solves "what is WIP / what shipped?" today. Defers revision-log pinning and anchor citation.

Recommendation: **all-in**. R1 and R4 are cheap and load-bearing for R2 to be trustworthy (a dashboard that points to unsearchable revision logs and uncited anchors is half-useful). R3 without R5 leaves the lifecycle underspecified to a new reader. The coherence is in the bundle; slicing it trades coherence for speed.

**One-line takeaway.** The harness routes work correctly today; this proposal makes it *legible* — predictable locations, durable views, and anchors a human can click.
