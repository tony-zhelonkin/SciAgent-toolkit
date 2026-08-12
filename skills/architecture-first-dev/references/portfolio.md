# Portfolio / Meta Specifications

Load this blob for portfolio-level routing. Each marked source section preserves the complete command specification and its original argument contract.

<!-- BEGIN SOURCE: commands/architect/meta-map.md -->
# Meta-map

Produce the portfolio-level cartography artifact `docs/_meta/map.md` from the per-feature maps + reviews + designs that already exist.

This is the entry point of the meta-layer. Use it when ≥2 features are in flight and you want one place that names their shared touch-points, convergent concerns, and inter-feature dependencies.

## Phase 0: Parse arguments

- `$ARGUMENTS` — optional comma-separated list of feature slugs (e.g. `normalisation,umap,theme-bundles`)
  - Empty / absent → auto-detect: every subdirectory of `docs/` that contains a `map.md`, except `docs/_meta/`
  - Explicit list → use exactly those slugs

Announce the resolution before proceeding:

```
Resolved meta-map scope:
  features = [{slug}, {slug}, ...]
  output    = docs/_meta/map.md
```

If a slug was passed that has no `docs/{slug}/map.md`, warn and ask whether to drop it or run `/map {slug}` first. WAIT.

## Phase 1: Verify prerequisites

The meta-layer is only useful when ≥2 features have been mapped. Count features in scope that have `docs/{slug}/map.md`.

- **<2 features have a map.md** → reject:

  ```
  Meta-layer requires ≥2 mapped features. Currently: {N}.

  Run /map <slug> on at least one more feature, then re-run /meta-map.
  ```

  Stop.

- **≥2** → continue.

## Phase 2: Check existing state

Always check whether `docs/_meta/map.md` already exists — this phase always runs, just with different outcomes:

- **Exists** → ask the user to choose:

  ```
  docs/_meta/map.md already exists ({date}, covers features [...]).

  Choose:
    1. Use existing — skip dispatch, surface what's there
    2. Regenerate — overwrite from current per-feature artifacts
    3. Append updates — re-dispatch with instruction to update only changed sections
  ```

  WAIT for choice. On (1), skip Phase 3 and go straight to Phase 4.

- **Does not exist** → state `No existing _meta/map.md — proceeding to meta-architect dispatch.` and continue. (Do NOT describe this phase as "skipped" — the check ran; there was simply nothing to prompt about.)

## Phase 3: Dispatch meta-architect

Use the Agent tool with `subagent_type: "meta-architect"`.

Prompt template:

```
Produce docs/_meta/map.md.

Features in scope: {comma-separated slugs}

For each feature read in this order:
  1. docs/{slug}/map.md (REQUIRED)
  2. docs/{slug}/synthesis.md (preferred if present)
  3. docs/{slug}/review/*.md (read all if no synthesis)
  4. docs/{slug}/design/*.md (if present, especially 03-decisions.md)
  5. docs/{slug}/audit.md or scope.md (if present)

Also read docs/_meta/design.md and docs/_meta/plan.md if they already exist (do not contradict accepted MADRs without flagging).

Output: docs/_meta/map.md following the _meta/map.md template in your system prompt. Cover:
  - Features in scope (status table)
  - Shared touch-points (files/modules touched by >1 feature, with file:line nature)
  - Convergent concerns (the cross-feature theme; ≥2 features must surface each)
  - Inter-feature dependencies (explicit blocking relationships)
  - Open cross-feature questions (OQ-M-N)
  - Reader map

Facts and citations only. Every claim cites a path inside docs/. Do not read source code.
```

The meta-architect agent has tools `Read, Grep, Glob, Write` — Grep and Glob are restricted to `docs/` (do not let it loose on the codebase). It will produce exactly one file.

## Phase 4: Surface the meta-map

After meta-architect completes (or after the user chose "use existing" in Phase 2):

1. Read `docs/_meta/map.md`
2. Print to user:

```
Meta-map ready: docs/_meta/map.md ({N} lines)

Features in scope:    {count}
Shared touch-points:  {count}
Convergent concerns:  {count}
Inter-feature deps:   {count}
Open questions (OQ-M): {count}
```

3. Quote the Convergent concerns headlines (titles only, not full entries).
4. Quote any `OQ-M-N` titles.
5. Recommend next step:

```
Recommended next steps:
  /meta-design                  — draft inter-feature ADRs (MADRs) from the convergent concerns
  /review {slug} --as <list>    — if a specific feature still needs review before meta-design
  Stop                          — if this was a snapshot only
```

## Rules

1. **`_meta/map.md` is produced only by `meta-architect`.** The main agent never writes it directly.
2. **Meta-architect reads artifacts, not code.** If meta-architect's output cites a path outside `docs/`, ask it to revise.
3. **The directory `docs/_meta/` is created if missing.** No manual setup required.
4. **Do not auto-progress to `/meta-design`.** Even if the meta-map looks complete, stop and let the user choose.
5. **`OQ-M-N` are distinct from per-feature `OQ-N`.** If meta-architect reuses a per-feature numbering scheme, ask it to renumber.
<!-- END SOURCE: commands/architect/meta-map.md -->

<!-- BEGIN SOURCE: commands/architect/meta-design.md -->
# Meta-design

Draft `docs/_meta/design.md` — the set of inter-feature ADRs (MADRs) that the convergent concerns surfaced in `_meta/map.md` resolve. These MADRs are inherited by every per-feature `/design` and constrain its local ADRs.

## Phase 0: Parse arguments

- `$ARGUMENTS` — optional comma-separated list of feature slugs
  - Empty → use the same scope as the most recent `_meta/map.md` (read its `features:` frontmatter)
  - Explicit list → use those slugs (must be a subset of `_meta/map.md`'s scope, or the user is broadening scope and should re-run `/meta-map` first)

Announce the resolution:

```
Resolved meta-design scope:
  features = [{slug}, {slug}, ...]
  output    = docs/_meta/design.md
```

## Phase 1: Verify prerequisites

1. **`docs/_meta/map.md` MUST exist.** If absent:

   ```
   /meta-design requires docs/_meta/map.md.

   Run /meta-map first.
   ```

   Stop.

2. **Per-feature input check.** For every feature in scope, prefer `synthesis.md` ≥ `review/*.md` ≥ `design/*.md` ≥ `map.md` alone. If a feature has only `map.md`, warn:

   ```
   Warning: feature {slug} has no review or design — meta-design will have thin input for it.
   Continue anyway? (y/n)
   ```

   WAIT.

## Phase 2: Check existing state

- **`docs/_meta/design.md` exists** → ask:

  ```
  docs/_meta/design.md already exists ({date}, status: {STATUS}, MADRs: {count}).

  Choose:
    1. Use existing — proceed to Phase 7 (downstream impact)
    2. Regenerate — overwrite from current artifacts
    3. Iterate — re-dispatch with redirects (you'll provide them)
  ```

  WAIT.

- **Does not exist** → state `No existing _meta/design.md — proceeding to meta-architect dispatch.` and continue.

## Phase 3: Dispatch meta-architect

Use the Agent tool with `subagent_type: "meta-architect"`.

Prompt template:

```
Produce docs/_meta/design.md.

Features in scope: {comma-separated slugs}

Read in this order:
  1. docs/_meta/map.md (REQUIRED — primary input)
  2. docs/_meta/design.md (if exists — prior MADRs you may iterate on but should not silently contradict)
  3. docs/_meta/plan.md (if exists — note dependencies you must respect)
  4. For each feature:
     - docs/{slug}/synthesis.md (preferred)
     - docs/{slug}/review/*.md (if no synthesis)
     - docs/{slug}/design/03-decisions.md (if exists — pay close attention to ADRs and "Deliberate departures")
     - docs/{slug}/design/01-architecture.md and 02-behavior.md (if exists)

Produce MADR-NNN entries covering EVERY convergent concern listed in _meta/map.md. Each MADR follows the template in your system prompt:
  - Status, Features affected, Source citation
  - Context, Options considered, Decision, Rationale, Trade-offs
  - Consequences per feature (for every affected feature, name what their per-feature design must do)

Also include:
  - Open inter-feature questions (carry from _meta/map.md, mark resolution state)
  - Deliberate meta-level rejections (concerns where we deliberately do NOT pin a portfolio decision)
  - Reader map

{If iterating: include user redirects here}

Status of all newly drafted MADRs: PROPOSED. The user will approve at Gate 1, after which the dispatching command updates status to ACCEPTED.
```

## Phase 4: Gate 1 — scope/MADR-set approval

Read `docs/_meta/design.md` and present:

```
Meta-design drafted at docs/_meta/design.md.

MADR headlines:
  MADR-001 — {title} (affects: {f1, f2})
  MADR-002 — {title} (affects: {f1, f3})
  ...

Open inter-feature questions: {count}
  OQ-M-1 — {one-line}
  ...

Deliberate meta-level rejections: {count}
  - {topic} — {reason}

Please review and:
  1. Approve all — I'll set status to ACCEPTED, append time-boxed items to _meta/deferred.md, and surface downstream impact (Phase 6/6.5/7)
  2. Approve some, revise others — name the MADRs to revise and the redirect for each
  3. Reject specific MADRs — name them and the reason (will move to "Deliberate meta-level rejections")
  4. Add a new MADR — name the concern; I'll re-dispatch
  5. Questions — ask about specific MADRs
```

WAIT for the user.

## Phase 5: Iteration loop

If the user requested changes (options 2, 3, or 4):

1. Capture the redirects exactly.
2. Re-dispatch meta-architect (Phase 3) with the redirects appended to the prompt.
3. Loop back to Phase 4. Repeat until the user picks option 1.

## Phase 6: Mark as approved

On approval, set the frontmatter `status:` of `docs/_meta/design.md` to `APPROVED` and update every MADR's `Status:` line from `PROPOSED` to `ACCEPTED` (unless the user explicitly rejected one — those become `Status: REJECTED` and stay listed for traceability).

This is a small edit; the main agent does it directly (do not re-dispatch the agent for status flips).

## Phase 6.5: Append portfolio-level deferrals to `_meta/deferred.md`

Mirror of `skills/architecture-first-dev/references/design.md` Phase 3b for the meta-layer. After flipping MADRs to ACCEPTED, walk `docs/_meta/design.md` and append every **time-boxed** portfolio item to `docs/_meta/deferred.md`. Permanent rejections stay in `_meta/design.md` only.

What counts as time-boxed (must be appended):

1. **Open inter-feature questions** — every `OQ-M-N` whose resolution line says `still-open — deferred to ...`, `parked-in-deferred`, or names a future round / future ticket. Items resolved by a MADR (`resolved-by-MADR-N`) are NOT appended.
2. **Deliberate meta-level rejections** — every entry that defers to a future round, "after Round 2 implement", "monitor and revisit", or names a Round-N ticket. Permanent rejections (no time horizon, no monitoring trigger) are NOT appended.
3. **MADR Consequences referencing a future round** — e.g., a Consequence line that says "no further keys in Round 2", "revisit when X ships", or "Round-3 follow-up". Append one row per distinct future-round commitment.
4. **MADR-internal "Time-boxed sub-items" subsection** — if the meta-architect populated this subsection in the `Deliberate meta-level rejections` block (per `agents/meta-architect.md`'s output template), each line is a deferred row.

Schema (full row, append-only — never overwrite existing rows):

| Field | Value for meta-level entries |
|-------|------------------------------|
| ID | auto-increment from max existing (`D-001` if file is empty / new) |
| Source feature | the literal string `_meta` |
| Source doc | `_meta/design.md § <section>` (e.g., `§ OQ-M-8`, `§ Deliberate meta-level rejections`, `§ MADR-006 Rationale`) |
| Description | one-line |
| Reason deferred | one-line |
| Cost | `S` / `M` / `L` heuristic |
| Depends on | explicit dependency or `—` |
| Status | `DEFERRED` |

If `docs/_meta/deferred.md` does not yet exist, create it with the canonical header from `skills/architecture-first-dev/references/design.md` Phase 3b before adding the first row. Read the file before appending to find the current max ID.

This is a small edit; the main agent does it directly (do not re-dispatch the agent).

After appending, surface to the user:

```
Appended {N} rows to docs/_meta/deferred.md (D-{first}..D-{last}):
  - D-{n} — {description} (source: {section})
  - ...

These will appear as Parking items when /meta-plan runs.
```

## Phase 7: Surface downstream impact

For each feature in scope that already has `docs/{slug}/design/`, check whether any MADR's "Consequences per feature" line conflicts with an existing per-feature ADR.

Heuristic (no need to be perfect — surface for human review):
- If a MADR consequence says "f1 must do X" and `docs/f1/design/03-decisions.md` already has an ADR whose Decision contradicts X, flag it.
- If a MADR is brand new and the per-feature design doesn't address it at all, flag it as "may need re-design".

Print:

```
Meta-design approved: docs/_meta/design.md ({M} MADRs ACCEPTED, {R} REJECTED)
Deferred items appended: docs/_meta/deferred.md (D-{first}..D-{last}, {N} new rows)

Downstream impact on existing per-feature designs (sorted by descending conflict count):
  - {feature}: {N} potential conflicts
      - MADR-001 vs docs/{feature}/design/03-decisions.md § ADR-005
      - ...
  - {feature}: {N} potential conflicts
      - ...
  - {feature}: 0 conflicts (edit-and-cite pass — only needs MADR citations added)
  - {feature}: no existing design — next /design will inherit MADRs cleanly

Recommended next steps (primary path = parallel propagation; fallback = conflict-density-ordered serial):

  /meta-apply                                  — PRIMARY. Propagate MADRs to every in-flight
                                                 design in one parallel batch. Dispatches N
                                                 feature-revisers simultaneously; one compound
                                                 gate on the diff; architect gate then fans out in
                                                 parallel per feature. Features flagged CONFLICT
                                                 in this Phase 7 report will be emitted as
                                                 MANUAL_REDESIGN_NEEDED by their feature-reviser
                                                 and cleanly skipped. Typical saving vs serial
                                                 /design: ~3× tokens, ~5× wall-clock, 1 compound
                                                 gate vs ~10 per-feature gates.

  /design <feature>                             — FALLBACK. Use for:
                                                   • features flagged CONFLICT in this report when
                                                     the conflict is structural and the sub-agent
                                                     correctly refused to speculate
                                                   • features with no `design/` yet (first-draft)
                                                   • features flagged NEEDS_ITERATION by architect
                                                     after /meta-apply (loop until READY)
                                                   • any feature you want a finer-grained,
                                                     hand-paced gate on

  /meta-plan                                   — once every per-feature design is READY, sequence
                                                 the portfolio.
```

**When to prefer `/meta-apply` over serial `/design`**: if Phase 7 reports mostly NEEDS_UPDATE
lines with ≤1 structural CONFLICT per feature, `/meta-apply` is strictly cheaper. If it reports
structural CONFLICTs that need re-architecting (e.g., a feature's core ADR must be re-thought
against a MADR, not just annotated), serial `/design` gives a finer-grained gate per feature. The
`/meta-apply` sub-agents auto-emit `MANUAL_REDESIGN_NEEDED` for such cases, so running
`/meta-apply` first is safe — it refuses to speculatively revise features it can't handle, and
produces no writes for those features. You lose nothing by trying it first.

**Cheap-first fallback ordering** (only when `/meta-apply` is inapplicable — e.g., no features
have `design/` yet, or every feature is flagged CONFLICT): order serial `/design` re-runs by
descending MADR-vs-existing-ADR conflict count from this Phase 7 report. Most-conflicted feature
first; zero-conflict features last. A MADR-vs-ADR conflict is a stop-gate — `/design` Phase 1
STOPs and asks the user to revise the MADR, override locally, or cancel. Surfacing the conflict
on the **most-conflicted design first** means any MADR revision happens *before* less-conflicted
features re-design and inherit the (now-revised) MADR. A heavy design with zero conflicts is
just an edit-and-cite pass; it produces no signal that could change a MADR and goes last
regardless of size.

**Computing the order**: count conflicts per feature directly from this Phase 7 report. Do NOT
order by descending `design/` size — that's only the fallback-to-fallback proxy when conflict
counts are unknown. A feature with `design/03-decisions.md` containing 12 ADRs that all align
with the MADRs has 0 conflicts and goes last; a feature with only `design/README.md` whose 3
planned ADR headlines all conflict has 3 conflicts and goes first.

## Rules

1. **Meta-design output is ADRs only.** No code-level specifics, no API signatures, no file:line edits.
2. **Every MADR names the features it affects.** No portfolio decision applies in a vacuum.
3. **Meta-design never overrides a feature-local decision in the feature's own doc.** It sets a constraint that the next per-feature `/design` must honour or explicitly flag as a meta-conflict.
4. **MADRs are PROPOSED until Gate 1.** Only the dispatching command may flip them to ACCEPTED.
5. **Convergent concerns map to MADRs or to "Deliberate meta-level rejections" (or both as a hybrid when honest)** — silence on a concern is a defect, but a concern that resolves a name (MADR) while deferring a semantic sub-question (rejection sub-item appended to deferred.md) is acceptable. Re-dispatch only if a concern is unaccounted for in either section.
6. **Time-boxed deferrals are appended to `_meta/deferred.md` by Phase 6.5.** Permanent rejections are not. The dispatching command does this append directly; the meta-architect agent only marks items for the append by listing them in the `Deliberate meta-level rejections § Time-boxed sub-items` subsection (see `agents/meta-architect.md`).
6. **Do not auto-progress to `/meta-plan`.** Surface the recommendation; let the user choose.
<!-- END SOURCE: commands/architect/meta-design.md -->

<!-- BEGIN SOURCE: commands/architect/meta-apply.md -->
# Meta-apply

Propagate every accepted MADR in `docs/_meta/design.md` to every in-flight per-feature design in a single parallel batch. Dispatches N `feature-reviser` subagents in parallel, surfaces one compound diff gate, then dispatches N `architect` subagents in parallel for the architecture-gate step.

Use after `/meta-design` Gate 1 approval, in place of N serial `/design` re-runs, whenever the propagation is mechanical (apply MADR consequences to existing per-feature design/*.md). Features flagged `MANUAL_REDESIGN_NEEDED` by a sub-agent are skipped cleanly — the user runs `/design <slug>` for those.

## Phase 0: Parse arguments

- `$ARGUMENTS` — optional comma-separated list of feature slugs.
  - Empty → auto-select every feature in `docs/_meta/design.md`'s frontmatter `features:` list that ALSO has a `docs/{slug}/design/` directory.
  - Explicit list → those slugs. Every slug MUST appear in `_meta/design.md`'s `features:` frontmatter (else reject: "slug {x} not in `_meta/design.md` scope; re-run /meta-design with updated scope or fix the slug").

Announce resolution:

```
Resolved meta-apply scope:
  features with design/       = [...]       ← feature-revisers will be dispatched for these
  features without design/    = [...]       ← skipped; run /design <slug> for first-draft
  output                      = edits to docs/{slug}/design/*.md (each in-scope feature)
                              + append to docs/_meta/deferred.md
```

If zero features have `design/`, reject: "No features in scope have a `design/` directory. Use `/design <slug>` for each — first-drafts are outside `/meta-apply`'s scope."

## Phase 1: Verify prerequisites

1. **`docs/_meta/design.md` MUST exist with `status: APPROVED`** in frontmatter. If missing or `status: DRAFT`:

   ```
   /meta-apply requires docs/_meta/design.md with status: APPROVED.

   Found: {absent | status: DRAFT | status: OTHER}

   Run /meta-design and approve at Gate 1 first.
   ```

   Stop.

2. **At least one feature in scope must have a `design/` directory** (enforced in Phase 0 already; restate if list is empty after filtering).

3. **No concurrent run in progress.** If any `docs/{slug}/design/*.md` in scope has been modified in the last 60 seconds (mtime-based sanity check), warn:

   ```
   Warning: docs/{slug}/design/*.md was modified recently. Another /meta-apply or /design may be running. Continue anyway? (y/n)
   ```

   WAIT. Otherwise continue.

## Phase 2: Allocate deferred-ID blocks (serialised, pre-dispatch)

**This phase runs in the main agent, BEFORE parallel dispatch.** Parallel appends to `_meta/deferred.md` would race; pre-allocation prevents that.

1. Read `docs/_meta/deferred.md` and find the current maximum ID `D-NNN`. If the file does not exist yet, the watermark is `D-000` (next available is `D-001`).
2. Assign each feature-reviser a contiguous 20-ID block starting above the watermark:
   - feature A (first in scope) → `D-{watermark+1}..D-{watermark+20}`
   - feature B → `D-{watermark+21}..D-{watermark+40}`
   - feature C → `D-{watermark+41}..D-{watermark+60}`
   - ...
3. Announce allocations before dispatch:

   ```
   Deferred-ID block allocations (exclusive per sub-agent):
     normalisation        → D-013..D-032
     umap                 → D-033..D-052
     color-encoding       → D-053..D-072
     theme-bundles        → D-073..D-092
     biological-workflow  → D-093..D-112

   Watermark after this run (if all blocks fill): D-112. Actual new rows will be ≤100.
   ```

The `20` block size is calibrated so a well-scoped sub-agent never exhausts it. A sub-agent reporting its block exhausted is a signal of overreach — the sub-agent will emit `MANUAL_REDESIGN_NEEDED` rather than steal from a sibling's block.

## Phase 3: Parallel dispatch of `feature-reviser` sub-agents

**CRITICAL:** Use a SINGLE assistant message with MULTIPLE `Agent` tool calls — this is what makes them run in parallel. Same pattern as `/review --as all` (see `skills/architecture-first-dev/references/review.md § Phase 3`). Sequential dispatch violates the design and defeats the cost savings.

For each feature in scope, derive these fields from `_meta/design.md`:

- **`MADRs to apply`**: every MADR section whose "Consequences per feature" block contains this slug on any bullet (grep-style scan; cite the MADR number + one-line consequence quote).
- **`Reader-map sections`**: entries under `_meta/design.md § Reader map` keyed by this slug.

Then, for each feature, issue one Agent call with `subagent_type: "feature-reviser"` and this exact prompt body (parameterised per feature):

```
Apply docs/_meta/design.md MADRs to feature `{slug}`.

Your assigned deferred-ID block: D-{first}..D-{last} (use ≤20 IDs; do not exceed).

MADRs to apply (and the `Consequences per feature: {slug}` lines that name this slug):
  - MADR-{N1} (§{section-header-1}) — consequence for {slug}: "{one-line quote}"
  - MADR-{N2} (§{section-header-2}) — consequence for {slug}: "{one-line quote}"
  - ...

Reader-map sections to re-read (from _meta/design.md § Reader map):
  - {path § ADR-N}
  - ...

Read in this order:
  1. docs/_meta/design.md (binding; all listed MADRs)
  2. docs/_meta/map.md (inter-feature dependencies — cite, never modify)
  3. docs/_meta/deferred.md (current state; your assigned ID block is exclusive to you)
  4. docs/{slug}/map.md, synthesis.md (if present), review/*.md, design/*.md

Produce:
  - Edits under docs/{slug}/design/*.md for every MADR consequence you can apply mechanically
  - Row appends to docs/_meta/deferred.md using your ID block
  - Diff summary as your final message (exact format in your agent definition)

If ANY listed MADR cannot be applied mechanically, emit `status: MANUAL_REDESIGN_NEEDED` in your
diff summary with a one-paragraph explanation per conflict. Do NOT write speculative revisions.
```

Dispatch all N calls in one message. Do not wait for one before sending the next.

## Phase 4: Collect diffs + surface compound report

After all parallel feature-revisers return, read each agent's final diff summary. Then, for each feature, read the post-dispatch state of `docs/{slug}/design/*.md` to confirm the reported edits actually landed (trust but verify — an agent's self-report describes intent, not necessarily effect).

Print the compound report. Table first, per-feature highlights second, MANUAL flags quoted prominently at the end:

```
Meta-apply dispatched {N} feature-revisers. Results:

┌──────────────────────┬─────────┬────────────────┬──────────────┬──────────────┬──────────────┬──────────┐
│ Feature              │ Status  │ Files touched  │ ±Lines       │ Revisions    │ Deferrals    │ Conflict │
├──────────────────────┼─────────┼────────────────┼──────────────┼──────────────┼──────────────┼──────────┤
│ normalisation        │ OK      │ 03, 05         │ +89 / -12    │ 2 ADRs       │ D-013..014   │ —        │
│ umap                 │ OK      │ README, 01, 03 │ +45 / -31    │ 3 ADRs       │ —            │ —        │
│ color-encoding       │ MANUAL  │ (none written) │ —            │ —            │ —            │ ADR-005  │
│ theme-bundles        │ OK      │ README, 03     │ +110 / -46   │ first-draft  │ D-073..075   │ —        │
│ biological-workflow  │ OK      │ 03, 04         │ +120 / -18   │ 3 ADRs       │ D-093..099   │ —        │
└──────────────────────┴─────────┴────────────────┴──────────────┴──────────────┴──────────────┴──────────┘

Per-feature highlights:
  - normalisation — ADR-006 Status OPEN → ACCEPTED (Option B), authorised by MADR-006;
                    ADR-013 opacity-for-degenerate re-homed to marker overlay per MADR-002.
  - umap — ADR-001/ADR-008 gene-set fields renamed per MADR-004; MADR-003 cited in ADR-003.
  - color-encoding — MANUAL_REDESIGN_NEEDED: feature has only design/README.md (no 03-decisions.md
                     yet); planned ADR-005 (opacity = padj-percentile) conflicts with MADR-002 at
                     portfolio level. First-draft + structural rework needed — not mechanical
                     edit. Run /design color-encoding.
  - theme-bundles — runtime column resolver dropped per MADR-001; AC revised for MADR-002 priority
                    chain; bundle-dim composition rule added.
  - biological-workflow — opacity priority chain added; ADR-010 hue revised for MADR-006;
                          ENTITY_PROFILES consumption captured in ADR-011.

MANUAL_REDESIGN_NEEDED flags (run /design <slug> for each — /meta-apply will not touch these):
  - color-encoding: {one paragraph explanation from sub-agent, verbatim}
```

## Phase 5: Compound approval gate (Gate 1 of /meta-apply)

Single gate. Five options:

```
Please review the compound diff above and pick:

  1. Approve all OK-status features
     → I'll flip DRAFT → APPROVED on each approved feature's design/*.md frontmatter
     → Dispatch architect gate in parallel on every approved feature's design/*.md
     → MANUAL-flagged features are left untouched; run /design <slug> manually for each

  2. Approve some, reject others
     → Name the features to reject; I'll revert those edits (git restore on
        docs/{slug}/design/*.md and strip the feature's deferred rows from _meta/deferred.md)
     → Remaining OK features proceed to architect gate

  3. Redo a specific feature
     → Name the feature + a redirect; I'll re-dispatch its feature-reviser
       (re-uses the same deferred-ID block — no re-allocation)

  4. Review a specific feature's diff
     → Name the feature; I'll show file-by-file diff

  5. Questions
```

WAIT for the user.

## Phase 6: On approval — parallel architect dispatch

For every feature that survived the compound gate with `Status: OK`:

1. Bump `status: DRAFT → APPROVED` in the frontmatter of every `docs/{slug}/design/*.md` file (main agent does this directly — trivial, not worth re-dispatching).
2. Dispatch `architect` subagent via Agent tool. **Use a SINGLE assistant message with MULTIPLE Agent tool calls** (same parallel-dispatch pattern as Phase 3 and `/review`). One architect per feature.

Architect prompt template (identical to `skills/architecture-first-dev/references/design.md § Phase 5`):

```
Review docs/{slug}/design/*.md against docs/{slug}/map.md.

Check internal consistency, completeness, pattern alignment with map.md, and design principles.
Pay special attention to the "Inherited meta-design constraints" block — the feature-reviser
just added or extended it; verify every MADR citation maps to a coherent ADR revision.

Write docs/{slug}/design/review.md with Verdict: READY FOR IMPLEMENTATION | NEEDS ITERATION | NEEDS DISCUSSION.
```

Dispatch all N calls in one message.

## Phase 7: Collect verdicts + emit next-step recommendation

After all parallel architects return, read each `docs/{slug}/design/review.md` and parse the Verdict line. Surface per-feature outcomes:

```
Meta-apply complete — {N} features revised, architect gate run in parallel:

  ✅ normalisation            READY FOR IMPLEMENTATION
  ✅ umap                     READY FOR IMPLEMENTATION
  ⚠️ biological-workflow      NEEDS ITERATION (see docs/biological-workflow/design/review.md)
  ✅ theme-bundles            READY FOR IMPLEMENTATION
  ⏭️ color-encoding           (skipped — MANUAL_REDESIGN_NEEDED at Phase 4)

Deferred rows appended: D-013..D-099 (N new rows total; _meta/deferred.md watermark now D-099).

Recommended next steps:
  /design color-encoding              — for the MANUAL-flagged feature (first-draft + re-architecting)
  /design biological-workflow --iterate — architect flagged NEEDS_ITERATION; loop to READY
  /meta-plan                          — once every feature has Verdict READY, sequence the portfolio
```

For each feature with a NEEDS ITERATION or NEEDS DISCUSSION verdict, **roll back the status bump**: flip `status: APPROVED → DRAFT` on that feature's `design/*.md` frontmatter. The verdict means the architect isn't satisfied; `APPROVED` would be a lie. A subsequent `/design <slug> --iterate` run will re-advance it when Verdict = READY.

For READY features, leave `status: APPROVED` in place — `/plan` requires it.

## Phase 8: Status-field housekeeping (trivial, main-agent)

After Phase 7, the state of every feature's frontmatter is:

| Verdict | design/*.md `status:` |
|---|---|
| READY FOR IMPLEMENTATION | `APPROVED` |
| NEEDS ITERATION | `DRAFT` |
| NEEDS DISCUSSION | `DRAFT` |
| (skipped, MANUAL) | unchanged from pre-run (whatever it was) |

Confirm in one line:

```
Frontmatter status set: {READY count} APPROVED, {NEEDS count} DRAFT, {MANUAL count} unchanged.
```

## Rules

1. **Parallel dispatch is mandatory.** Phase 3 (feature-revisers) and Phase 6 (architects) each use a single assistant message with multiple Agent tool calls. Sequential dispatch defeats the entire point of this command.
2. **Compound gate ONCE, architect gate PER FEATURE.** Gate 1 approves the batch of edits. Then architect runs in parallel on each approved feature and returns per-feature verdicts — no second compound gate on architect output.
3. **MANUAL features are skipped, not failed.** A sub-agent emitting `MANUAL_REDESIGN_NEEDED` is behaving correctly — it refused to speculate. The user's cue is to run `/design <slug>` for that specific feature.
4. **Deferred-ID blocks are pre-allocated and exclusive.** Sub-agents never allocate their own IDs. A sub-agent reporting its block exhausted is overreach → that feature is MANUAL.
5. **Sub-agents never modify `_meta/design.md` or `_meta/map.md`.** Those belong to `meta-architect`. `_meta/deferred.md` is the only shared meta-layer file sub-agents may append to.
6. **Do not auto-progress to `/meta-plan`.** Surface the recommendation; the user starts when ready.
7. **Do not self-loop on NEEDS_ITERATION.** Surface the verdict; the user runs `/design <slug> --iterate` (or re-runs `/meta-apply <slug>` with a redirect) manually. Self-looping would hide architect signal behind another compound turn.
8. **Large diffs (>200 lines touched in a single feature) are a smell, not a defect.** Surface them; suggest the user drill in via option 4 or, if the edits feel aggressive, re-run `/design <slug>` manually for a careful hand pass.
<!-- END SOURCE: commands/architect/meta-apply.md -->

<!-- BEGIN SOURCE: commands/architect/meta-plan.md -->
# Meta-plan

Sequence per-feature implementation phases across the portfolio so dependencies land before their consumers, and so simultaneously active phases don't collide on the same files.

## Phase 0: Parse arguments

- `$ARGUMENTS` — optional comma-separated list of feature slugs
  - Empty → use the same scope as `_meta/design.md` (read its `features:` frontmatter)
  - Explicit list → must be a subset of `_meta/design.md`'s scope

Announce:

```
Resolved meta-plan scope:
  features = [{slug}, {slug}, ...]
  output    = docs/_meta/plan.md
```

## Phase 1: Verify prerequisites

1. **`docs/_meta/design.md` MUST exist with `status: APPROVED`.** If absent or DRAFT:

   ```
   /meta-plan requires docs/_meta/design.md with status APPROVED.

   Current: {missing | DRAFT}
   Run /meta-design (and approve at Gate 1) first.
   ```

   Stop.

2. **Per-feature plan input.** For each feature in scope, check whether `docs/{slug}/plan/phase-NN.md` files exist. If a feature has no plan files:

   ```
   Warning: feature {slug} has no plan/ directory — its phases cannot be sequenced.

   Choose:
     1. Continue — meta-plan will note "{slug}: no phases yet" and skip sequencing for it
     2. Stop — run /plan {slug} first, then re-invoke /meta-plan
   ```

   WAIT.

3. **At least one feature must have plan files.** If every feature in scope lacks plans, reject — there is nothing to sequence.

## Phase 2: Check existing state

- **`docs/_meta/plan.md` exists** → ask:

  ```
  docs/_meta/plan.md already exists ({date}, covers {N} portfolio phases).

  Choose:
    1. Use existing
    2. Regenerate
    3. Iterate with redirects
  ```

  WAIT.

- **Does not exist** → state `No existing _meta/plan.md — proceeding to meta-architect dispatch.` and continue.

## Phase 3: Dispatch meta-architect

Use the Agent tool with `subagent_type: "meta-architect"`.

Prompt template:

```
Produce docs/_meta/plan.md.

Features in scope: {comma-separated slugs}

Read in this order:
  1. docs/_meta/design.md (REQUIRED — MADR consequences drive ordering)
  2. docs/_meta/map.md (for inter-feature dependencies)
  3. For each feature with a plan/:
     - docs/{slug}/plan/README.md (phase summary table)
     - docs/{slug}/plan/phase-NN.md (per-phase Files-to-Create / Files-to-Modify)
     - docs/{slug}/design/03-decisions.md (for ADR cross-references)
  4. For each feature WITHOUT a plan/:
     - Note "{slug}: no phases yet" in the output; do not invent phases

Produce docs/_meta/plan.md following the _meta/plan.md template in your system prompt. It MUST contain:
  - Dependency graph (Mermaid) — feature-phases as nodes, blocking relationships as edges
  - Portfolio-phase order — table assigning every feature-phase to a portfolio phase P0..PN
  - Shared-file collision check — for each portfolio phase, list files modified by every active phase and flag overlaps
  - Parking items — items from _meta/deferred.md flagged for a future wave
  - Reader map

A phase that depends on another (per MADR consequence or _meta/map.md inter-feature dependency) MUST land in a later portfolio phase than its dependency.

If two simultaneously active feature-phases touch overlapping line ranges in the same file, flag it as "overlap — split or sequence" rather than silently sequencing.

Read but do not modify per-feature plan.md or design.md files.
```

## Phase 4: Gate — present sequencing for approval

Read `docs/_meta/plan.md` and present:

```
Meta-plan drafted at docs/_meta/plan.md.

Portfolio phase order:
  P0 — {features and their phases}    rationale: {one line}
  P1 — ...
  ...

Dependency graph:
  {render the Mermaid graph as ASCII or quote it raw}

Collision check:
  - {N} clean phases
  - {M} overlap warnings (these are blockers):
      P{N}: {file} — touched by {fA phase-X, fB phase-Y}

Parking items: {count from deferred.md activated for future waves}

Please review and:
  1. Approve — sequencing is correct
  2. Request resequencing — name the swap (e.g., "move biological-workflow phase-01 to P1")
  3. Resolve a collision — name the file and which feature should land first
  4. Questions
```

WAIT.

## Phase 5: Iteration

If the user requests changes, capture redirects and re-dispatch (Phase 3) with redirects appended. Loop until approved.

## Phase 6: Surface downstream

```
Meta-plan approved: docs/_meta/plan.md

Portfolio-phase assignments:
  - {feature}: phase-01 → P0, phase-02 → P2
  - {feature}: phase-01 → P2
  - {feature}: no plan/ — run /plan {slug} before participating

Recommended next steps:
  /plan {feature}                     — for any feature whose plan does not yet exist
  /implement {feature} 1              — start P0 work (per portfolio sequence)
```

## Rules

1. **Meta-plan is sequencing, not phase decomposition.** Per-phase decomposition is `/plan`'s job.
2. **No collision is silently sequenced.** A flagged overlap requires user resolution before any P{N} containing it goes live.
3. **Dependency edges come from `_meta/design.md` (MADR consequences) and `_meta/map.md` (inter-feature dependencies).** Do not invent dependencies; cite the source.
4. **A feature without a plan is noted, not faked.** Do not generate placeholder phases.
5. **Do not auto-progress to `/implement`.** Surface the recommendation; let the user start when ready.
<!-- END SOURCE: commands/architect/meta-plan.md -->
