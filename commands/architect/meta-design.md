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

Mirror of `commands/design.md` Phase 3b for the meta-layer. After flipping MADRs to ACCEPTED, walk `docs/_meta/design.md` and append every **time-boxed** portfolio item to `docs/_meta/deferred.md`. Permanent rejections stay in `_meta/design.md` only.

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

If `docs/_meta/deferred.md` does not yet exist, create it with the canonical header from `commands/design.md` Phase 3b before adding the first row. Read the file before appending to find the current max ID.

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
