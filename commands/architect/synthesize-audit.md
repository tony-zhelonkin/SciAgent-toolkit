---
status: probationary
prune_review_at: 2026-11-20
---

# Synthesize-audit

Integrate N audit slices (concerns) plus the `/components-extract` substrate into the FULL, renderable `components.json` — the step that turns the deterministic substrate into an architectural model. Also emits a `synthesis-vN.md` narrative and a decision menu.

This is **distinct from `/synthesize`.** Per-feature `/synthesize` collapses N reviews (perspectives on one map) into a consensus. `/synthesize-audit` operates on N slices (concerns traced through code) PLUS the extractor's physical substrate, and produces a different artifact: a logical-component model, classifications, judgment edges, prune verdicts, and an OSS-core boundary. Reviews give perspectives; slices give concern-depth; the substrate gives ground truth.

## The deterministic-sandwich position

```
/components-extract        →  substrate components.json   (DETERMINISTIC: physical files,
                                                            static edges, metrics, git_sha)
/audit-slice               →  slices/*.md                 (JUDGMENT: findings, connascence,
                                                            decoupling proposals)
/synthesize-audit (Ph 1-2) →  FULL components.json        (THIS STEP — merges substrate +
                                                            judgment edges into logical model)
rollup_logical_edges.py    →  components.json (updated)   (DETERMINISTIC: derives logical
                                                            edges from physical static graph)
/synthesize-audit (Ph 3)   →  validates manifest          (GATE: --strict orphan check)
/architecture-treemap      →  treemap.html                 (DETERMINISTIC render)
```

**Do not ask the agent for what a function can compute — the agent authors judgment; the tooling derives structure; the validator gates the result.**

The extractor emits ONLY statically-derivable facts. The synth agent adds everything that is NOT statically detectable: it groups physical files into `logical_components`, assigns classifications, and authors judgment edges (shared-state, background-knowledge, forward contracts). Structural import edges between logical components are derived DETERMINISTICALLY by `rollup_logical_edges.py` from the physical static graph — the agent does not author these. The substrate's static direct-call edges (`evidence_class: static`) and per-component metrics are carried through UNCHANGED — never invent or overwrite a static fact.

## Phase 0: Parse arguments

- `$ARGUMENTS[0]` (optional) — audit date dir. Defaults to the most recent `docs/_meta/architecture-audit/{date}/` containing a `slices/` subdir.
- `$ARGUMENTS` may also name an explicit substrate path (the extractor's output `components.json`).

If no audit dir resolves, tell the user to run `/audit-slice` first and stop.

## Phase 1: Locate inputs

1. **Slices:** `docs/_meta/architecture-audit/{date}/slices/*.md` — must have ≥1. If zero, stop and route to `/audit-slice`.
2. **Substrate:** the extractor's `components.json`. Look in the audit dir, then ask the user for the path. If absent, WARN:
   ```
   No /components-extract substrate found. The synthesis can proceed on slices alone,
   but the result will lack metrics badges, static edges, and the deterministic trust
   scaffold — every component will render "judgement only". Strongly recommend running
   /components-extract first. Proceed without substrate? (y/N)
   ```
   Default is to wait. The substrate is the trust scaffold; without it the treemap shows judgment with no ground truth behind it.

## Phase 1.5: Archive prior synthesis on re-run

If `docs/_meta/architecture-audit/{date}/synthesis-vN.md` already exists for the highest N, and a slice file is newer than it, this is a re-synthesis. Determine the next version: let `M = max({i : synthesis-vi.md exists})` (or `0`), produce `synthesis-v{M+1}.md`. Do NOT overwrite an existing `synthesis-vi.md` — versions are immutable; a new run adds a new version. State which version is being produced.

## Phase 2: Dispatch synth (audit variant)

Use the Agent tool with `subagent_type: "synth"`. The synth agent's default brief is per-feature consensus; the audit variant below overrides it.

Prompt template:
```
Produce an AUDIT synthesis, not a per-feature consensus. You are integrating N
concern-slices PLUS a deterministic substrate into a renderable model.

Read:
  - docs/_meta/architecture-audit/{date}/slices/*.md   (the concern slices)
  - {substrate components.json path}                    (the /components-extract output:
                                                          physical_components, static edges,
                                                          metrics, git_sha, extractor block)
  - {prior synthesis-v{M}.md, if a re-run — refine against it}

You read these artifacts + the substrate. You MAY spot-verify a slice's file:line
claim, but the slices already did the traversal — do not re-run the audit.

Produce TWO outputs:

1. docs/_meta/architecture-audit/{date}/components.json — the FULL manifest, validating
   against skills/architecture-treemap/components.schema.json. Specifically:
   - Carry the substrate's physical_components, static direct-call edges
     (evidence_class:"static"), per-component metrics, git_sha, snapshot_id, and
     extractor block through UNCHANGED.
   - Group physical files into logical_components. Each MUST have a classification
     (core | seam | removable), a description, physical_files refs, size_estimate_loc,
     and SHOULD carry the metrics rolled up from its physical files. Tag each with the
     audit_slices (slice ids) and smoke_findings that bear on it.
   - Set an optional `epistemic_source` on each classification and judgment edge:
     `measured` (deterministic), `metric-anchored` (judgment biased by metrics —
     typical for core), or `requires-your-intent` (extrinsic — typical for seam and
     removable). A `removable` verdict MUST trace to a slice/finding; never
     auto-assert removable from structure alone (the renderer marks an ungrounded
     removable "requires your judgment"). The renderer falls back to a sensible
     default when the field is absent, so it is optional but recommended.
   - Add ONLY the judgment edges the slices found: shared-state edges and
     background-knowledge edges (runtime contracts, invariants, forward contracts),
     each evidence_class:"audit-asserted" (renders dashed), with a `smell` and the
     slice/finding it came from. DO NOT author logical↔logical structural import
     edges — those are derived deterministically by rollup_logical_edges.py (Phase
     2.6) from the physical static graph. A logical pair may carry BOTH a derived
     solid edge (a real import was detected) AND an authored dashed edge (a runtime
     contract or design concern is also asserted) — that is meaningful, not a
     duplicate; do not try to reconcile them.
   - Write prune_candidates from the slices' decoupling proposals (R-NN): logical_component,
     bounded (true|false), loc_to_remove, unblocks_simplification_of, latent_bug_fixed.
   - Write core_boundary: rationale, core[], seam[], removable[], open_questions[].
   - Optionally write finding_index mapping each smoke/slice tag to a title+severity.

2. docs/_meta/architecture-audit/{date}/synthesis-v{N}.md — the narrative + a DECISION MENU.
   Sections: Headline; Substrate-vs-judgment summary (what the metrics flagged vs what the
   slices found); Logical-component model; Edge enumeration (static vs audit-asserted);
   Prune verdicts; Core/seam/removable boundary with rationale; then a `## Decision menu` —
   each decision as: options + recommended default + rationale + consequences. The decision
   menu is the artifact the user grills; make it the load-bearing section.

Rules:
  - Never invent or overwrite a static fact (metrics, static edges, git_sha). Those are
    the substrate's; carry them through.
  - Every classification, smell, audit-asserted edge, and prune verdict must trace to a
    slice finding (E-NN / R-NN) or a smoke finding. No synthesis-original architecture.
  - A logical_component with no slice and no finding gets classification by best judgment
    but note it as low-confidence in the narrative — the renderer shows it "judgement only".
  - The components.json MUST validate. After writing, the command validates it (Phase 3).
```

## Phase 2.6: Derive logical edges deterministically

After the synth agent writes the manifest (Phase 2) and before validation (Phase 3),
run the rollup script to derive logical↔logical structural import edges from the
physical static graph:

```
python3 skills/architecture-treemap/scripts/rollup_logical_edges.py \
    docs/_meta/architecture-audit/{date}/components.json
```

This script is IDEMPOTENT — it removes any previously derived edges (marked
`derived: true`) and re-derives from scratch, so re-runs accumulate nothing.
It prints one summary line to stderr: `rollup_logical_edges: derived N logical
edge(s) from M physical static import(s)`. Surface that line in the Phase-4 output.

The derived edges are tagged `"derived": true`, `evidence_class: "static"`,
`epistemic_source: "measured"`. They render as solid lines in the logical view —
the same visual weight as directly-detected physical imports, because that is what
they are (rolled up).

## Phase 3: Validate the manifest

Run the renderer's validator on the produced manifest:

```
python3 skills/architecture-treemap/scripts/validate_components.py \
    docs/_meta/architecture-audit/{date}/components.json --strict
```

If validation FAILS, surface the errors to the user and re-dispatch synth with the
specific schema/referential-integrity errors quoted. A non-leaf orphan (a logical
component with no logical-id edges after rollup) is an ERROR under --strict — it
means the rollup could not connect the node (no physical imports to roll up) and no
judgment edge was authored either; the agent must either author a judgment edge or
justify the node as a leaf (classification: removable, or explicit `"leaf": true`).
Do NOT route to `/architecture-treemap` on a failing manifest — the renderer refuses
to render invalid input anyway, and a half-valid manifest wastes the render pass.

## Phase 4: Surface and route

```
Audit synthesis complete (v{N}):
  docs/_meta/architecture-audit/{date}/components.json   (validates ✓)
  docs/_meta/architecture-audit/{date}/synthesis-v{N}.md

Headline:
  {quote the headline}

Boundary:
  core:      {core list}
  seam:      {seam list}
  removable: {removable list}

Prune candidates:
  {logical_component} — {loc_to_remove} LOC, bounded={bounded}

Decision menu: {count} decisions awaiting your call (see synthesis-v{N}.md § Decision menu).

Recommended next step:
  /architecture-treemap docs/_meta/architecture-audit/{date}/components.json
                                  — render the self-contained treemap.html from this model
```

## Robust authoring checklist

Hard-won lessons. The synth agent MUST run this checklist before the manifest is
considered done. These failures are silent — the manifest validates and renders
even when they are present — so the checklist is the only guard.

1. **No orphan logical nodes.** Structural import edges between logical components
   are derived DETERMINISTICALLY by `rollup_logical_edges.py` (Phase 2.6) — the
   agent is not responsible for hand-authoring them. The agent's responsibility is:
   (a) author judgment edges (shared-state, background-knowledge, forward contracts)
   where the slices found runtime coupling not visible in static imports; and
   (b) for any node that the rollup cannot connect (no physical imports, no judgment
   edge), justify it as a leaf via `classification: "removable"` or `"leaf": true`.
   The validator's no-non-leaf-orphan check under `--strict` is the enforcing gate —
   it fails loudly if any non-leaf logical node ends up edgeless after rollup.
2. **Ground every edge in a real source citation.** Each edge's `evidence`
   string must quote the import/call as written (file:line + the exact
   statement, aliases included). Set `evidence_class` honestly: `static` ONLY
   for a real detectable import/call (verify the line before claiming it);
   `audit-asserted` for runtime contracts and background-knowledge invariants.
3. **Cross-cutting nodes need representation.** Config, shared kernels, and util
   sinks attract high fan-in. Author their incoming edges (the renderer dims a
   high-fan-in / zero-fan-out / config node automatically so it does not
   hairball). Do not leave a shared dependency floating because "everyone uses
   it" — that is exactly the node whose edges must be drawn.
4. **The tool must not audit its own output.** Generated snapshots and artifacts
   (anything under `architecture-audit/`) are excluded by the extractor by
   default; never re-add them as components.
5. **Metrics are deterministic substrate, never judgment.** Read them to bias
   attention — high `refactor_pressure` / churn / fan-in = look there first —
   but never fabricate or override a metric. Metrics bias the audit queue; they
   do not deliver verdicts. The verdict is the classification, and it traces to
   a slice finding, not to a number.

## Rules

1. **Requires ≥1 slice.** Zero slices → route to `/audit-slice`. This step integrates judgment; it does not generate it.
2. **The substrate is the trust scaffold.** Proceed without it only on explicit user override, and warn that every component renders "judgement only".
3. **Never overwrite a static fact.** Metrics, static edges, git_sha come from the extractor and pass through unchanged. Synthesis adds the judgment layer; it does not edit the ground truth.
4. **Every judgment traces to a slice or smoke finding.** No synthesis-original classifications, smells, or edges.
5. **The manifest MUST validate before routing to render.** Validation gate is mandatory (Phase 3).
6. **Versions are immutable.** A re-run produces `synthesis-v{N+1}.md`; it never overwrites a prior version. (Same immutability discipline as the extractor's snapshots.)
7. **This command is probationary** (`prune_review_at: 2026-11-20`), paired with `/audit-slice`. The decision-menu format is on trial as a candidate for promotion to every synthesis.
