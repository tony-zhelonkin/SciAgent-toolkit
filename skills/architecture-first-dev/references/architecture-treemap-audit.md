# Architecture-Treemap Audit Specifications

**Status: probationary. Prune review: 2026-11-20.** This retrospective stack stays off the navigation autopilot. Each marked source section preserves the complete command specification and its original argument contract.

<!-- BEGIN SOURCE: commands/architect/components-extract.md -->
# Components-extract

Produce the **deterministic substrate** of an architecture audit: a `components.json` whose every entry traces back to the AST, the git history, or a metric tool — never to judgment. This is the foundation the rest of the `architecture-treemap` stack rests on.

The substrate biases attention toward hotspots (high churn, high fan-in, high complexity) and lets a reviewer sanity-check the synthesis's judgment against ground truth. It emits NO classifications, smells, background-knowledge edges, prune verdicts, or core/seam/removable boundary — those are not statically detectable and are authored later by `/audit-slice` and `/synthesize-audit`.

Pairs with `/architecture-treemap`, which renders a synthesized manifest.

## Phase 0: Parse arguments

- `$ARGUMENTS[0]` — path to the repository root to scan (required). Reject if missing: "Usage: /components-extract <repo-root>".
- `--out <path>` — output path (default: `<repo-root>/components.json`).
- `--since-days <N>` — churn window in days (default: 90).
- `--lang <python>` — source language for the import graph (only `python` is supported today; other values warn and skip the graph).
- `--include-nonsource` — include non-source scaffolding (default: OFF; see below).
- `--include-docs` — include documentation files (default: OFF).
- `--include-audit-output` — include the skill's own prior snapshots (default: OFF).

## Phase 1: Run the extractor

```bash
python skills/architecture-treemap/scripts/extract_components.py <repo-root> \
    [--out <components.json>] [--since-days 90] [--lang python]
```

The script is a pure function of the tree state at the current git HEAD plus the optional tools installed. It prints a one-line summary to stderr (files scanned, python files, static edges, unparseable files skipped, metric tools used, schema-validation result).

## What it emits

- **`physical_components[]`** — one per *architectural* source/config/test file under the repo. `.gitignore` is honored (via `git ls-files`); vendored and build dirs (`node_modules`, `dist`, `.venv`, `__pycache__`, `*.egg-info`, ...) are skipped. Each carries `id` (collision-free kebab-case of the path), `path`, `size_loc`, `kind` (module / config / test / doc / asset / other, by extension + location), `logical_owners: []` (judgment — left empty), a `metrics{}` block, and (if it sits in an import cycle) `cycle_id` + `cycle_size`.
- **`cycles[]`** — import cycles found by Tarjan strongly-connected-components over the static import graph (deterministic). Each is `{id, members[]}` for an SCC of size > 1. Teaches the Acyclic Dependencies Principle.

### Non-source inventory filter (default ON)

`git ls-files` honors `.gitignore` for UNTRACKED files only. This filter is the complement: it drops TRACKED scaffolding that is not architecture and would otherwise drown the signal. **Excluded by default:**

- license files (`LICENSE*`, `COPYING*`), lockfiles (`package-lock.json`, `poetry.lock`, `*.lock`, ...);
- root-level tooling/dotfiles (`.gitignore`, `.mcp.json`, `package.json`, `pyproject.toml`, ...);
- harness/editor/CI dirs (`.claude/`, `.github/`, `.vscode/`, `.idea/`, ...);
- data files by extension (`.csv`, `.tsv`, `.parquet`, `.rds`, `.npy`, `.pkl`, ...);
- **everything under a top-level `docs/` directory, by path-prefix** (this catches `docs/**/*.csv` that an extension-only filter misses).

**Kept:** source modules and tests (including `scripts/verify_*.py` and one-off scripts — those are real); and config/schema/asset files that live UNDER a source root (`src/` or a discovered package dir), e.g. `src/**/plugin_payload.schema.json` — a wire-contract schema belongs with the code it governs. Root-level tooling JSON is dropped; src-resident schema JSON is kept.

Pass `--include-nonsource` to restore the excluded set (e.g. when auditing a tooling repo where the configs ARE the architecture).
- **`edges[]`** — the Python intra-repo import graph at the physical level. Each edge is `type: "direct-call"`, `evidence_class: "static"` (renders SOLID), with an `evidence` note of the form `path:line -> import module`. Stdlib and third-party imports are dropped; only imports that resolve to a file in the repo become edges. `src/<pkg>/` layout is handled.
- **`logical_components: []`** — intentionally EMPTY. Logical grouping and classification are judgment; `/synthesize-audit` fills this. The schema requires the key but allows an empty array.
- **Top-level provenance** — `schema_version: "1.0"`, `audit_date` (today), `project` (repo dir name), `git_sha` (current HEAD), `snapshot_id`, and an `extractor{}` block recording `tool_version`, `extracted_at`, and which metric tools actually ran.

## The five locked metrics

Attached per physical component in `metrics{}`. All optional in the schema — a metric is omitted (never fabricated) when its tool is absent or it cannot be computed for that file.

| Metric | Source | Notes |
|---|---|---|
| `loc` | stdlib counter | Python: non-blank, non-comment-only lines. Other kinds: raw line count. |
| `fan_out` | `ast` import graph | Distinct intra-repo modules this file imports. |
| `fan_in` | `ast` import graph | Count of repo files importing this one. |
| `cyclomatic` | `radon` (optional) | Mean cyclomatic complexity across blocks. **Omitted gracefully if `radon` is not importable** — there is no hard dependency. |
| `churn_90d` | `git log --since` | Commits touching the file in the window. Skipped when the tree is not a git repo. |
| `test_ratio` | rollup | A repo-level test-LOC / source-LOC ratio attached to test files' metrics. Per-file test_ratio is not statically attributable; the true per-component value is a synth-level rollup once logical grouping exists. |

Two more deterministic metrics are derived from the import graph:

| Metric | Source | Notes |
|---|---|---|
| `instability` | derived | `I = fan_out / (fan_in + fan_out)`, in [0,1]. Undefined (omitted) for isolated nodes. Teaches the Stable-Dependencies Principle. |
| `refactor_pressure` | derived | `churn * cyclomatic * max(fan_in,1) / max(test_ratio, 0.05)`. **Absent when radon is** (cyclomatic missing). |

### Optional dependency — LOUD degradation

`radon` is optional (`pip install radon`) and runs only at extract time — it does not affect the portable HTML. But its absence is announced LOUDLY and honestly, never silently:

- The extractor prints a `WARNING` to stderr naming the unavailable metrics, and stamps them in `extractor.unavailable_metrics` (`["cyclomatic", "refactor_pressure"]`) with a per-metric `degradation_notes`.
- When radon is absent, a **clearly-labelled degraded proxy** `refactor_pressure_loc_proxy` is emitted (LOC/100 substituting for complexity) so the "look here first" lens still ranks something. The renderer shows it marked as a degraded proxy and, if even the proxy is unavailable, greys the lens with "install radon to enable complexity & refactor-pressure" rather than showing an empty lens.
- `jsonschema` is also optional; without it self-validation is SKIPPED with a warning, and the stdlib `validate_components.py` is the authoritative gate.

`cloc`/`tokei` are NOT depended on; the stdlib LOC counter is the deterministic floor.

### Scaling-regime hint

The stderr summary prints a regime hint keyed to component + edge counts: *small* (one `/audit-slice` + `/synthesize-audit` pass holds the whole graph) vs *large* (scope concerns narrowly, expect multiple passes, ~1M tokens). The deterministic substrate does not degrade with size — only the judgment layer does.

## Immutability and snapshot semantics

Each run is an **immutable, git-SHA-stamped snapshot**. `git_sha` records the HEAD the substrate was extracted at; `snapshot_id` and `audit_date` datestamp it. Re-running produces a NEW `components.json` — it never mutates a prior snapshot in place. To keep a history, write each run to a dated path (`--out audit/<date>/components.json`).

## Determinism

The extractor is deterministic modulo the wall-clock `extracted_at` field: same tree state + same tool set + same flags produces the same `physical_components`, `edges`, and metric values. File ordering is sorted; id disambiguation suffixes are assigned in sorted order.

## Self-validation

After writing, the script validates its output against `components.schema.json` if `jsonschema` is importable (best-effort; reported on stderr). The authoritative gate is `validate_components.py`, which `/architecture-treemap` runs before rendering and which refuses to render on any schema or referential-integrity error.

## Next steps

```
/audit-slice <hunch>     — author one judgment slice (human-triggered; refuses auto-piped input)
/synthesize-audit        — integrate slices + substrate into the full components.json (fills logical_components, classifications, smells, prune verdicts, core_boundary)
/architecture-treemap <path>  — render the synthesized manifest to a self-contained treemap.html
```

## Rules

1. **Substrate is facts, not assessment.** If asked to emit a classification, smell, or background-knowledge edge, refuse — those are judgment, authored by `/synthesize-audit`. Emitting a guess would defeat the design's trust scaffold.
2. **Omit, never fabricate.** A metric whose tool is absent is left out of the `metrics` block; it is never invented or zero-filled.
3. **One file out.** The extractor writes exactly one `components.json`.
4. **Never crash on a bad file.** A single unparseable source file is skipped and counted in the stderr summary; the run completes.
<!-- END SOURCE: commands/architect/components-extract.md -->

<!-- BEGIN SOURCE: commands/architect/audit-slice.md -->
---
status: probationary
prune_review_at: 2026-11-20
---

# Audit-slice

Dispatch N concern-driven slicers in parallel over EXISTING code. Each slicer takes one cross-cutting concern and traces it through code + tests + docs + ADRs in execution-time order, tagging coupling by Page-Jones connascence type. Output is one slice file per concern under `docs/_meta/architecture-audit/{date}/slices/`.

This is **Shape γ — a sibling to `/map`, not a replacement.** `/map` is feature-driven, one agent, all surfaces, cheap. `/audit-slice` is concern-driven, N agents, depth along each concern, expensive. Pick the one that matches the work shape.

## Off-the-autopilot gate (read before Phase 0)

`/audit-slice` is deliberately OFF the navigation autopilot. There is no nav-contract auto-fire for it; no other command recommends it as an automatic next step; it **REFUSES auto-piped input** from another agent. The trigger is USER-AUTHORED: a hunch, an endpoint, a feature name, a smell description — typed by a human. The act of a human naming the concerns IS the judgment gate. If the concern list arrives synthesized from another command's output rather than typed by the user, refuse:

```
/audit-slice takes user-authored concerns only. The act of naming the concerns is
the judgment gate — it must not be auto-piped from another agent's output. Please
type the concern(s) you want traced.
```

## Invocation banner (always print first)

Before doing anything, print:

```
/audit-slice is a retrospective architecture-audit tool. A full audit (N slices +
/synthesize-audit) typically runs ~0.8-1M tokens — comparable to a 6-reviewer panel,
but it processes N CONCERNS, not N perspectives. It earns its keep only with an
empirical anchor (a smoke/demo pass, an in-vivo regression, a pre-pivot OSS-core
decision). If you are scoping a NEW feature, run /map instead — it is ~20× cheaper.

Proceed with the audit?
```

WAIT for confirmation before Phase 1. If the user is clearly mid-audit (an audit dir for today already exists), state that and proceed without re-prompting.

## Phase 0: Parse arguments

- `$ARGUMENTS` — one or more concern names / hunches. Multiple concerns may be space-separated bare slugs (`lasso-state theme-bundle-surface discoverability`) or quoted phrases (`'the /export route still calls the legacy writer'`).

If no concern is given, ask:
```
Name the concern(s) to slice. A concern is a cross-cutting thread — a hunch, an
endpoint, a piece of shared state, a smell, or an emergent concept — NOT a feature
and NOT a file. Examples:
  /audit-slice lasso-state theme-bundle-surface
  /audit-slice 'the export path is dead weight'
```

Normalise each concern to a kebab-case slug for filenames; keep the original phrasing to pass to the slicer.

## Phase 1: Verify / create the audit directory

The audit parent is `docs/_meta/architecture-audit/{date}/` where `{date}` is today (`YYYY-MM-DD`).

- If it does not exist, create it plus a `slices/` subdir.
- If it exists (a re-run today), state that and continue — slices are additive; do not clobber existing slice files. If a slice file for a given concern slug already exists, ask the user whether to overwrite or pick a new slug.

Announce:
```
Audit dir: docs/_meta/architecture-audit/{date}/
Slices dir: docs/_meta/architecture-audit/{date}/slices/
Concerns ({N}): {list}
```

## Phase 2: Dispatch slicers in parallel

**CRITICAL:** Use a SINGLE assistant message with MULTIPLE Agent tool calls (`subagent_type: "slicer"`), one per concern — this is what makes them run in parallel. Do not dispatch sequentially.

Default model per slicer is the agent's default (`opus`). For a concern the user flags as shallow (a single dead route, a narrow rename), note in that slicer's prompt that sonnet-class depth is acceptable — keep the traversal to the spine.

Number slices `01..N` in the order concerns were given.

Prompt template (one per concern):
```
Trace the concern `{original concern phrasing}` through the codebase.

Output path: docs/_meta/architecture-audit/{date}/slices/{NN}_{concern-slug}.md

Trace in EXECUTION-TIME order: input → transform → second-writer/reconciliation →
consumer. Re-explore the code; verify line numbers at HEAD. Cover code + tests +
docs + ADRs that lie ON this concern's spine; exclude everything off-spine.

Tag every smell with its Page-Jones connascence type (Name/Type/Convention/
Execution-Order/Location/Identity/Value/Algorithm) and state which edge type
(direct-call / shared-state / background-knowledge) it maps to for synthesis.

Produce: traversal, branching points, findings (E-NN with file:line + connascence
tag, including negative-space findings), cross-slice touchpoints (name the other
concern slugs), decoupling proposals (R-NN), open questions for synthesis.

Diagnose, do not redesign. Do not assign core/seam/removable — that is the
synthesis's judgment. One file only.
```

## Phase 3: Write the slice-index README

After all slicers complete, write `docs/_meta/architecture-audit/{date}/README.md` (overwrite if a re-run today added slices):

```markdown
---
date: YYYY-MM-DD
git_sha: {HEAD short SHA}
status: probationary
---

# Architecture audit — {date}

Concern-driven slices. Each traces one cross-cutting concern in execution-time order.
This audit is the JUDGMENT layer; pair it with `/components-extract` (the deterministic
substrate) and feed both to `/synthesize-audit`.

## Slices
| # | Concern | Findings | Decoupling proposals | Connascence flagged |
|---|---------|----------|----------------------|---------------------|
| 01 | <a href="slices/01_{slug}.md">{concern}</a> | {E-count} | {R-count} | {tags} |
| ...|

## Cross-slice touchpoints
{Aggregate the touchpoints each slicer reported — where two concerns collide.}
```

## Phase 4: Route to synthesis

```
Audit slices complete ({N} files) under docs/_meta/architecture-audit/{date}/slices/.
Slice index: docs/_meta/architecture-audit/{date}/README.md

Per-slice headlines:
  01 {concern}: {first finding's headline}
  ...

Recommended next step:
  /synthesize-audit               — integrate these N slices + the /components-extract
                                    substrate into the full renderable components.json

If you have NOT yet run /components-extract on this repo, run it first — the synthesis
merges the deterministic substrate (physical files, static edges, metrics) with the
slices' judgment. Without the substrate the synthesis has no trust scaffold.
```

## Rules

1. **User-authored trigger only.** Refuse auto-piped concern lists. The human typing the prompt is the gate.
2. **Parallel dispatch is mandatory.** One message, multiple Agent calls. Sequential dispatch violates the design.
3. **Slicers re-explore; reviewers do not.** This is the distinguishing capability — do not constrain a slicer to read-only artifacts the way `/review` constrains reviewers.
4. **One concern per slicer, one file per slice.** If a slicer writes elsewhere or merges concerns, flag it.
5. **Additive, never clobbering.** A re-run today appends slices; never overwrite an existing slice without asking.
6. **This command is probationary** (`prune_review_at: 2026-11-20`). It is the right tool ~once per quarter per major project, not per-feature. Police the trigger; if you find yourself reaching for it to scope a new feature, you wanted `/map`.
<!-- END SOURCE: commands/architect/audit-slice.md -->

<!-- BEGIN SOURCE: commands/architect/synthesize-audit.md -->
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
<!-- END SOURCE: commands/architect/synthesize-audit.md -->

<!-- BEGIN SOURCE: commands/architect/architecture-treemap.md -->
# /architecture-treemap

Render a `components.json` audit manifest into a self-contained interactive
`treemap.html`. Pure transform: no LLM, no network, no synthesis. Same input →
same output.

**Pair with `/diagram` for per-feature Mermaid; use this command for cross-
feature audit-level visualisation.**

## Usage

```
/architecture-treemap <path-to-components.json> [--output <path-to-html>]
```

- `<path-to-components.json>` — required. Path to a valid `components.json`
  that conforms to the architecture-treemap schema.
- `--output <path>` — optional. Default: `<input-dir>/treemap.html`.

If the path is omitted the command looks for the most recent directory under
`docs/_meta/architecture-audit/` matching `YYYY-MM-DD` and uses its
`components.json`.

## Phase 0: Resolve input path

1. If `$ARGUMENTS[0]` is supplied, use it as `<input>`.
2. If omitted, glob `docs/_meta/architecture-audit/????-??-??/components.json`.
   Pick the lexicographically largest directory (ISO date sorts as chronological).
   If none found:

   ```
   No components.json found. Tried:
     - docs/_meta/architecture-audit/YYYY-MM-DD/components.json

   The components.json must conform to the schema at
     skills/architecture-treemap/components.schema.json

   Run /synthesize-audit to generate one, or author it manually.
   ```

   Stop.

## Phase 1: Validate

Run the validator before rendering. Print results and stop on error.

```bash
python skills/architecture-treemap/scripts/validate_components.py <input>
```

On schema error → print errors and stop. Half-rendered treemaps are worse than none.

Referential-integrity warnings are printed but do not block rendering.

## Phase 2: Render

```bash
python skills/architecture-treemap/scripts/render_treemap.py <input> [--output <output>]
```

The renderer:
- Injects the validated JSON into the HTML template.
- Inlines D3 v7 (vendored asset — no CDN required).
- Inlines the `treemap.bundle.js` renderer.
- Writes one self-contained `.html` openable via `file://`.

## Phase 3: Report

Report to the user:

```
Rendered: <absolute-path-to-html>
Open with: file://<absolute-path>

Views available:
  Logical  — treemap by size_estimate_loc, coloured by classification
  Physical — treemap by size_loc, coloured by primary logical-owner classification
  Graph    — force-directed layout with three classification-gravity clusters

Edge overlay (Logical + Graph):
  Click a cell/node → edges for that component light up.
  Shift-click adds. Background-click clears.
  "Show all edges (faint)" toggle in header.

Shift-click nodes in Graph view to compare edge sets.
The force-graph seed is in the URL hash (#<integer>) for reproducible screenshots.
```

If the manifest has no `logical_components` (empty array), note:

```
NOTE: logical_components is empty. The Logical and Graph views will show a
"No logical model yet" notice. Run /synthesize-audit to populate them.
The Physical view renders from physical_components alone.
```
<!-- END SOURCE: commands/architect/architecture-treemap.md -->
