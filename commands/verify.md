# Verify

Mechanical drift check: compare what `plan/phase-NN.md` and `design/03-decisions.md` said would exist against what's actually in the code. Produces `docs/{feature}/verify.md` with a Verdict (`CLEAN | INCOMPLETE | DRIFT | NEEDS REVIEW`).

Main-agent only; no sub-agent dispatch. Uses `Read`, `Grep`, `Glob`, `Bash` (for mtime / file existence) against the live codebase. Designed to catch the bug class you can't eyeball across thousands of lines: "did we build all N phases, and does each ADR have at least a grep-able footprint?"

Use after `/implement` phases land. Not meant to replace code review — meant to catch **incompleteness** and **scope drift** that human review misses because it's boring and mechanical.

## Phase 0: Parse arguments

- `$ARGUMENTS[0]` — feature slug (required). Reject if missing: "Usage: /verify <feature>".
- `$ARGUMENTS[1+]` — optional phase range (e.g., `1-3`, `5`). If present, verify scopes to those phases only. Empty → verify every phase in `plan/`.

## Phase 1: Verify prerequisites

1. **`docs/{feature}/plan/` MUST exist with at least one `phase-NN.md`.** If absent:

   ```
   docs/{feature}/plan/ not found. /verify compares plan vs code; run /plan {feature} first.
   ```

   Stop.

2. **`docs/{feature}/design/03-decisions.md` MUST exist.** If absent:

   ```
   docs/{feature}/design/03-decisions.md not found. /verify checks ADR compliance; run /design {feature} first.
   ```

   Stop.

3. **`docs/{feature}/map.md` MUST exist.** Used to scope the "feature area" for scope-drift detection. If absent, warn and proceed with phase-named files only as the scope set.

**Scribe-on-latest read rule.** When reading `map.md`, `synthesis.md`, or `design/03-decisions.md`, also read any `## User resolutions`, `## Binding chat decisions`, or `## Post-READY amendments` section present. A decision in these sections that names a grep-able anchor (constant, function, column, error phrase) is treated like a MADR consequence: derive an anchor and verify it appears in code.

## Phase 2: Inventory what the plan promised

Read every `docs/{feature}/plan/phase-NN.md` in scope. From each, extract:

- **Files to Create** — every path under an `## Files to Create` or `### Files to Create` heading.
- **Files to Modify** — every path under an `## Files to Modify` or `### Files to Modify` heading.
- **Verification items** — lines under `## Verification` heading (test file paths, invariants named).

Parse defensively: markdown headings vary slightly across features. Match case-insensitively; accept `H2` or `H3`; accept variants `files to create`, `files-to-create`, `Creates`, etc. If a phase file has neither section (unusual), note "phase-NN has no files section; skipped" and continue.

Build a normalised in-memory table:

```
phase | path                          | action  | source
------|-------------------------------|---------|----------
1     | pathway_explorer/schema.py    | create  | plan/phase-01.md
1     | pathway_explorer/data_loader.py | modify | plan/phase-01.md
2     | pathway_explorer/similarity.py | modify | plan/phase-02.md
...
```

## Phase 3: Inventory what the design promised

Read `docs/{feature}/design/03-decisions.md`. For each `## ADR-NNN` heading:

1. Extract the ADR number, title (from the `## ADR-NNN — {title}` line), and the **Decision** section body (the one-or-two-sentence declarative statement).
2. Derive a **checkable anchor** — the most grep-able string the Decision implies. Heuristic:
   - If the Decision names a function signature (e.g., `classify_member(...)`), anchor = that function name.
   - If the Decision names a constant (`BUNDLE_Z_THRESHOLD = 2.0`), anchor = the constant name.
   - If the Decision names a column/field name (`genes_full_set`, `within_method_z`), anchor = that exact string.
   - If the Decision names a file (`config.py § ENTITY_PROFILES`), anchor = both the file and the symbol.
   - If the Decision is purely behavioral (e.g., "hard-fail if X missing"), anchor = the phrase that would appear as an assertion or error message (`raise KeyError`, `assert within_method_z`, etc.).
3. If no mechanical anchor can be derived, mark the ADR `NEEDS_HUMAN_REVIEW` — these are behavioral-only Decisions (like "colour palette should be perceptually uniform") that a grep can't verify.

Also read any "Inherited meta-design constraints" section and treat each inherited MADR bullet as an additional row (with anchor derived from the MADR's consequence-for-this-feature line).

Build the ADR table:

```
adr     | title                                  | anchor                | grep target        | verifiable?
--------|----------------------------------------|-----------------------|--------------------|-----
ADR-001 | Similarity on genes_full_set           | genes_full_set        | similarity.py      | yes
ADR-003 | Hard-fail on within_method_z           | within_method_z + raise | data_loader.py   | yes
ADR-014 | Concordance ship-gate r ≤ 0.7          | signal_concordance + 0.7 | similarity.py   | yes
ADR-009 | Top-K=12 stratified table              | TOP_K or top_k = 12   | cluster panel code | yes
MADR-001| within_method_z canonical column       | within_method_z       | (all consumers)    | yes
ADR-004 | 4 guardrails on cluster-union export   | (behavioral cluster)  | —                  | NEEDS_HUMAN_REVIEW
...
```

## Phase 4: Check file existence + mtime

For every row in the phase table, check the file on disk:

- **`create` action**: file should exist now. If missing → `MISSING`. If present → `OK`.
- **`modify` action**: file should exist AND have been modified since the plan was drafted (approximation — compare mtime against plan's frontmatter `date:`). If unchanged → `NOT YET IMPLEMENTED`. If modified later → `OK`. If file is absent entirely → `MISSING` (the plan assumed it existed to modify, so absence is a pre-existing data issue).

Use `Bash` for mtime:

```bash
stat -c '%Y' <path>   # Linux
```

Convert to ISO date; compare against `plan/phase-NN.md` frontmatter `date:` (parse YAML frontmatter). Tolerate ±1 day slop (filesystem timestamps can be stale after git operations).

## Phase 5: Check ADR anchors

For every ADR/MADR row with `verifiable: yes`, run `Grep` for the anchor in the feature's scope files (as named by `map.md § Key files` and the phase tables). Record:

- **FOUND** — anchor appears at least once in a file the plan named.
- **MISSING** — anchor nowhere in the scope. Red flag.
- **AMBIGUOUS** — anchor appears but only in a comment/docstring/test fixture (not in executable code). Flag as weak evidence.

For `NEEDS_HUMAN_REVIEW` rows, don't grep; they go to the "Open questions for human review" section of the output.

## Phase 6: Scope-drift detection

Collect the set of files:

- **Plan scope**: every path named in the phase tables (Phase 2).
- **Map scope**: every path under `map.md § Key files` and `map.md § Blast radius` (best effort — parse the markdown headings).

Union these = **expected feature area**.

For every file in the expected feature area that's been modified recently (mtime within the plan's drafting window to now), check that the file appears in the plan's Files-to-Create OR Files-to-Modify lists. If a file was modified but is NOT in any phase's lists → **scope drift** candidate.

Use `Bash(find)` with a date predicate for the mtime scan; exclude `.git/`, `docs/`, `tests/` fixtures, build artifacts (`__pycache__`, `*.pyc`, `dist/`, `build/`), and lock files (`uv.lock`, `poetry.lock`).

Note: scope drift isn't automatically a defect — helper functions often get added to existing files during implementation. The goal is to *surface* them so the user can decide whether to backfill the plan or revert.

## Phase 7: Test inventory

Parse every phase's `## Verification` section for test-file references (e.g., `tests/test_similarity.py`, `tests/test_schema.py`). Check each:

- File present? If no → `MISSING`.
- If present, count tests (`grep -c "^def test_"` for pytest; adapt for other frameworks if named in map.md).
- Compare against any "N unit tests" count promised in the plan. Missing tests → flag.

## Phase 8: Compute verdict

Apply the following precedence (first match wins):

1. **INCOMPLETE** — any `create` file is MISSING, OR any `modify` file is `NOT YET IMPLEMENTED`, OR any promised test file is missing.
2. **DRIFT** — scope drift detected (≥1 unexplained file in the feature area), OR ≥1 ADR anchor is MISSING (the code doesn't contain something the design said it would).
3. **NEEDS REVIEW** — INCOMPLETE/DRIFT both clean, BUT >30% of ADRs were flagged `NEEDS_HUMAN_REVIEW` (i.e., the design is behavior-heavy and verify couldn't reach conviction on most of it).
4. **CLEAN** — everything else.

## Phase 9: Write verify.md

Write `docs/{feature}/verify.md`. One file; overwrite if present. Format:

```markdown
---
feature: {slug}
date: YYYY-MM-DD
verdict: CLEAN | INCOMPLETE | DRIFT | NEEDS REVIEW
phases_scoped: 1-N (or comma list from $ARGUMENTS[1])
---

# Verify — {feature}

## Verdict: {VERDICT}

{One-sentence state summary, e.g., "8 of 12 phases complete; no scope drift; 3 ADRs need hand-verification."}

## Phase completion

| Phase | File | Action | Status | Note |
|---|---|---|---|---|
| 1 | `pathway_explorer/schema.py` | create | ✅ OK | modified 2026-04-18 |
| 1 | `pathway_explorer/data_loader.py` | modify | ✅ OK | modified 2026-04-18 |
| 2 | `pathway_explorer/similarity.py` | modify | ⚠️ NOT YET IMPLEMENTED | mtime pre-dates plan |
| 3 | `pathway_explorer/ubiquity.py` | create | ❌ MISSING | — |
| ... | ... | ... | ... | ... |

**Summary**: {N} of {M} files at expected state. {P} missing, {Q} not yet implemented.

## ADR compliance (mechanical anchors)

| ADR | Title | Anchor | Found | Status |
|---|---|---|---|---|
| ADR-001 | Build similarity on genes_full_set | `genes_full_set` | ✅ `similarity.py:47, 82, 113` | OK |
| ADR-003 | Hard-fail on missing within_method_z | `raise KeyError.*within_method_z` | ✅ `data_loader.py:92` | OK |
| ADR-014 | Concordance ship-gate r ≤ 0.7 | `signal_concordance` | ❌ not found | MISSING |
| MADR-001 | within_method_z canonical column | `within_method_z` | ✅ 14 matches across 5 files | OK |
| ADR-004 | 4 guardrails on cluster-union export | (behavioral) | — | NEEDS_HUMAN_REVIEW |
| ... | ... | ... | ... | ... |

**Summary**: {N}/{M} mechanically-verifiable ADRs OK; {P} MISSING; {Q} need human review.

## Scope drift

{If any files modified in the feature area are not in the plan:}
| File | First modification | Note |
|---|---|---|
| `pathway_explorer/utils.py` | 2026-04-18 | Added helper function `_normalise_scale`; not in any plan phase. |
| ... | ... | ... |

**Recommendation**: {backfill plan / revert / accept as incidental helper}.

{If no drift, write "No scope drift detected — every modified file in the feature area was named by a plan phase."}

## Test inventory

| Test file | Promised | Present | Count (promised / actual) | Status |
|---|---|---|---|---|
| `tests/test_similarity.py` | 8 unit tests | ✅ exists | 10 / promised 8 | OK (exceeds) |
| `tests/test_schema.py` | 6 unit tests | ❌ missing | — / 6 | MISSING |

## Open questions for human review

ADRs that can't be verified mechanically. Each needs a hand-read or a reviewer pass:

- **ADR-004 — 4 guardrails on cluster-union export.** The guardrails are surfaced as UI panels + inline text; whether all four appear correctly and are legible requires a browser pass. Recommended: smoke-test one dashboard HTML manually, or dispatch `/review {feature} --as graphic,divergent`.
- **ADR-009 — stratified column layout with Top-K=12.** K is present in the code (`TOP_K = 12`), but whether the stratified breakdown semantically matches the ADR's intent requires a domain read. Dispatch `/review {feature} --as bioinf` if in doubt.
- ...

## Recommended next action

{one sentence, based on verdict:}

- **CLEAN** → `/implement {feature} {next-phase}` to continue, or proceed to the next feature's `/plan`.
- **INCOMPLETE** → `/implement {feature} <N>` for each missing phase. Re-run `/verify` after.
- **DRIFT** → hand-read the drifted files; either (a) backfill the plan with a new phase, (b) revert the unexplained changes, or (c) promote the drift to an ADR if it's a real design decision that was made during implementation.
- **NEEDS REVIEW** → `/review {feature} --as divergent,code-reviewer` for the behavioural-quality lens, or hand-read the ADRs flagged in "Open questions."
```

## Rules

1. **Read-only.** `/verify` writes exactly one file: `docs/{feature}/verify.md`. It reads the codebase and the docs; it does not edit source. Any suggested fixes are surfaced for the user to apply.
2. **Mechanical, not semantic.** This command checks "does the anchor exist in the code?" not "is the code correct?" Correctness belongs to `/audit` (not yet built) or the existing `/review` panel pointed at the code.
3. **Be honest about coverage.** Every ADR flagged `NEEDS_HUMAN_REVIEW` must appear in the "Open questions for human review" section. Never silently skip an ADR by failing to derive an anchor.
4. **False positives are acceptable; false negatives are not.** If `/verify` flags a drift that turns out to be a helper function the user is fine with, that's cheap to dismiss. If it *misses* a missing file, the user may ship INCOMPLETE thinking it's done. Err toward flagging.
5. **No auto-fixes.** `/verify` never modifies the code or the design to bring them into alignment. That's a human call — either update the design (the reality won) or revert the code (the design won).
6. **Do not re-dispatch mapper.** Mapper produces a full map.md; that's expensive and duplicative. `/verify` does targeted greps directly. If the user wants a fresh whole-codebase map, `/map {feature}` is the right command.
7. **Scope honestly.** Files under `docs/`, `tests/` fixtures, and build artifacts are never scope-drift candidates. Only source files the plan would have scoped.
8. **Frontmatter is load-bearing.** The `verdict:` field in the output is machine-parseable so `/status` can aggregate it across features. Use exactly one of the four values.