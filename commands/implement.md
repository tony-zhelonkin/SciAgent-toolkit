# Implement

Execute implementation phases for a feature. Default: one phase per invocation (hand-paced, matches the hand-edit workflow). Opt-in `--auto` mode: run end-to-end through all remaining phases, stopping only at failure, uncertainty, or completion.

## Phase 0: Parse arguments

Accepted shapes:

| Invocation | Mode | Behaviour |
|---|---|---|
| `/implement <slug> <N>` | **single-phase** (default) | Executes exactly phase N, reports completion, stops. Today's behaviour — unchanged. |
| `/implement <slug> --auto` | **end-to-end** | Runs from the first unfinished phase onward until failure, stop-at-uncertainty, or all phases complete. One-liner per completed phase. |
| `/implement <slug> <N> --auto` | **end-to-end from N** | Runs from phase N onward with the same end-to-end semantics. |
| `/implement <slug>` (no phase, no flag) | **clarify** | Lists available phases from `docs/{feature}/plan/` and asks whether to run `1` (single) or `--auto`. Do NOT silently pick a mode. |

`--auto` and a phase number can combine (`<slug> <N> --auto`) but `<slug> --auto <N>` is also accepted (order-independent flag parsing).

If the slug is missing, reject: `Usage: /implement <slug> [<N>] [--auto]`.

## Phase 1: Load phase context

Read in order:
1. `docs/{feature}/plan/phase-NN.md` — what this phase does (for single-phase mode, only N; for `--auto`, load the starting phase's doc; subsequent phase docs load lazily per iteration).
2. `docs/{feature}/plan/README.md` — phase ordering, verify prior phases are done.
3. `docs/{feature}/design/01-architecture.md`, `02-behavior.md` — structural + behavioral reference.
4. `docs/{feature}/design/03-decisions.md` — constraints to respect.
5. `docs/{feature}/map.md` — existing patterns to follow.

**Do NOT re-grep for orientation.** The map is your source of truth for where things are.

**Scribe-on-latest read rule.** For every upstream artifact above, also read any `## User resolutions`, `## Binding chat decisions`, or `## Post-READY amendments` section present. Treat those as equal-priority to design ADRs and plan checklists; cite them explicitly when they shape an implementation choice.

## Phase 2: Verify prior phases are done

For phases `1..N-1` (where `N` is the phase about to run), check the phase file for a "Completed" marker or verify key files listed as "Files to Create" exist and are non-empty.

If prior phases are incomplete, stop:
```
Phase {N-1} does not appear complete ({reason}). Complete previous phases first.
```

In `--auto` mode: the starting phase is the first unfinished one, so Phase-2 gating reduces to "confirm phases 1..(start-1) are complete." If they aren't, stop and surface which.

## Phase 3: Execute

### Single-phase mode

Announce:
```
Starting Phase {N}: {focus}
Files to create: {list}
Files to modify: {list}
```

Then:
1. Create/modify files as specified in the phase doc.
2. Follow existing code patterns from `map.md` (naming, structure, imports).
3. Keep the surface area minimal — only what the phase specifies.
4. Run the project's tests if test infrastructure exists (check `map.md` for test locations).
5. Update the phase file: mark each "Files to Create" / "Files to Modify" as done.

Then Phase 4 (verify) and Phase 5 (handoff) below.

### `--auto` mode

Announce once:
```
Starting /implement {slug} --auto
  Plan: {N} total phases, starting at phase {start} ({focus})
  Stop conditions: phase failure | stop-at-uncertainty | all phases complete
```

Then loop `phase = start; phase ≤ N; phase++`:

1. Load `phase-NN.md` (if not already loaded).
2. Execute the phase: create/modify files, follow patterns, run phase-scoped tests.
3. Run the phase's Verification checklist.
4. **Decision point** — exactly one path applies:
   - **All Verification items pass** → emit a one-liner, `phase += 1`, continue loop:
     ```
     Phase {N} ✅ — {created_count} files created, {modified_count} modified, {test_count} tests pass
     ```
   - **Any Verification item fails** → emit a one-liner with the failure, stop the loop, go to Phase 4 (verify) and Phase 5 (handoff) below:
     ```
     Phase {N} ❌ — {failure_reason}
     ```
   - **Stop-at-uncertainty triggers** (phase doc is ambiguous, pre-existing state unexpected, etc.) → same treatment as failure; surface the ambiguity and stop:
     ```
     Phase {N} ⚠️ stopped — {ambiguity}
     ```
5. Between phases in `--auto`, do NOT emit the full Phase-4 block (see below). The one-liner IS the per-phase report. Full detail stays in the phase file for audit.

When the loop completes (all phases green):
```
/implement {slug} --auto complete — {N} phases shipped.

Per-phase summary:
  Phase 1 ✅ — ...
  Phase 2 ✅ — ...
  ...

Next: /verify {slug}   (mechanical drift check against plan + design)
```

## Phase 4: Verify against checklist (single-phase OR post-stop in --auto)

For each check in the phase doc's Verification section, confirm it passes. Report:
```
Phase {N} complete.
✅ {check 1}
✅ {check 2}
{⚠️ or ❌ for anything that didn't pass cleanly}
```

If anything failed, do NOT mark the phase complete. Stop and ask the user how to proceed.

In `--auto` mode this block fires exactly once — when the loop stops (either at completion of the final phase, or at a failure/uncertainty mid-run). It does NOT fire between phases.

## Phase 5: Hand off to next phase (single-phase mode only)

```
Phase {N} done. Next: Phase {N+1}: {focus}.

Run /implement {feature} {N+1} when ready.
  (or /implement {feature} {N+1} --auto to run from here to the end.)
```

**WAIT** for user to trigger the next phase. Do NOT continue automatically.

In `--auto` mode: do NOT emit this block. The `--auto` loop already decided what's next.

## Rules

1. **Default is one phase per invocation.** `--auto` opts into end-to-end. Even in `--auto`, stop conditions (failure, uncertainty, all-complete) still fire — the flag removes gates between successful phases, not safety checks.
2. **Follow the plan.** If a phase doc seems wrong, stop and ask — do NOT silently deviate. Re-run `/design` and `/plan` if the design is actually off.
3. **Minimal surface area.** Do not add "while we're here" refactors or "nice to have" helpers beyond what the phase specifies.
4. **Match existing patterns.** Code style, naming, structure all from `map.md`.
5. **No speculative error handling.** Handle the error paths that behavior docs specify. Do not add defensive code for scenarios that can't happen.
6. **Stop at uncertainty.** If the phase doc is ambiguous, ask the user rather than guessing. This rule is unchanged by `--auto` — ambiguity still halts the loop.
7. **`--auto` output discipline.** Per-phase output is a single line: status glyph, phase number, brief count/failure reason. Do NOT emit the full Phase-4 verification block between phases; the compressed line IS the report in `--auto` mode. Full detail lives in the phase file.
8. **`--auto` never auto-runs `/verify`.** The final recommendation is to run `/verify` as a separate step. `--auto` is the implementation loop, not an end-to-end pipeline runner.
