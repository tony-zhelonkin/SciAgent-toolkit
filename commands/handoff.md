# Handoff

Write the session's record: why the work was undertaken, what it was designed
around, what it achieved, and where it stopped.

**You write this yourself, in this session.** The reasoning that makes a handoff
worth reading — the premise you started from, the failure you were designing
around, the reading you rejected — exists only in this context. It cannot be
recovered from a diff, so it is not something to delegate to a fresh reader.

## When to use

- Ending a session, or reaching a point where the work would survive a break.
- After a stage runs to completion, a decision is made, or a finding lands.

## When not to use

- Recording a single decision and its evidence → `reasoning-trace` writes the
  topic note. A handoff cites those notes; it does not restate them.
- Narrating work still in flight. Write it when there is an outcome to state.

## What it writes

One live file, plus the archives it has superseded:

| Scope of the work | File |
|---|---|
| One stage | `docs/_internal/<stage-stem>/session.md` |
| Spans or precedes stages | `docs/_internal/_project/session.md` |

Resolve the scope from what the session touched. If the work genuinely spans
one stage and part of another, ask which scope owns it rather than guessing —
the wrong scope buries the record where nobody looks for it.

**Keep the copy you supersede.** Before writing, move the current `session.md`
to `session-<date>.md`, where `<date>` is the last day that copy was actual —
`date +%F` at the time of the handoff. A second handoff the same day appends a
counter: `session-2026-08-24-2.md`. Then write the new `session.md`.

```bash
d=$(date +%F); n="session-$d"
[ -e "docs/_internal/<scope>/$n.md" ] && n="$n-2"      # and -3, …
mv docs/_internal/<scope>/session.md "docs/_internal/<scope>/$n.md"
```

The live file is what a reader opens; the archives are the history of intent,
each one readable on its own terms and stamped with the period it described. Git
holds the diffs, but a diff does not tell a reader when a claim stopped being
true, and nobody reconstructs a superseded premise from a patch.

`docs/_internal/` is its own git repository — commit there as part of the
handoff. If it is not a repository yet, say so and stop: the archive would have
nowhere durable to live, and losing the previous session's record is the one
failure this command must not cause.

Nothing else is modified. This command does not commit the parent repository.

## Shape

Sixty lines is the ceiling, and most sessions need fewer. A handoff nobody
finishes reading has failed at its only job. Cut history the git log already
holds, cut speculation, cut restatement of the topic notes you are citing.

```markdown
# <scope> — <one line: what this session was for>

## Premise
What was believed to be true when the work started, and what prompted it.

## Designing around
The failure or constraint the work had to accommodate. Name it concretely:
the measurement, the limit, the thing that broke last time.

## Intended result
What success looked like when the work began — and, if it moved, what it
became and why.

## Where it stands
Current point, last thing completed, immediate next action. Three lines.

## What changed
- `02_analysis/stages/NN_<name>.R` — purpose → `03_results/<stage>/...`
- Flag `[UNCOMMITTED]` on anything run but not committed.
- Flag artifacts under `03_results/` that no `README.md` captions yet.

## Decisions
- **<call made>:** one sentence → `docs/_internal/<scope>/<topic>.md`
- `[NO RECORD]` on any non-trivial call whose reasoning was never written.
  That flag is the point: an unrecorded decision is not reproducible.

## Traps
Only what cost time and would cost it again. Exact numbers, exact paths.

## Next
1. Actionable task, with the command to run.
```

Omit a section that has nothing in it. An empty heading is noise.

## Before finishing

- Every path exact; every number measured rather than recalled.
- Premise, Designing around, and Intended result are filled from this session's
  reasoning. If you cannot state them, say so explicitly instead of inventing
  them — a fabricated premise is worse than an absent one.
- Each non-trivial decision either cites a topic note or carries `[NO RECORD]`.
- Under sixty lines.
- The superseded `session.md` kept as `session-<date>.md`, not overwritten.
- `docs/_internal/` committed.

Report back: the file written, the next action, and the count of `[UNCOMMITTED]`
and `[NO RECORD]` flags raised.
