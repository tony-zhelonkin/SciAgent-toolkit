# _project — scio hardening, the dual-name mechanism, three granted decisions

## Premise

Opened to finish plan #47 phase 02a. Became a sweep for one defect class: **an
instruction naming something the mechanism does not guarantee.** Late on, the owner
granted three decisions previously blocked.

## Designing around

A claim with no mechanism behind it drifts silently, because nothing checks. And a
fix that reads correct can still be measured wrong: phase 01 recorded
`01_scripts/` as existing nowhere, and acting on that made eight projects right and
two wrong.

## Where it stands

24 commits on `dev`, 62 ahead of `origin/dev`. Tree clean, **72 tests passing**,
`lint --check toolkit --strict` clean. One warning left, `docs-layout` on
`docs/_internal/`, which decision 1 clears. Nothing half-done.

Headline change: one list, `lib/scio/common.sh::_SCIO_TOOLKIT_DIRS`, with staged
discovery so evidence outranks the name. **The rename is no longer a flag day** —
that is what made it look larger than it was. Also: delegate-cli rebuilt around
`status.sh`, `--bg` gone; five new checks; `session-history/` as the memory shape.

## THE THREE DECISIONS — granted, NOT executed

Each has a note with measured state, steps and traps. Read it before touching.

1. `decision-1-memory-repo.md` — `docs/_internal/` becomes its own repo. **Do
   first**, it shares a tree with decision 2. One owner answer needed: delete the
   107 MB of codex logs, or move them out?
2. `decision-2-rename.md` — rename to `scio`. **Four names, not one**: GitHub repo,
   local bare hub, this working copy, consumer vendor dirs. `gh` is authed as
   `tony-zhelonkin`. Do NOT bump scbio-docker's pin `3ab37688`.
3. `decision-3-fleet-repin.md` — push to 8 named consumers. First authorised write
   outside this repo. Pilot `DC-nexus`, review, then fan out.

## Decisions recorded

- Interaction stance → `interaction-stance-home.md`: belongs in dev-env, scio
  changes nothing.
- CRAFT character budget → declined (#36); the line cap bounds shape and says so.

## Traps

- `craft` refuses a hand-edited block rather than clobbering it. Read the local
  edits; that refusal is not a cue for `--force`.
- Editing a hook template obliges `tools/gen-template-provenance.sh`, or every
  consumer holding the new bytes is read as user-authored and never refreshed.
- `13036-DM_DMlab_summer_2025` vendors at `01_Scripts/`, capital S.
- `handoff/session-state.md` is a second continuity record; migrate it in
  decision 1.

## Next

1. Answer the log question, then execute decisions 1, 2, 3 in that order.
2. **Reply to the bulkiRNA agent** — `bulkirna-agent-reply.md`. Still waiting,
   **not** unblocked, and its premise is false: bulkiRNA is absent from the image
   running 7 of 8 containers, so `coresh_*` reaches none of them. `#47` ph2 and
   `#50` stay gated on `#49`, which the owner deferred.
3. Unaudited for the defect class: `templates/project/` and the 20 agent files.
