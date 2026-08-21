# Consultation — is the `docs/_internal/` structure carrying its own weight?

Read-only. Produce a written recommendation. **Modify nothing.**

The owner's question, in his words: *"consult if we are not over-complicating the
`_internal` structure of documents, so that we would be sure that all folders we
have are actually doing some form of work or bear reason, and at least intent,
and not just there for no matter."*

Take that literally. **The default answer for any directory should be "delete
it".** Argue each one back in, or don't.

## Read first

- `docs/proposals/2026-08-11-offline-distribution/50_ADRs.md` — **ADR-D10**, just
  written. It is the thing under review.
- `lib/scio/lint.sh` — the `internal-memory` check (search
  `_lint_check_internal_memory`). This is the only enforcement, and it defines
  operationally what the structure *is*.
- `docs/_internal/research/2026-08-20-internal-skeleton/FINDINGS_field.md` — what
  24 real projects actually contain. Cite it; do not re-derive it.
- `templates/project/_common/docs/_internal/README.md.template`,
  `templates/project/analysis/AGENTS.md.template`,
  `templates/project/software/AGENTS.md.template` — the instructions as they
  now read, after today's deletions.

## The named locations, as of today

| Location | Claimed job |
|---|---|
| `docs/_internal/<stage-stem>/session.md` | where work on that stage stopped, updated in place |
| `docs/_internal/<stage-stem>/reasoning/<topic>.md` | one decision, written before proceeding past it |
| `docs/_internal/_project/session.md` | the same, for work spanning stages |
| `docs/_internal/_project/reasoning/<topic>.md` | the same, for reasoning spanning stages |
| `docs/_internal/_project/plans/<date-slug>/` | a phased plan: `00_INDEX.md` + `NN_<slug>.md` |
| `docs/_internal/scientific-context.md` | stable framing, changes rarely (analysis template only) |
| `docs/_internal/README.md` | the local contract, a few lines |
| `docs/_internal/<work-stem>/` | the same grammar in a repo with no stages (software template) |

Plus, in this toolkit's own tree (which is not a consumer project and runs no
`link`): `handoff/`, `plans/`, `research/`.

## What is settled — do not reopen

- **No new hook.** Not Claude, not Codex, not git.
- **Never scaffold empty.** A directory arrives with its first real file.
  `.gitkeep` shipped in five projects and was ignored in all five.
- **No hand-maintained metadata.** The owner rejected intent-carrying config
  knobs from lived experience: they produced churn and variables nobody could
  justify. Do not propose a manifest, an index file agents must maintain, a
  front-matter field, or a status enum.
- **Do not name a location nothing creates.** Root `_scratch/` was named as "the
  only sanctioned throwaway zone" and exists in no project; the claim was deleted
  today rather than relocated.
- **Stage-keying is the key** because it is the number the human already reasons
  in. Anything requiring a new vocabulary is worse.
- Enforcement is `scio lint --check internal-memory`, opt-in, bash only, reading
  the tree.
- Simplicity binds. The owner: *"I really don't want this make more complex."*
  **Deleting counts as progress.**

## Questions

1. **Which of the eight locations earn their existence?** For each, name the
   failure that occurs if it does not exist. If you cannot name one, say delete.
2. **Does `reasoning/` earn being a subdirectory** rather than
   `<stage-stem>/<topic>.md` files sitting beside `session.md`? A directory that
   usually holds one file is a promise, and today's whole exercise was about
   deleting promises.
3. **Does `_project/` earn itself?** Alternatives: no such directory, with
   cross-stage notes living at the `docs/_internal/` root; or a `00_` stem, since
   the stage grammar already reserves `00`–`09` for preliminaries.
4. **Is `plans/<date-slug>/{00_INDEX,NN_slug}.md` over-structured?** The owner
   just ruled it lives at `_project/plans/<date-slug>/` and that `00_STATE.md`
   and `session.md` are one form under two names, so one spelling wins. Which
   spelling, and does a plan need a directory at all rather than one file that
   grows? Evidence both ways: 14782-DM's six-phase plan is the fleet's best
   working record; this toolkit's own eight earlier plan directories sat
   untracked and unread until today.
5. **Two grammars or one?** The software template says `<work-stem>` where the
   analysis template says `<stage-stem>`. Is that one grammar with two names
   (fine) or two conventions (a smell)? The lint evaluates the stem predicate
   only where a stage directory exists to check against.
6. **What does the lint check flag that it should not, and miss that it should
   catch?** Be concrete about predicates. Do not propose one needing semantic
   judgement.
7. **The strongest argument that this is still too much.** Make it properly. If
   the right answer is "one `session.md` per stage and nothing else", say so.

## Deliverable

1. A keep/delete verdict per location, each with the failure it prevents.
2. Your answers to 2–5, with a recommendation, not a menu.
3. Concrete lint changes (add / remove / loosen).
4. What to delete from ADR-D10 and from the three templates.
5. Where you disagree with this brief.

Distinguish what you verified by reading from what you inferred.
