# Consultation — choosing at publication time whether project memory goes public

Read-only. Produce a written recommendation. **Modify nothing.**

## The requirement, in the owner's words

> Ideally we need a mechanism where, when we make the repo public, we decide and
> can choose whether to make its `_internal` public as well, with a durable
> link — be it a submodule or some other less churnful topology.

So: the memory must be **separable at publication**, the choice must be **the
owner's, made late**, and if the answer is yes there must be a **durable link**
from the public repository to the memory. If the answer is no, the public
repository must not leak it.

## What exists today, verified

- `docs/_internal/` is written into the project's `.gitignore` by
  `scio link` (a managed `SCIO:GITIGNORE` block). The parent repository does not
  track it and cannot see it.
- Measured across five sampled projects: **8,317 files, 2 tracked.**
- **Two projects independently made `docs/_internal/` a nested git repository** —
  126 and 210 commits, 195 and 708 tracked files — and **neither has a remote**.
  So they have history and no reachability.
- ADR-D10 decision 6 records that topology (nested repo in place, plus a
  parent-tracked `docs/internal-memory.md` naming it) and deliberately does not
  implement it. Read it:
  `docs/proposals/2026-08-11-offline-distribution/50_ADRs.md`.
- The grammar is `docs/_internal/<stage-stem>/session.md` plus flat topic notes,
  and `_project/` for what spans stages, with plans at
  `_project/plans/<date-slug>/`. `session.md` is **updated in place**.
- Enforcement is one opt-in check, `scio lint --check internal-memory`
  (`lib/scio/lint.sh`, search `_lint_check_internal_memory`). Read it.
- Full field evidence:
  `docs/_internal/research/2026-08-20-internal-skeleton/FINDINGS_field.md`.
  Structural verdict just accepted:
  `docs/_internal/research/2026-08-21-internal-minimalism/CONSULT_minimalism.report.md`.

## Two objections already on the record

**Churn.** A submodule was rejected in ADR-D10 because durability would depend on
pin bumps, and 20 of 24 fleet copies currently sit at one stale pin. Parent
tracking was rejected because agent write volume would flood the human's commit
log and restore the publication hazard.

**The in-place rule needs an archive.** ADR-D10 justifies overwriting
`session.md` with "git history is its archive". That is only true where a
repository exists. In an ignored, non-nested tree an overwrite destroys the prior
state with nothing behind it. The owner has just ruled that the absolute claim is
deleted and that `lint --check internal-memory` should instead **observe** the
condition: a populated memory tree with no history mechanism gets a finding.
`link` creates nothing. Take that ruling as settled input, not as a question.

## Questions

1. **Enumerate the topologies** that could satisfy "decide late, durable link, no
   churn in the parent's log", and rank them. At minimum consider: nested repo
   with a remote added at publication; submodule whose pin is bumped *only* at
   publication; an orphan branch in the same repository; a separate repository
   named by a tracked pointer file; `git subtree`; a snapshot committed into the
   public repo at release. For each: what it costs during daily work, what it
   costs at publication, and how it fails.
2. **Does "bump the pin only at publication" dissolve the churn objection?** That
   was the objection's whole basis. If it does, say what breaks instead — a
   submodule pointing at a repository that does not exist yet, a pin that is
   stale by construction between publications, `git clone` behaviour for someone
   without access.
3. **The orphan-branch option specifically.** `git checkout --orphan memory` in
   the same repository: one repo, separate history, the parent's default log
   unpolluted. But publishing the repository publishes every branch, which is the
   wrong default when the answer is "no". Is it salvageable, and how would the
   "no" case be enforced rather than remembered?
4. **What is the minimum mechanical support scio should provide?** Constraints,
   all settled: **no new hook** of any kind; **bash only**, no `jq`/`yq`/
   `python3`; **never scaffold empty**; **no hand-maintained metadata** a human
   must keep accurate; `lint` is the enforcement surface; and **do not name a
   location nothing creates**. A documented procedure plus one lint predicate may
   well be the right answer — say so if it is. If you propose anything that
   writes, say which verb owns it and why that verb.
5. **Redaction is the real hazard.** Publishing reasoning notes can leak
   credentials, absolute lab paths, unpublished data, collaborator names, and
   patient-adjacent identifiers. What can be checked *mechanically* before a
   publication step, and what is irreducibly a human read? Do not propose a check
   that needs semantic judgement.
6. **Recovery.** Someone clones the public parent in three years. What must the
   pointer file contain for the memory to be actually obtainable — and what
   should it say when the answer was "memory stays private"? An unresolvable
   pointer is worse than no pointer.
7. **The lint predicate for the ruling above.** Specify it: how does a bash check
   decide that a populated memory tree "has no history mechanism", without
   `jq`, and without false-positiving on a tree that is legitimately brand new?

## Deliverable

1. Ranked topologies with the failure mode of each (Q1–Q3).
2. One recommendation, stated as what changes in scio and what stays a procedure.
3. The publication-time sequence, concretely, for both answers — yes and no.
4. The redaction map: mechanical vs human (Q5).
5. The pointer file's required content for both answers (Q6).
6. The lint predicate, in enough detail to implement (Q7).
7. Where you disagree with this brief.

Distinguish what you verified by reading from what you inferred. Keep the
recommendation short enough to act on.
