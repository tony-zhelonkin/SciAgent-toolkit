# Plan INDEX — stop claiming what the mechanism cannot deliver

**Date:** 2026-08-20 · **Scope:** IN — remove the toolkit's false claims, add the
one check that replaces them, restore the scio/dev-env seam, write ADR-D10 /
OUT — the fleet sweep, the ADR-D9 rename, memory remotes and hub syncing, any
new hook, mass-scaffolding `docs/_internal/` into the 10 projects that lack it

**Rule:** one phase == one bounded brief == one implementer. Phases are
numbered by slug, not reading order. Read §4 before fanning out: two phases
touch the same file and must serialize.

Evidence for every claim below: `docs/_internal/research/2026-08-20-internal-skeleton/`
— four read-only field inspections across 24 provc-managed projects
(`inspect-*.report.md`), a digest (`FINDINGS_field.md`), and two consults
(`CONSULT_skeleton.report.md` on layout, `CONSULT_enforcement.report.md` on
mechanism and ownership). Do not re-derive; cite.

---

## 1. The one finding this plan acts on

Every defect found on 2026-08-20 is the same defect: **an instruction naming
something the mechanism does not guarantee.**

| The toolkit says | The mechanism does | Verified at |
|---|---|---|
| durable memory lives in tracked files, in `docs/_internal/` | `link` writes `docs/_internal/` into `SCIO:GITIGNORE` | `lib/scio/link.sh:278` |
| `_scratch/` is "the only sanctioned throwaway zone" | no project has one at its root; one has an empty one | `FINDINGS_field.md` §6 |
| "`docs/_internal/` missing — run: `scio link`" | `link` never creates that directory | `lib/scio/lint.sh:604` |
| AGENTS.md: read the skill | activation mounted 60 of 86 skills, silently | 14761-DM, measured |

The last one cost a live investigation in 14761-DM: `Error: Unknown skill:
delegate-cli` while the skill sat on disk in that project's own pinned
checkout, one of **26 unreachable of 86**. The current `link` deletes that
failure class (membership *is* directory membership), so it needs no fix here —
it needs delivery, which is the blocked fleet sweep.

**The review question this plan installs, for every future claim: what enforces
this, and can it actually deliver?**

## 2. What we decided, compactly

- **Layout.** Memory mirrors the analysis: `02_analysis/stages/30_grn.R` →
  `03_results/30_grn/` → `docs/_internal/30_grn/`. `_project/` holds memory
  that spans or precedes stages. **Two durable forms only:** `reasoning/<topic>.md`
  and `session.md`, the latter updated in place with git history as the archive.
- **No mirror for `02_analysis/helpers/`.** Category error: a helper serves
  several stages, so its rationale files with the consuming stage or under
  `_project/reasoning/`.
- **No new hook.** Not Claude, not Codex, not git `pre-commit`. Ranked
  surfaces: populated filesystem > `scio lint` > git hook > harness hooks.
  A `pre-commit` calling lint stays available as a later escalation *if*
  evidence shows agents skip lint.
- **Example-first.** Structure propagates when the tree holds an authentic
  instance whose form can be copied. `.gitkeep` is not an example — it is an
  unsupported claim, and agents correctly ignored it in every STING child
  while writing real handoffs to `sessions/`. So: never scaffold empty; create
  a directory with its first real file.
- **The seam.** scio owns the project's path grammar, examples, always-on
  router text and lint predicates. dev-env owns user-global habits. Test: *if
  it must follow the repo across humans and machines it is scio's; if it must
  follow this human across unrelated repos it is dev-env's.* Corollary that
  settles the smearing: **dev-env may name a path for a tool (watching,
  mounting, credentials); it may not name a path for an agent.**
- **Topology** (recorded in ADR-D10, not implemented here): `docs/_internal/`
  as a nested git repo in place — the parent already ignores it, so nesting
  needs no wiring — plus a **parent-tracked** `docs/internal-memory.md` naming
  where the memory lives, so a nested repo is declared rather than hidden.
  14761-DM and 14782-DM already invented this independently (126 and 210
  commits, 195 and 708 tracked files, **no remote on either**).

## 3. Phase table

| # | Slug | What it ships | Repo | Parallel? |
|---|---|---|---|---|
| 01 | `remove-false-claims` | delete the 3 unbacked claims in lint/craft/hook text | scio | after none; **blocks 02** |
| 02 | `internal-memory-lint` | the `internal-memory` check that replaces them | scio | **after 01** |
| 03 | `scaffold-deletions` | remove `.gitkeep` + default category dirs from templates | scio | yes |
| 04 | `dev-env-seam` | strip project paths from the user-global template | scio | yes |
| 05 | `scbio-docker-stale-verb` | 7 docs citing a deleted verb | scbio-docker | yes |
| 06 | `adr-d10` | the decision record for all of the above | scio | yes |
| 07 | `delegation-seam` | cite scbio-docker for the codex sandbox cause | scio | yes |

Suggested fan-out: **{01→02} ∥ 03 ∥ 04 ∥ 05 ∥ 06 ∥ 07**, six lanes, the first
serial internally. 06 is the decision record — if an implementer disagrees with
a deletion, 06 is where the disagreement is recorded and resolved, not the code.

Phase 07 also carries an **owner decision** (the `seccomp=unconfined` posture)
that no implementer may take. See `07_delegation-seam.md`.

## 4. Conflict map — read before delegating

- **01 and 02 both edit `lib/scio/lint.sh`.** 01 deletes a predicate, 02 adds
  a check plus dispatch, help text and header docs. Serialize them, or give
  both to one implementer. Do not run them in parallel worktrees and merge.
- 03 edits `templates/project/**` only. No overlap with 01 (which edits
  `templates/project/_common/.claude/hooks/no_ephemeral.sh.template` — a
  different file, but same tree: confirm no `_common/docs/_internal/` scaffold
  is touched by both).
- 04 edits `docs/_internal/handoff/dev-env/global-AGENTS.md.template` only.
- 05 is a different repository. It has its own pin drift (`3ab3768` recorded vs
  checked out) — **do not bump the submodule pin** as part of 05.
- 06 edits `docs/proposals/2026-08-11-offline-distribution/50_ADRs.md` only.

## 5. Global gates — every phase, before reporting done

```bash
bash tests/run-all.sh                      # 62 passing, 0 failing at f234b2c
bin/scio lint --check toolkit --strict     # clean
```

Plus, for any phase touching `craft.yaml`:

```bash
# The rendered body must stay byte-identical unless the change is INTENDED to
# alter it — every consumer stores a SHA1 of this block and a change reads as
# drift in all 24 on re-pin. Phase 01 DOES intend to change it; diff and read.
git show HEAD:craft.yaml > /tmp/before.yaml
SCIO_TOOLKIT=$PWD bash -c '. lib/scio/block.sh; . lib/scio/craft.sh; _craft_render_body'
```

The CRAFT body is 17 lines against a `max_lines: 25` cap enforced by
`lint --check toolkit` since `287e8e4`. Phase 01 removes lines; nothing here may
add net lines.

## 6. Out of scope — do not do these

- **The fleet sweep** (#30). Blocked on the owner's review of `106f59f..dev`.
  It is what actually delivers the 26 unreachable skills to 14761-DM, but it is
  not this plan.
- **The ADR-D9 rename** (repo → `scio`, vendor path → `01_modules/scio/`).
  Owner's call, and it would make phase 05's wording wrong twice.
- **Memory remotes / hub syncing.** The owner explicitly deferred this:
  durability of *structure and intent* is the subject; where bytes get pushed
  is not.
- **Touching the 8,317 existing files** across 24 projects. No relocation
  exercise. Existing memory is handled opportunistically when a project next
  becomes active.
- **Creating `docs/_internal/` anywhere.** 10 of 24 projects lack it and that
  is fine; absence is the majority state and phase 01 exists precisely to stop
  the toolkit pressuring them.
- **Any new hook.**

## 7. Blocked on the owner, carried forward

1. Review `106f59f..dev` — now 4 commits, and the owner intends to tweak
   wording during review (see `nvim` note in the session handoff).
2. Commit 14761-DM (was blocked on a live `codex exec` writing
   `02_analysis/helpers/grn_viz/`; re-check).
3. ADR-D9 rename: run before or after the sweep.
4. Task #36 — whether the CRAFT block gets a token-weight budget. The cap
   binds shape (17/25 lines) but not cost (4,658 chars, longest bullet 693).
5. The `seccomp=unconfined` posture for codex delegation inside the containers
   (phase 07). Recommendation is to enable it; it is a real security loosening
   on a shared machine, so it is not an implementer's call.

## 8. State at plan time

`dev` = `f234b2c`, **5 ahead** of `hub/dev` = `origin/dev` = `106f59f`.
`main` = `c83dfe2` both remotes. 62 tests passing. `lib/scio` 2,619 LOC.
Catalog: 85 skills, 21 agents, 7 commands. 9 ADRs (D10 will be the tenth).
Fleet: 24 real toolkit copies — 20 at `5e5347e` (52 behind; 15 fleet-managed,
5 frozen by owner decision 2026-08-18), 1 at `106f59f`, 1 canonical, 1 PanSci
pinned, 1 excluded.
