# Phase 02 — the `internal-memory` check

**Repo:** scio · **Blocked by:** phase 01 (same file: `lib/scio/lint.sh`)
**Read first:** `00_INDEX.md` §2, `CONSULT_enforcement.report.md` §2

## Goal

One new opt-in lint check that is the *only* enforcement for the memory
skeleton. No hook accompanies it — that is the decision, not an omission.

## Contract

`scio lint --check internal-memory [--strict]`

- **Absent `docs/_internal/` → clean no-op, zero output.** Non-negotiable: 10 of
  24 projects lack it and absence is legitimate. Use the same absent-subject
  guard shape as `docs-layout` and `results-layout`.
- When present, the permitted immediate children are `_project/` and stage
  stems matching a real `02_analysis/stages/NN_<stem>.*` (accept the
  `02_analysis/scripts/` spelling too, for the pre-rename window — see
  `test_lint_stage_thinness.sh` test 6 for the precedent).
- Every stage directory present must hold a non-empty `session.md` **or** at
  least one non-empty `reasoning/*.md`. A directory with neither is the
  empty-scaffold failure this check exists to catch.
- Flag: `.gitkeep` anywhere under `docs/_internal/`; the retired flat
  namespaces `handoffs/`, `sessions/`, `plans/`, `reports/`, `research/` as
  immediate children; dated session piles (more than one `*session*.md` in one
  directory).
- Flag non-memory payloads: caches, virtualenvs, checkpoints, bytecode, logs,
  parquet. The field precedent for why: one project's `docs/_internal/` is
  568 MB because a `.venv` lives in it (3,813 files, 91.5% of the tree);
  another is 231 MB around a 196 MB model checkpoint.
- Exempt a nested `docs/_internal/.git/` — that topology is recommended in
  ADR-D10 and must not be reported as junk.
- WARN by default, exit 1 under `--strict`, matching every sibling check.

Register in the dispatch `case`, the `--check` validity list, the `-h` usage
text, and the file-header check list. It is **opt-in like `toolkit`** — do NOT
add it to `all`. Rationale: it would fire in every consumer before any project
has adopted the skeleton, and the always-on `all` set must stay quiet on
legitimate absence.

Selectable names go 11 → 12; `all` stays at 10.

## What it cannot do — state this in the header comment

It cannot judge whether reasoning is sound, whether anything important went
unrecorded, or whether a note is current. It cannot see work that never left
`/tmp`. Claiming otherwise would repeat the defect this plan removes.

## Verify

```bash
bash tests/run-all.sh
bin/scio lint --check internal-memory --help   # name appears in usage
bin/scio lint --check bogus                    # error lists the new name
```

New `tests/test_lint_internal_memory.sh` covering, at minimum: absent tree is
silent; a valid `_project/` + one stage dir passes; a stage dir with neither
`session.md` nor `reasoning/*.md` fails; a `.gitkeep` fails; `handoffs/` as an
immediate child fails; a stage stem with no matching stage file fails; a
committed `.venv`/`.parquet` fails; a nested `.git/` does **not** fail;
`--strict` promotes; non-strict exits 0.

## Do not

- Do not add a hook of any kind.
- Do not add this check to `all`.
- Do not create `docs/_internal/` in any template.
- Do not require any hand-maintained metadata field — the owner rejected that
  approach explicitly (`BRIEF_consult.md` §2).
