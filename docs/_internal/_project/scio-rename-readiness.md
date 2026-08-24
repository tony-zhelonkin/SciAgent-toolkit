# The scio rename — what is done and what the fleet still has to do

**Date:** 2026-08-24 · **Scope:** ADR-D9

## Premise it started from

ADR-D9 decided the rename in August and sequenced it into the fleet sweep. The
sweep is deferred until scio is parked, so the rename looked blocked. It was not:
what actually blocked it was that the rename was a **flag day**.

## What was wrong

Six places in the mechanism hardcoded `SciAgent-toolkit`, and two shipped
executable templates did too. The day a project moved to `01_modules/scio/`:

- `scio link` would not find the project's own toolkit
- `scio lint --check freshness` would stop checking the pin
- the pre-commit hook would stop linting
- the caption sweep's fallback would stop running

All four silently. That is why the sweep had to be atomic, and why it kept
looking bigger than it was.

## What is done

- `lib/scio/common.sh::_SCIO_TOOLKIT_DIRS=(scio SciAgent-toolkit)` — one list,
  sourced first by every verb. `bin/scio` decides module order; modules do not
  source each other.
- Discovery is staged: declared in `.gitmodules`, then `01_modules/`, then one
  level down. **Evidence outranks the name**; the name only breaks ties inside a
  stage. A candidate must carry `bin/scio` or `craft.yaml` to count.
- Freshness requires a candidate with its own `.git`, because `git -C` walks
  upward and would otherwise compare the project's repository and report nothing.
- The two templates find the toolkit by globbing `01_modules/*/bin/scio`, so they
  are name-independent and cannot desynchronize from the list.
- The caption hook called `validate --check captions`. That verb was deleted with
  validate.sh; its fallback path could never have worked. It calls `lint`.
- Shipped prose, `README.md`, and a new project's `analysis_config.yaml` name
  `01_modules/scio/`.

**A project may sit at either name. The sweep no longer has to be atomic.**

## What the fleet still has to do

1. `git mv 01_modules/SciAgent-toolkit 01_modules/scio` per consumer, plus the
   `.gitmodules` path edit. ~25 copies, four vendor-dir spellings live
   (`01_modules`, `01_Modules`, `01_scripts`, `01_Scripts`).
2. Rename the GitHub repository to `scio`. **Owner action** — it needs GitHub
   access, and GitHub redirects the old URL so existing fetches keep working.
3. Rename the submodule path in scbio-docker (`toolkits/SciAgent-toolkit`). This
   edits the parent's `.gitmodules` and index; it is not a pin bump.
4. `.git` gitdir pointers follow the directory and need no separate edit.

## Deliberately not renamed

- `docs/extractor-import-fidelity/**` and `docs/implementation-kickoff.md` carry
  `repo_root:` frontmatter recording where that work happened. A record of the
  past does not move.
- Tests that build the legacy name are fixtures, not claims.
- `docs/_internal/**` is memory.

## Trap worth carrying

Editing a hook template obliges `tools/gen-template-provenance.sh`. Without the
rehash, every consumer holding the new bytes is read as user-authored and never
refreshed again. The manifest is append-only, so older hashes stay recognized.
