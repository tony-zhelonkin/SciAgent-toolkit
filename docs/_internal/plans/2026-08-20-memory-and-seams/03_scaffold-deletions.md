# Phase 03 — stop shipping empty directories

**Repo:** scio · **Parallel-safe** · **Read first:** `00_INDEX.md` §2 (example-first)

## Goal

The scaffold taught agents the wrong thing. Remove it. A directory arrives with
its first real file or it does not arrive.

## Evidence this is not a preference

All five STING child projects ship `handoffs/` containing only `.gitkeep`. In
**every one**, the real handoffs were written to `sessions/` instead. Each child
carries ~10 scaffold files before any content exists. Ten spellings of "handoff"
exist across the fleet — because the sanctioned location was empty, so every
agent invented its own. Detail: `FINDINGS_field.md` §3, §4.

An empty directory is not an example. It is an unsupported claim, and agents
discounted it correctly.

## Changes

In `templates/project/**`, delete from the scaffolded tree:

- every `.gitkeep` under a `docs/_internal/` path
- the default category directories `handoffs/`, `sessions/`, `plans/`,
  `reports/`, `research/`
- any per-directory `README.md` that only describes a naming convention for a
  directory with no content

**Keep** `docs/_internal/README.md` if and only if it states the local contract
in a few lines. A short true contract is worth keeping; a naming spec for empty
directories is not.

Audit what remains against one question: *would an agent reading this tree see
an instance to copy, or a promise?* Delete the promises.

## Verify

```bash
bash tests/run-all.sh
bash tests/test_mount_layout_identity.sh   # source shape == mount shape invariant
bash tests/test_plan_templates.sh
bash tests/test_template_provenance.sh     # if any MANAGED template moved

# A fresh scaffold must contain no empty directory and no .gitkeep:
find <scratch-project> -type d -empty
find <scratch-project> -name .gitkeep
```

Both `find` calls must print nothing.

## Do not

- Do not touch `lib/scio/lint.sh` — phases 01 and 02 own that file.
- Do not remove `templates/skill/` or anything under the three mounted
  categories: a directory symlink cannot filter, so source shape must equal
  mount shape (`test_mount_layout_identity.sh` enforces this).
- Do not delete anything from an existing consumer project. This phase changes
  what new scaffolds produce, nothing already in the field.
