# `architect-mentor` vs `cs101` — drift investigation

## Question

The toolkit ships exactly one system-prompt file today —
`system-prompts/cs101.md` — whose frontmatter declares
`name: architect-mentor`. The `architect` role pins
`output_style: architect-mentor`. Downstream projects had dangling
symlinks to `cs101_v0.1.md` and `cs101_v0.2.md`.

Were `cs101_v0.1`, `cs101_v0.2`, and `architect-mentor` ever distinct
styles, or are they renamed/superseded versions of one another?

## Findings

**They were distinct styles in the initial prototype, then collapsed
during the sciagent refactor to a single style.**

### Timeline

| Commit | Subject | Effect on system prompts |
|---|---|---|
| `8023cb7` | Prototype for a compbio-oriented DevOps agentic harness | Added `system-prompts/cs101_v0.1.md` (frontmatter `name: cs-101 style`, 110 lines) and `system-prompts/cs101_v0.2.md` (frontmatter `name: architect-mentor`, 141 lines). Two distinct files with different `name:` and different content. |
| `ac0e7aa` | restructure dirs for sciagent refactor | Deleted `cs101_v0.1.md`; renamed `cs101_v0.2.md` → `output-styles/cs101.md`. Result: one file remains, named `architect-mentor`. |
| `4c80959` | cs101: add Architectural Axioms section | Extended `output-styles/cs101.md` with a new section. Same identity. |
| `8e15ccb` | scrub personal name, profession, and email from docs and skills | Edited `output-styles/cs101.md` (46 lines changed). Same identity. |
| `8c3380a` | refactor: rename output-styles/ to system-prompts/ as source dir | Moved `output-styles/cs101.md` → `system-prompts/cs101.md`. Same identity. |

### Identity of each historic style

- **`cs101_v0.1` (`name: cs-101 style`)** — a Socratic concept-first
  teaching style for a computational biologist learning CS fundamentals.
  Body title: "Mentor style Instructions". Genuinely different from v0.2
  in tone and structure. **No longer exists in the tree.**
- **`cs101_v0.2` (`name: architect-mentor`)** — a fuller rewrite
  targeted at a clinician-scientist building architectural fluency.
  Body title: "Architect-Mentor Style". Cites the "three interaction
  patterns" learning-literature framing. **This is what today's
  `system-prompts/cs101.md` descends from.**

### Today

`system-prompts/cs101.md` is the v0.2 lineage plus the
Architectural-Axioms extension and the scrub. The `cs101` filename is
historical — the `name:` field (`architect-mentor`) is the canonical
identifier. The two dangling symlinks `.claude/output-styles/cs101_v0.1.md`
and `cs101_v0.2.md` were never updated when `ac0e7aa` restructured the
toolkit; they pointed at files that no longer exist.

## Implications

1. **Single coherent style today.** No salvage work needed for v0.1 —
   it was a deliberately different prototype that didn't survive the
   refactor. If a Socratic-teacher style is wanted again, recreate it as
   a new file with a fresh `name:` (e.g., `cs101-socratic`).
2. **Frontmatter `name:` is the authority.** The resolver added in
   commit `7bcba8f` (`feat(activate): resolve output_style by
   frontmatter, validate before mutation`) treats the frontmatter as the
   single source of truth — filenames are now incidental.
3. **`cs101` as a filename should probably be retired.** Once a second
   system-prompt file exists, the name is misleading (the file isn't
   "v0.1" of anything anymore). A rename like
   `system-prompts/architect-mentor.md` would match the frontmatter
   identity and reduce future drift. Deferred — the new resolver makes
   this a no-op refactor whenever convenient.
