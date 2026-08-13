# Commands

Canonical slash-command definitions live as flat files. `scio link`
exposes the directory at `.claude/commands` and `.agents/commands`. Source
filenames match their consumer-visible names. Basename uniqueness is enforced
by `tests/test_no_duplicate_basenames.sh`.

| Command | Purpose | Bound by |
|---------|---------|----------|
| `/commit` | Create a git commit with a curated message | every role |
| `/decompose` | Multi-agent, multi-pass planning campaign over a large multi-artifact problem (planning only) | all roles |
| `/implement` | Multi-agent, multi-pass implementation campaign driven by a SOLIDIFIED plan — Opus conductor + per-phase implementer+reviewer pairs + seam windows (continues `/decompose`) | all roles |
The architect workflow is routed through `skills/architecture-first-dev/SKILL.md`, which loads its complete stage specifications from that skill's `references/` directory.

Lane provenance and composition live in `docs/lanes.md`.
