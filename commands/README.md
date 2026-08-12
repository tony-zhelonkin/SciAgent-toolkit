# Commands

Canonical slash-command definitions. The resolver walks this tree recursively
and plants flat symlinks under `.claude/commands/`, so subfolders here are
purely organizational. Basename uniqueness across the whole tree is enforced
by `tests/test_no_duplicate_basenames.sh`.

| Command | Purpose | Bound by |
|---------|---------|----------|
| `/commit` | Create a git commit with a curated message | every role (top level) |
| `/decompose` | Multi-agent, multi-pass planning campaign over a large multi-artifact problem (planning only) | all roles |
| `/implement` | Multi-agent, multi-pass implementation campaign driven by a SOLIDIFIED plan — Opus conductor + per-phase implementer+reviewer pairs + seam windows (continues `/decompose`) | all roles |

Subfolders:

- `decompose/` — the multi-agent planning-campaign command (bound by all roles)
- `implement/` — the multi-agent implementation-campaign command (bound by all roles; continues `/decompose`)

The architect workflow is routed through `skills/architecture-first-dev/SKILL.md`, which loads its complete stage specifications from that skill's `references/` directory.

`commit.md` lives at the top level because it is bound by every role.
