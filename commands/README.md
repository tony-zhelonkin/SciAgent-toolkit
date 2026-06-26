# Commands

Canonical slash-command definitions. The resolver walks this tree recursively
and plants flat symlinks under `.claude/commands/`, so subfolders here are
purely organizational. Basename uniqueness across the whole tree is enforced
by `tests/test_no_duplicate_basenames.sh`.

| Command | Purpose | Bound by |
|---------|---------|----------|
| `/commit` | Create a git commit with a curated message | every role (top level) |
| `/decompose` | Multi-agent, multi-pass planning campaign over a large multi-artifact problem (planning only) | base, architect |
| `/map` | Map a feature into a design surface | architect |
| `/design` | Produce a design doc from a map | architect |
| `/diagram` | Render an architecture diagram | architect |
| `/plan` | Build a phased implementation plan | architect |
| `/implement` | Multi-agent, multi-pass implementation campaign driven by a SOLIDIFIED plan — Opus conductor + per-phase implementer+reviewer pairs + seam windows (continues `/decompose`) | all roles |
| `/review` | Run a domain reviewer over the map | architect |
| `/synthesize` | Synthesize reviewer verdicts | architect |
| `/verify` | Verify phase doc against actual state | architect |
| `/status` | Render phase-doc status | architect |
| `/architect` | Consistency gate over design output | architect |
| `/meta-map` | Meta-level: map a meta-feature | architect |
| `/meta-design` | Meta-level: design a meta-feature | architect |
| `/meta-plan` | Meta-level: plan a meta-feature | architect |
| `/meta-apply` | Meta-level: apply a meta-change | architect |

Subfolders:

- `architect/` — the architect-pipeline commands (map, review, design, plan, verify, …); see `architect/_superseded/` for the former single-phase implement
- `decompose/` — the multi-agent planning-campaign command (bound by all roles)
- `implement/` — the multi-agent implementation-campaign command (bound by all roles; continues `/decompose`)

`commit.md` lives at the top level because it is bound by every role.
