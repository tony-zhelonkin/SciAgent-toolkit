# Commands

Canonical slash-command definitions. The resolver walks this tree recursively
and plants flat symlinks under `.claude/commands/`, so subfolders here are
purely organizational. Basename uniqueness across the whole tree is enforced
by `tests/test_no_duplicate_basenames.sh`.

| Command | Purpose | Bound by |
|---------|---------|----------|
| `/commit` | Create a git commit with a curated message | every role (top level) |
| `/map` | Map a feature into a design surface | architect |
| `/design` | Produce a design doc from a map | architect |
| `/diagram` | Render an architecture diagram | architect |
| `/plan` | Build a phased implementation plan | architect |
| `/implement` | Execute a phase from a plan | architect |
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

- `architect/` — the 14 architect-pipeline commands

`commit.md` lives at the top level because it is bound by every role.
