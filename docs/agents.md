# Agents

Canonical sub-agent definitions live as flat `agents/<name>.md` files. The
resolver binds those exact filenames into `.claude/agents/` and
`.agents/agents/`. Each file carries YAML frontmatter (`name`, `description`,
`model`, `color`); basename uniqueness is enforced by
`tests/test_no_duplicate_basenames.sh`.

| Agent | Purpose | Invoked by |
|-------|---------|------------|
| `architect` | Consistency gate for design output | `/architect` (architect role) |
| `bioinf` | Bioinformatics design reviewer | `/review --as bioinf` |
| `bio-interpreter` | Biological mechanism research | base role, on-demand |
| `captions` | Figure-legend generator | base role, on-demand |
| `code-reviewer` | Refactor review pass | base role, on-demand |
| `divergent` | Adversarial reviewer | `/review --as divergent` |
| `doc-curator` | Repo documentation cleanup | base role, on-demand |
| `docs-librarian` | Tool/package documentation lookup | base role, on-demand |
| `feature-reviser` | Revise feature map per review | `/synthesize` follow-up |
| `figure-audit` | Inspect rendered figures against the visual contract | base role, on-demand |
| `graphic` | Visualization design reviewer | `/review --as graphic` |
| `handoff` | Session handoff doc generator | base role, end-of-session |
| `insight-explorer` | Data-file exploration with skepticism | base role, on-demand |
| `mapper` | Map a feature into design surface | `/map` |
| `meta-architect` | Meta-level architecture coordinator | architect pipeline |
| `ml` | Machine-learning design reviewer | `/review --as ml` |
| `stat` | Statistical design reviewer | `/review --as stat` |
| `status-reporter` | Phase-doc status renderer | `/status` |
| `slicer` | Review architectural slices and dependency seams | architect pipeline |
| `synth` | Synthesize review verdicts | `/synthesize` |
| `wetlab` | Wet-lab feasibility reviewer | `/review --as wetlab` |

Lane provenance and composition live in `docs/lanes.md`. See
`docs/agents-architect.md` and `docs/agents-analysis-base.md` for the two agent
family orientations.
