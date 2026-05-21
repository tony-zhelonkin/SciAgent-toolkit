# Agents

Canonical sub-agent definitions. The resolver walks this tree recursively and
binds agents into `.claude/agents/` and `.agents/agents/` as flat symlinks, so
the subfolder layout here is purely organizational — it groups agents by the
role family they serve. Each `<name>.md` carries YAML frontmatter (`name`,
`description`, `model`, `color`); basename uniqueness across the whole tree is
enforced by `tests/test_no_duplicate_basenames.sh`.

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
| `graphic` | Visualization design reviewer | `/review --as graphic` |
| `handoff` | Session handoff doc generator | base role, end-of-session |
| `insight-explorer` | Data-file exploration with skepticism | base role, on-demand |
| `mapper` | Map a feature into design surface | `/map` |
| `meta-architect` | Meta-level architecture coordinator | architect pipeline |
| `ml` | Machine-learning design reviewer | `/review --as ml` |
| `stat` | Statistical design reviewer | `/review --as stat` |
| `status-reporter` | Phase-doc status renderer | `/status` |
| `synth` | Synthesize review verdicts | `/synthesize` |
| `wetlab` | Wet-lab feasibility reviewer | `/review --as wetlab` |

Subfolders:

- `architect/` — agents that serve the architect role (12)
- `analysis-base/` — agents that serve the base bioinformatics role (7)

See `architect/README.md` and `analysis-base/README.md` for one-paragraph
orientations on each family.
