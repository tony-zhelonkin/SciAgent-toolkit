# architect/ agents

The twelve agents that the `architect` role binds. They implement the
design-before-code pipeline: a feature is *mapped* into a design surface,
*reviewed* by domain specialists, the verdicts are *synthesized*, and a
*feature-reviser* threads the consensus back into the map. `meta-architect`,
`status-reporter`, and `architect` (the consistency gate) coordinate phase
transitions and emit phase docs. See
`../../docs/workflows/architect/00-quickstart.md` for the canonical cadence.

| Agent | Role in the pipeline |
|-------|----------------------|
| `mapper` | Builds the initial feature map |
| `architect` | Consistency gate over design output |
| `meta-architect` | Coordinates meta-level architecture decisions |
| `synth` | Synthesizes reviewer verdicts |
| `feature-reviser` | Revises the map per synthesized review |
| `status-reporter` | Renders phase-doc status |
| `bioinf` | Bioinformatics design reviewer |
| `wetlab` | Wet-lab feasibility reviewer |
| `graphic` | Visualization design reviewer |
| `stat` | Statistical design reviewer |
| `divergent` | Adversarial reviewer |
| `ml` | Machine-learning design reviewer |
