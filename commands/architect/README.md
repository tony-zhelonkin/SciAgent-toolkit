# architect/ commands

The slash commands that drive the architect role's design-before-code
pipeline. They cover the standard cadence (`/map` → `/review` → `/synthesize`
→ `/design` → `/plan` → `/implement` → `/verify`) plus meta-level variants
for changes to the pipeline itself, the `/diagram` rendering helper, the
`/status` phase-doc renderer, the `/architect` consistency gate, and the
retrospective architecture-treemap audit stack.

For the canonical sequence and worked examples see
`../../docs/workflows/architect/00-quickstart.md`. For the rationale behind
the pipeline shape see `../../docs/workflows/architect/01-architecture.md`.

## Commands at a glance

**Per-feature pipeline**

| Command | What it does |
|---|---|
| `/map <slug>` | Read-only codebase cartography → `docs/{slug}/map.md` |
| `/review <slug> --as <spec>` | N reviewer lenses in parallel → `review/*.md` |
| `/synthesize <slug>` | Collapse ≥2 reviews into a consensus → `synthesis.md` |
| `/design <slug>` | Staged drafting + architect gate → `design/*.md` |
| `/architect <slug>` | Standalone architect consistency/completeness gate |
| `/plan <slug>` | Phased decomposition → `plan/phase-NN.md` |
| `/implement <slug> <N>` | Implement one phase (hand-paced; `--auto` for end-to-end) |
| `/verify <slug>` | Mechanical drift check: plan + design vs code → `verify.md` |
| `/diagram <slug>` | Collect Mermaid from `design/` → `diagrams.md` (read-only) |

**Portfolio / meta**

| Command | What it does |
|---|---|
| `/status [slug,...]` | One-page portfolio-state snapshot (read-only, <30s) |
| `/meta-map` | Cross-feature portfolio inventory → `_meta/map.md` |
| `/meta-design` | Inter-feature ADRs (MADRs) → `_meta/design.md` |
| `/meta-apply` | Propagate accepted MADRs to in-flight designs in parallel |
| `/meta-plan` | Sequence portfolio phases → `_meta/plan.md` |

**Architecture-treemap audit stack** (retrospective; probationary, prune review 2026-11-20)

| Command | What it does |
|---|---|
| `/components-extract <repo>` | Deterministic AST/git/metrics substrate → substrate `components.json` |
| `/audit-slice <concern...>` | Dispatch N concern-driven slicers in parallel (user-authored trigger; off-the-autopilot) |
| `/synthesize-audit` | Integrate slices + substrate → full renderable `components.json` + decision menu |
| `/architecture-treemap [path]` | Render a self-contained `treemap.html` from `components.json` (deterministic; validates first) |
