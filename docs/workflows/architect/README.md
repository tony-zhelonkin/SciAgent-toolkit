# architect role — Workflow documentation

An architecture-first software development harness for Claude Code. Composable expert review panel + staged pipeline + human gates at every boundary.

---

## TL;DR

You're building a non-trivial feature. You want to think about structure, behavior, and decisions before writing code — and you want multiple expert lenses (statistical, ML, wet-lab, UX, contrarian) applied before design. This harness gives you six slash commands that walk through that process, each stage producing a markdown artifact the next stage reads.

```
/map <feature>                   → docs/{feature}/map.md
/review <feature> --as <spec>    → docs/{feature}/review/<reviewer>.md (parallel)
/synthesize <feature>            → docs/{feature}/synthesis.md (if ≥2 reviewers)
/design <feature>                → docs/{feature}/design/*.md (gated by architect)
/plan <feature>                  → docs/{feature}/plan/phase-NN.md
/implement <feature> <phase>     → code, one phase at a time
```

**Meta-layer (optional, for ≥2 features in flight):**

```
/meta-map                        → docs/_meta/map.md (cross-feature inventory + touch-points)
/meta-design                     → docs/_meta/design.md (inter-feature ADRs — MADRs)
/meta-plan                       → docs/_meta/plan.md (sequences per-feature phases across portfolio)
```

---

## Document index

| File | View | What it answers |
|------|------|------------------|
| [00-quickstart.md](./00-quickstart.md) | **Start here** | Cadence, decision tree, per-stage tips, prototypical meta-design walkthrough |
| [01-architecture.md](./01-architecture.md) | Structural | What files exist, how they connect, what the role YAML binds together |
| [02-behavior.md](./02-behavior.md) | Behavioral | What each stage does, data flow between stages, dispatch patterns |
| [03-decisions.md](./03-decisions.md) | Decisions | Why this shape — the token-efficiency contract, composable panel design, naming conventions, tradeoffs |
| [04-extending.md](./04-extending.md) | How-to | Adding a new reviewer, a new command, or adapting the role for a different project |
| [05-meta-approach.md](./05-meta-approach.md) | How-to | The meta-layer — when, how, and why to invoke `/meta-map` / `/meta-design` / `/meta-plan` for cross-feature work |

---

## Quick start

### Activate the role

From the toolkit root:

```bash
scripts/activate-role.sh architect --project-dir /path/to/your/project
```

This symlinks 10 agents, 1 skill, and 9 commands into `.claude/` of the target project. Reload Claude Code to pick them up.

### Run the pipeline on a feature

```
/map umap-fix "Audit UMAP primitive: neighbour-preservation, entity mixing, contrast dependence"
/review umap-fix --as all
/synthesize umap-fix
/design umap-fix
/plan umap-fix
/implement umap-fix 1
```

Each command stops at a human gate. Review the artifact, redirect if needed, advance when satisfied.

### Invocation shapes for `/review`

```
/review <feature> --as all                           — all 6 reviewers in parallel
/review <feature> --as all --but stat,divergent      — all minus listed
/review <feature> --as bioinf,ml                     — specific subset
/review <feature>                                    — interactive; asks which to run
```

---

## What this is not

- **Not for one-line fixes** — use the generic Claude Code flow
- **Not for throwaway scripts** — the gate overhead dominates small work
- **Not a bioinformatics role** — intentionally minimal; use `base` or `pathway-signature-agent` alongside if you need anndata/scanpy/scvi skills

---

## At a glance — what was built

| Artifact | Location | Purpose |
|----------|----------|---------|
| Role YAML | `roles/architect.yaml` | Manifest binding 10 agents + 1 skill + 9 commands |
| Reviewer agents | `agents/{bioinf,wetlab,graphic,stat,divergent,ml}.md` | 6 reviewer perspectives |
| Pipeline agents | `agents/{mapper,architect,synth}.md` | 3 per-feature pipeline agents |
| Meta agent | `agents/meta-architect.md` | Cross-feature orchestrator (writes only `docs/_meta/*`) |
| Per-feature commands | `commands/{map,review,synthesize,design,plan,implement}.md` | The 6 per-feature slash-command templates |
| Meta commands | `commands/{meta-map,meta-design,meta-plan}.md` | The 3 cross-feature slash-command templates |
| Skill | `skills/architecture-first-dev/SKILL.md` | Methodology reminder loaded with the role |
| Activation | `scripts/activate-role.sh` (extended) | Now symlinks `commands:` list in addition to agents/skills |
| Runtime per-feature artifacts | `docs/{feature}/{map,review,synthesis,design,plan}/` | Where the pipeline writes per feature |
| Runtime meta artifacts | `docs/_meta/{map,design,plan,deferred}.md` | Where the meta-layer writes (sibling of per-feature dirs) |

See [01-architecture.md](./01-architecture.md) for the complete component map.
