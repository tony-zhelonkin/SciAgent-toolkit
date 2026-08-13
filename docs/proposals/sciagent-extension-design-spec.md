# SciAgent-toolkit — Extension Design Spec

> **Outcome note (superseded).** Superseded by the three-verb toolkit
> architecture documented in [`docs/architecture.md`](../architecture.md).
> Retained as a dated design record; the original body below is unchanged.

**Date**: 2026-05-22
**Status**: Proposed, awaiting decision
**Author**: synthesized from a conversation about Antigravity Science Skills / ERA / Co-Scientist / Robin, against the current `SciAgent-toolkit` architecture (`docs/architecture.md`)

This document is a single-pass design proposal covering ADRs, new skills, sub-agents, slash commands, roles, CLI subcommands, a benchmarking framework. Each section is independently actionable — you can accept ADR-001 without accepting ADR-004. 

---

## 0. Executive summary

Three threads, one toolkit:

1. **Spec hygiene** — your custom SKILL.md frontmatter (`requires`, `contraindications`, `scope`, etc.) is good design but lives in unnamespaced keys, risking future collision with the agentskills.io spec. Fix the namespace; formalize the dependency graph; add a validator. **(ADRs 001, 002, 003, 007.)**

2. **Optimizer overlay** — add an ERA-style tree-search capability as an *overlay role* (`scorable`), not a baseline. Activate when you have a scorable problem; deactivate otherwise. Search space = your existing skill graph. Ablation studies become a single command. **(ADRs 004, 009; skills 2.1–2.5; agents 3.1–3.3; command 4.1; role 5.1.)**

3. **Loop closure** — `workflow-skill-creator` analog via session trace recording (`sciagent record`), and a benchmark harness (`sciagent bench`) modelled on the Antigravity Skills paper's three-tier eval pattern. These compound — the more skills you ship, the more you need both. **(ADRs 005, 006; agents 3.4; commands 4.3, 4.4; CLI subcommands 6.2, 6.3, 6.4.)**

A shipping template overlay, `lab-loop`, demonstrates how to wire comp-bio prediction to a wet-lab assay readout without committing to a specific disease, organ system, or assay format. Users fork it into a project-specific overlay when ready. **(ADR 008; role 5.2.)** Cloud-lab MCP integration is explicitly deferred. **(ADR 010.)**

---

## 1. Architecture Decision Records

ADRs are numbered, dated, statused. Use them as the audit trail for `docs/decisions/NNN-slug.md`.

---

### ADR-001 — Namespace custom metadata under `metadata.sciagent.*`

**Status**: Proposed
**Date**: 2026-05-22

#### Context
The agentskills.io spec defines `metadata:` as a free-form extension point but recommends "reasonably unique key names to avoid accidental conflicts." Current `SciAgent-toolkit` skills use top-level keys inside `metadata`: `requires`, `scope`, `complementary-skills`, `contraindications`, `tier`, etc. If the spec evolves and adopts any of these with different semantics, every skill silently breaks.

#### Decision
Move all custom keys under `metadata.sciagent.*`. Example diff:

```diff
 metadata:
-  scope: orchestrator
-  requires: [foo, bar]
-  contraindications: ["…"]
+  sciagent:
+    scope: orchestrator
+    requires: [foo, bar]
+    contraindications: ["…"]
   skill-author: SciAgent-toolkit
   version: 2.0.0
```

Keep `skill-author`, `version`, `last-reviewed`, `upstream-docs`, `tags`, `category`, `license` at the `metadata.*` level since these have natural conventions in the wider ecosystem.

#### Consequences
- **Positive**: future-proof against spec evolution; explicit toolkit attribution; easy `grep` for sciagent-specific extensions.
- **Negative**: one-time migration across all existing skills (~50+ files). Mitigated by a one-shot codemod.
- **Neutral**: harnesses that ignore `metadata` still work; nothing changes at runtime.

#### Alternatives
- Keep flat (status quo) — risky.
- Use a top-level non-spec key like `x-sciagent:` — diverges further from spec, harder to upstream.

---

### ADR-002 — Formalize `requires` as a hard dependency graph

**Status**: Proposed
**Date**: 2026-05-22

#### Context
`requires:` is currently informational. Orchestrator skills list atomic skills they compose, but nothing enforces that the listed skills exist, no cycles are checked, no version pinning. As the library grows past ~50 skills, drift will produce silently-broken orchestrators.

#### Decision
- `requires:` is a hard contract: every listed skill must exist in `skills/` at activation time.
- `complementary-skills:` stays soft (recommendation only).
- Add semver-compatible version constraints (`requires: [{name: harmonypy, version: ">=0.1.0"}]` or plain `[harmonypy]` resolved to latest).
- Cycle detection is mandatory.
- `sciagent doctor` (ADR-007) enforces the contract.

#### Consequences
- **Positive**: orchestrator integrity guaranteed; safe refactoring of atomics; enables ablation searches (ADR-004).
- **Negative**: stricter authoring discipline; slower skill creation.
- **Neutral**: tier 1 of the dependency graph is shallow today, so migration is cheap now and expensive later.

#### Alternatives
- Soft requires (status quo) — defers the problem.
- External manifest (`skills.lock`) — adds a second source of truth.

---

### ADR-003 — Two-tier scope: `orchestrator` vs `atomic`

**Status**: Proposed (already mostly in place; this formalizes it)
**Date**: 2026-05-22

#### Context
Current skills informally distinguish "orchestrators" (compose other skills) from "atomic" (wrap one tool/library). Some skills currently set `scope:` and `tier:` redundantly. The redundancy creates ambiguity.

#### Decision
- One field: `metadata.sciagent.scope`, values `orchestrator | atomic`. Drop `tier`.
- Orchestrators must have non-empty `requires:`.
- Atomics must have empty `requires:` (no skill composes another skill *and* is itself atomic).
- The two-tier limit is intentional. Three+ tiers create routing pain identical to Claude Code's three-tier resolution problem already flagged in `docs/architecture.md`.

#### Consequences
- **Positive**: clean mental model; clear authoring rule.
- **Negative**: occasional gymnastics when a skill feels mid-tier ("mostly atomic but composes one other skill"). Resolution: it's an orchestrator. Conservatively expand the orchestrator set, not the atomic set.
- **Neutral**: aligns with how Anthropic's official skills repo is implicitly structured.

#### Alternatives
- Three+ tiers — rejected (see Context).
- Single tier with `requires` doing the work — possible but loses authoring guidance.

---

### ADR-004 — Optimizer as overlay role (`scorable`), not baseline

**Status**: Proposed
**Date**: 2026-05-22

#### Context
ERA (Aygün et al., Nature 2026) demonstrates that LLM + tree search over code mutations beats hand-written methods on multiple biology benchmarks (OpenProblems.bio, GIFT-Eval, CDC ensemble for COVID forecasting). The pattern generalizes to any *scorable* task: define a metric, define a baseline, run tree search.

Three implementation choices:
1. Bake into base role — always loaded.
2. Standalone tool, separate from `sciagent`.
3. Overlay role, activated only when needed.
4. 
#### Decision
Implement as an overlay role `scorable` (ADR-008 will mirror this pattern for project-specific overlays). Activate via `sciagent activate <base> scorable`.


Rationale: tree search consumes ~20K tokens per session in agent context (mutator + scorer + selector prompts) and is irrelevant for most interactive work. Making it opt-in keeps the base role clean. The overlay model in `docs/architecture.md` (last-wins shadowing, depth-2 cap) handles this elegantly without new mechanism.

#### Consequences
- **Positive**: zero overhead when not searching; clean activation semantics; reuses existing architecture.
- **Negative**: users must remember to activate before searching; one extra step.
- **Neutral**: composes naturally with any project overlay — `sciagent inject scorable` adds it on top.
#### Alternatives
- Baseline integration — rejected, see Context.
- Separate `eraly` CLI — rejected, fragments tooling.

---

### ADR-005 — `workflow-skill-creator` analog via session trace recording

**Status**: Proposed
**Date**: 2026-05-22

#### Context
Antigravity Skills includes a `workflow-skill-creator` skill that distills a user-handheld interaction into a reusable skill. This is the highest-leverage authoring tool in the bundle: it converts ad-hoc workflows into shipping skills without manual SKILL.md drafting.

Current `SciAgent-toolkit` has `sciagent new skill <name>` (template scaffold) but nothing that learns from a session.

#### Decision
Add session trace recording + skill distillation:
- `sciagent record start` / `sciagent record stop` — write a JSONL trace to `.sciagent/traces/<timestamp>.jsonl` capturing user prompts, agent responses, tool calls.
- `sciagent new skill <name> --from-trace <path>` — invoke a `trace-distiller` sub-agent (section 3.4) that produces a draft SKILL.md + `scripts/` skeleton.
- User reviews and edits before commit.

#### Consequences
- **Positive**: skill velocity 5–10x; users can author by demonstration; the distilled SKILL.md is dependency-aware (writes correct `requires:` based on observed tool calls).
- **Negative**: trace format coupling — Claude Code, Pi, and other harnesses log differently. Implement a normalizer.
- **Neutral**: requires `trace-distiller` agent (3.4) and one new command (4.3).

#### Alternatives
- Manual authoring only — slow, rejected.
- Full LLM-only auto-skill creation (no recording, just "describe a workflow") — produces weaker skills; the recorded trace is the ground truth.

---

### ADR-006 — First-class benchmarking via `sciagent bench`

**Status**: Proposed
**Date**: 2026-05-22

#### Context
Antigravity Skills' core contribution isn't the skills themselves — it's the *evidence the skills pass tests*. Internal capability benchmark (67 tasks, 49% → 93% reliability) + external benchmark (BioReason VEP, 41% → 61%). Without this, skills are unaudited.

`SciAgent-toolkit` has 50+ skills and no benchmarking. The first skill that silently breaks will not be discovered until a user complains.

#### Decision
Add three-tier benchmark harness (full spec in section 7):
- **Tier 1 — internal capability tests** (40–60 tasks, hand-curated, LLM-judged binary pass/fail).
- **Tier 2 — external benchmarks** (BixBench, OpenProblems.bio, BioReason as available).
- **Tier 3 — end-to-end workflow benchmarks** (2–3 full projects reproduced from raw inputs).

Invocation: `sciagent bench <tier> --role <name>`. Output: markdown report + JSON for diffing across versions.

#### Consequences
- **Positive**: regression detection; publishable evidence; concrete reliability claim for any role.
- **Negative**: significant upfront authoring (~2 weeks of work to seed Tier 1); ongoing maintenance.
- **Neutral**: opens path to a methods paper modelled on Antigravity Skills' own.

#### Alternatives
- Skip — accept silent drift.
- External CI only (GitHub Actions) — possible but disconnects from the role-activation workflow.

---

### ADR-007 — `sciagent doctor` validator

**Status**: Proposed
**Date**: 2026-05-22

#### Context
ADRs 001, 002, 003 introduce contracts the system must enforce. Currently no validator exists; broken skills are detected at runtime when an agent fails.

#### Decision
Add `sciagent doctor` subcommand with checks:
- Every SKILL.md parses as valid YAML frontmatter + markdown body.
- `metadata.sciagent.*` keys are namespaced correctly (ADR-001).
- `requires:` graph: all references resolve, no cycles, orchestrators have non-empty requires, atomics have empty requires (ADR-002, 003).
- Frontmatter passes agentskills.io reference validator (`uvx --from git+https://github.com/agentskills/agentskills validate`).
- Symlink hygiene: dangling symlinks in `.claude/skills/`, `.agents/skills/`.
- Manifest consistency: `.sciagent/manifest.json` matches actual filesystem state.

Output modes: `--json` for CI, plaintext for humans, `--fix` for safe auto-repair (re-create missing symlinks, remove dangling ones).

#### Consequences
- **Positive**: drift caught at author-time; CI-friendly; safe migrations.
- **Negative**: another command surface to maintain.
- **Neutral**: ~300 LOC of bash + the YAML parser already in `lib/sciagent/roles.sh`.

---

### ADR-008 — PathwayExplorer as overlay role, not separate repo

**Status**: Proposed
**Date**: 2026-05-22

#### Context
A postdoc or small group is unlikely to have *one* research project for the toolkit's lifetime. Projects pivot. Diseases change. Assay platforms get replaced. Hardcoding any single project's assumptions into the toolkit architecture is premature commitment.

But every comp-bio-driven scientific project does follow a recognizable shape:
- A **structured hypothesis space** (pathways, gene sets, candidate compounds, perturbations, …).
- An **experimental protocol selector** (which assay best probes the hypothesis).
- An **assay-output analyzer** (the dry side of the wet-lab loop — read raw data, extract signals).
- A **hypothesis ranker** (decide what to test next given accumulated evidence).
- An **orchestrator** that closes the loop.

This shape is genus-level invariant across virtually all empirical biology projects. The species-level differences (which pathway DB, which organoid model, which assay readout) are what change per project.

#### Decision
Two-part decision:

1. **First-class support for user-defined project overlays.** A `roles/<project-name>.yaml` overlay file authored by the user is treated identically to shipped overlays. The toolkit supports unlimited project overlays; the user picks one or more to activate per session.

2. **One shipping template overlay called `lab-loop`** that demonstrates the closed-loop pattern with deliberately under-specified, swappable skills. `lab-loop` is *not* a complete project — it's a starting point. Users fork it via `sciagent new role --template lab-loop --name <project>` and customize.

The template ships with five skill stubs (section 2.6) that illustrate the five-component shape above. Users replace each stub's body with project-specific implementation. The orchestrator skill remains as-is; the four atomic stubs are swap targets.

#### Consequences
- **Positive**: zero commitment to a specific project at toolkit level; the same toolkit serves dermatology, neuroscience, immuno-oncology, plant biology with different overlays; users can maintain multiple parallel projects (one overlay each) without forking the toolkit; community contributions of project-specific overlays become possible.
- **Negative**: the shipping template is intentionally incomplete — a new user might expect it to "just work" out of the box, which it doesn't.
- **Neutral**: matches how Anthropic ships `examples/` skills (illustrative, not production), and how the FutureHouse stack composes Crow/Falcon/Finch into Robin (the orchestrator carries the project shape; the agents are pluggable).

#### Alternatives
- Hardcode a specific project (e.g., the previously-proposed `organoquerry` overlay) — rejected, premature commitment.
- No template, only docs — possible but loses the working-example value that makes new overlays cheap to author.
- Ship multiple templates (e.g., `lab-loop`, `dry-only`, `multi-omic-loop`) — rejected for v1; revisit if `lab-loop` adoption shows clear gaps.


---

### ADR-009 — uv-pinned execution environment for sandboxed code

**Status**: Proposed
**Date**: 2026-05-22

#### Context
The Antigravity Skills paper reports significant reproducibility loss when sandboxed code runs in user-variable Python environments. They adopted `uv` to pin environments per skill execution. Same issue will hit any tree-search-based optimizer running LLM-generated code (ADR-004).

#### Decision
- `sandbox-execution` skill (2.3) executes code in a `uv`-managed venv defined per-skill or per-role.
- Default `pyproject.toml` template ships with `roles/scorable.yaml`.
- Each tree-search node runs in an isolated venv; the venv is cached and reused for sibling nodes with the same dependency set.

#### Consequences
- **Positive**: reproducible results; cache-friendly for tree search; matches DeepMind's stated best practice.
- **Negative**: `uv` is an additional system dependency; first-time setup adds friction.
- **Neutral**: `uv` is already the de-facto Python tooling in 2026; near-universal in the comp-bio Python world.

#### Alternatives
- conda envs — slower, larger, no improvement over `uv`.
- No pinning (rely on user env) — same problem Antigravity warned about.

### Personal take 
I tend to work with my analysis repos within Docker containers in an interactive manner and the Docker container has everything pre-installed most of the packages at least that I work with they are pre-installed into the container so I would not really let UV manage things entirely. Uh not that I don't see value in it, it's just that m for most part it seems like I have this pre configured for my personal use case. So I'm not sure. Does it give additional reproducibility if working within a container with all these tools pre-installed? I also don't know, so I would need additional judgment and thinking and grilling on that. 

---

### ADR-010 — Defer cloud-lab MCP integration

**Status**: Proposed
**Date**: 2026-05-22

#### Context
Emerald Cloud Lab is the only major commercial cloud lab still operating (Strateos folded 2023). ECL pricing (~$250K/year general access) puts it outside academic budgets. No production MCP server for ECL exists publicly. Wet-lab-in-the-loop is the field's next frontier but isn't accessible at small-lab scale yet.

#### Decision
- Document the *interface* the toolkit would expose if/when an ECL MCP server becomes available — a `cloud-lab-execution` skill stub with documented but unimplemented contracts.
- Do not implement actual cloud-lab integration in 2026.
- For wet-lab loop closure in OrganoQuerry, integrate with locally-operated assays via file-drop conventions (.FCS files, gene count matrices), not remote execution.

#### Consequences
- **Positive**: no wasted effort on integration that costs $250K/year to test; toolkit stays usable; preserves option value for 2027 when MCP wrappers likely emerge.
- **Negative**: OrganoQuerry's "transfer-acceleration" claim is measured against human-paced lab work, not robot-paced.
- **Neutral**: matches Robin's actual implementation pattern (file-drop, not robot integration).

#### Alternatives
- Build an OpenTrons MCP wrapper now — possible (hobbyist-grade ones exist), but OT-2 is too limited for flow cytometry + organoid culture. Defer.

---

## 2. New skills

All skill specs use ADR-001 namespacing. Each skill ships with `SKILL.md` plus optional `scripts/`, `references/`, `assets/`.

### 2.1 `tree-search-controller` (atomic)

```yaml
---
name: tree-search-controller
description: "ERA-style tree search over candidate code solutions for a scorable task. Use when you have (1) a baseline implementation, (2) a deterministic scoring function returning a float, and (3) a budget of 100-1000 candidate evaluations. Implements UCB1 selection, tracks tree state, and orchestrates mutator/scorer sub-agents. For non-scorable tasks use a regular agentic workflow; for hyperparameter sweeps under fixed code use a grid/random search."
license: MIT
metadata:
  skill-author: SciAgent-toolkit
  version: 0.1.0
  last-reviewed: 2026-05-22
  sciagent:
    scope: atomic
    requires: []
    complementary-skills: [quality-score-template, sandbox-execution, embedding-similarity]
    contraindications:
      - "Do not use without a deterministic scoring function — use quality-score-template first."
      - "Do not use for budgets <50 nodes — overhead dominates."
---
```

Body covers: tree representation, UCB1 formula (`score + c·sqrt(ln(N_total)/N_node)`), expansion policy, back-propagation of visit counts, termination criteria. Includes `scripts/run_search.py` reference implementation derived from `github.com/google-research/era`.

### 2.2 `quality-score-template` (atomic)

```yaml
---
name: quality-score-template
description: "Template and validation rules for writing a scoring function consumable by tree-search-controller. A valid score function is deterministic, returns a single float (higher = better), completes in <5min, and handles exceptions by returning -inf. Use before invoking tree-search-controller on a new task. For benchmarks already in the benchmark-loader catalogue, the score function is auto-supplied."
license: MIT
metadata:
  skill-author: SciAgent-toolkit
  version: 0.1.0
  last-reviewed: 2026-05-22
  sciagent:
    scope: atomic
    requires: []
    complementary-skills: [tree-search-controller, benchmark-loader]
    contraindications:
      - "Do not write non-deterministic score functions — tree search will not converge."
---
```

Body: rule-by-rule validation, three worked examples in `references/` (regression MSE, classification F1, custom comp-bio metric — e.g. scIB overall score), template `scripts/score_template.py`.

### 2.3 `sandbox-execution` (atomic)

```yaml
---
name: sandbox-execution
description: "Safely execute LLM-generated Python in a uv-pinned, resource-limited environment. Use as the execution layer beneath tree-search-controller, but also standalone when running any untrusted code. Enforces CPU/memory/wallclock limits, captures stdout/stderr/return value, isolates filesystem writes to a per-execution temp directory."
license: MIT
metadata:
  skill-author: SciAgent-toolkit
  version: 0.1.0
  last-reviewed: 2026-05-22
  sciagent:
    scope: atomic
    requires: []
    complementary-skills: [tree-search-controller]
    contraindications:
      - "Do not use for production deployments — this is a local sandbox, not a security boundary against adversarial code."
---
```

Body: `uv venv` per execution, `resource.setrlimit` for CPU/memory caps, `subprocess` with timeout, temp dir isolation, return value capture. `scripts/sandbox_run.py` reference.

### 2.4 `benchmark-loader` (atomic)

```yaml
---
name: benchmark-loader
description: "Uniform interface to public scientific benchmarks: BixBench, OpenProblems.bio (batch-integration, label-projection, etc.), BioReason (VEP-Coding, VEP-Non-SNV), GIFT-Eval. Each benchmark exposes load_data(), score(prediction), and metadata (split, metric, baseline). Use when running sciagent bench tier-2, or when invoking /optimize with --benchmark=<name>."
license: MIT
metadata:
  skill-author: SciAgent-toolkit
  version: 0.1.0
  last-reviewed: 2026-05-22
  sciagent:
    scope: atomic
    requires: []
    complementary-skills: [tree-search-controller, quality-score-template]
    contraindications: []
---
```

Body: per-benchmark loader specs in `references/<benchmark>.md`. Each `references/<benchmark>.md` includes dataset URL, schema, scoring metric formula, citation.

### 2.5 `embedding-similarity` (atomic)

```yaml
---
name: embedding-similarity
description: "Embed text or code snippets and compute pairwise similarity. Used by tree-search-controller to detect semantically duplicate branches and prune redundant exploration. Backend-agnostic — supports OpenAI, Gemini, Voyage, sentence-transformers embeddings via a unified API."
license: MIT
metadata:
  skill-author: SciAgent-toolkit
  version: 0.1.0
  last-reviewed: 2026-05-22
  sciagent:
    scope: atomic
    requires: []
    complementary-skills: [tree-search-controller]
    contraindications: []
---
```

### 2.6 OrganoQuerry domain skills

Five new skills under the `organoquerry` overlay. Stubs only — fill in based on disease/organ commitment (ADR-008, section 5.2).

- **`pathway-topology-query`** (atomic) — wraps KEGG REST, Reactome ContentService, OmniPath, SIGNOR. Returns subgraphs queryable by gene set, pathway ID, or disease ID.
- **`organoid-protocol-selector`** (atomic) — given a tissue type and disease, returns ranked protocol candidates from `references/protocols/` (curated by you; not auto-generated).
- **`flow-fcs-analyzer`** (atomic) — reads `.FCS` files, runs configurable gating (singlets, live/dead, marker-positive), computes MFI / fold-change vs control. Reuse `finch` code patterns from `github.com/Future-House/finch` (MIT-licensed).
- **`hypothesis-ranker`** (atomic) — given a pathway subgraph + readout, ranks candidate perturbation nodes by expected information value (Bayesian OED).
- **`organoquerry-orchestrator`** (orchestrator) — chains the four atomics into an end-to-end workflow. `requires: [pathway-topology-query, organoid-protocol-selector, flow-fcs-analyzer, hypothesis-ranker]`.

---

## 3. New sub-agents

Claude Code / Pi sub-agent format. Each ships as `agents/<name>.md`.

### 3.1 `mutator`
Stateless. Input: `(current_code, mutation_type, research_idea_or_skill_pointer)`. Output: rewritten code. Tools: `Read` (to load `complementary-skills` from the active role for the mutation alphabet). No persistent context across calls.

### 3.2 `scorer`
Wraps `sandbox-execution` + the active `score()` function. Input: `(code, dataset_path)`. Output: `float`. Handles timeouts, OOM, exceptions — returns `-inf` per the `quality-score-template` contract.

### 3.3 `selector`
Pure-Python implementation invoked as a sub-agent for tool-use parity. Input: tree state JSON. Output: `(node_id_to_expand, mutation_type)`. Implements UCB1 with semantic-similarity pruning via `embedding-similarity`.

### 3.4 `trace-distiller`
Reads a JSONL session trace. Identifies invariant operations (constant commands, API endpoints) vs variant (user inputs, dataset paths). Drafts a `SKILL.md` with invariants in the instructions, variants as parameters. Drafts a `scripts/` skeleton from observed tool calls. Outputs to a fresh `skills/<name>/` directory; user reviews and edits before commit.

---

## 4. New slash commands

### 4.1 `/optimize`
```
/optimize --metric=<path.py> --baseline=<path.py> [--budget=N] [--benchmark=<name>] [--mutation-mode=normal|ablation]
```
Launches tree search via `tree-search-controller`. Streams progress (current best score, expansion count) to chat. On termination, writes `output/best.py`, `output/tree.json`, `output/diff.md`.

`--mutation-mode=ablation` restricts the mutator's alphabet to "remove a required skill" or "swap for a complementary skill" — turns `/optimize` into an ablation runner.

### 4.2 `/explain-best`
Reads the most recent `output/tree.json`. Walks from baseline (root) to winner along the path of score-improving nodes. For each mutation, generates a human-readable explanation. Output: `output/explanation.md`. This is the methods-section auto-writer.

### 4.3 `/record start` / `/record stop`
Toggles session trace recording. Writes to `.sciagent/traces/<timestamp>.jsonl`. Recording is opt-in — never auto-on, never silent.

### 4.4 `/bench`
Triggers `sciagent bench` from inside the chat without dropping to terminal. Useful for "run tier-1 against current role" mid-session.

---

## 5. New roles

### 5.1 `scorable` (overlay)

```yaml
# roles/scorable.yaml
name: scorable
description: "ERA-style optimizer overlay. Activate on top of a base role when you have a scorable problem."
metadata:
  sciagent:
    role-type: overlay
    requires-base: true   # cannot activate solo
skills:
  - tree-search-controller
  - quality-score-template
  - sandbox-execution
  - benchmark-loader
  - embedding-similarity
agents:
  - mutator
  - scorer
  - selector
commands:
  - optimize
  - explain-best
output_style: empirical   # terse, score-focused
```

### 5.2 `organoquerry` (overlay)

```yaml
# roles/organoquerry.yaml
name: organoquerry
description: "Pathway-topology-driven hypothesis generation with flow-cytometry validation loop. Activate on top of compbio base."
metadata:
  sciagent:
    role-type: overlay
    requires-base: true
    disease: TBD          # pin once you commit (dermatology suggestion: atopic dermatitis or melanoma)
skills:
  - pathway-topology-query
  - organoid-protocol-selector
  - flow-fcs-analyzer
  - hypothesis-ranker
  - organoquerry-orchestrator
agents:
  - hypothesis-proposer    # reuses or specializes generation pattern from Co-Scientist / Robin
commands:
  - propose                # propose next experiment given current data
  - close-loop             # run hypothesis-rank → assay-suggest → analysis once new data arrives
```

Stacking: `sciagent activate compbio organoquerry` activates OrganoQuerry. Adding optimizer on top via `sciagent inject` once that's supported, or by redefining a deeper-stacked composite role.

**Note**: the current architecture (`docs/architecture.md` §2) caps stack depth at 2. Two overlays simultaneously (`scorable` + `organoquerry`) would require either (a) raising the cap to 3 — discouraged in the spec, (b) creating a composite role `organoquerry-scorable` that pre-merges both, or (c) using `sciagent inject` to add scorable's skills atop the organoquerry overlay. **Recommendation: option (c)** — keeps depth at 2 and matches the spec's intent for ad-hoc additions.

---

## 6. New CLI subcommands

### 6.1 `sciagent doctor` (ADR-007)
Validates skills, role YAML, symlinks, manifest. Modes: default, `--json`, `--fix`. Exit code 0 = healthy, 1 = warnings, 2 = errors.

### 6.2 `sciagent bench <tier>` (ADR-006)
- `tier-1` — runs internal capability tests (section 7.1).
- `tier-2` — runs external benchmarks (section 7.2).
- `tier-3` — runs end-to-end workflow benchmarks (section 7.3).
- `--role <name>` — restricts to a specific role's effective skill set.
- `--baseline <commit-or-tag>` — diffs against a previous run.
- Output: `bench-reports/<timestamp>/report.md` + `summary.json`.

### 6.3 `sciagent record start|stop|list` (ADR-005)
- `start` — begin trace capture; non-blocking, writes to `.sciagent/traces/<timestamp>.jsonl`.
- `stop` — close current trace, emit path.
- `list` — show recent traces with timestamp, duration, prompt count.

### 6.4 `sciagent new skill <name> --from-trace <path>` (ADR-005)
Invokes `trace-distiller`. Creates `skills/<name>/SKILL.md` + scripts skeleton. Opens in `$EDITOR` for review before commit.

---

## 7. Benchmarking framework

ADR-006 specifies the three tiers. Concrete specs follow.

### 7.1 Tier 1 — internal capability tests

40–60 hand-curated tasks. Each lives in `bench/tier1/<task-id>/` with:
- `prompt.md` — the user request as sent to the agent.
- `expected.yaml` — accepted outputs / patterns / known correct answers.
- `judge.md` — LLM-as-judge prompt scoring binary pass/fail against `expected.yaml`.
- `tags:` — domain (`single-cell`, `flow`, `proteomics`, `literature`), difficulty (`basic`, `intermediate`, `expert`).

Examples relevant to your stack:
1. "Annotate cell types in this 10K-cell AnnData using markers from CellMarker." → judge checks ≥80% of top-5 cell types present.
2. "Run GSEA on this DE result against KEGG pathways; report top 5 enriched terms." → judge checks specific term IDs are present in output.
3. "Given these `.FCS` files from a phagocytosis assay, produce gated MFI per well." → judge checks numeric values within ±10% of ground truth.
4. "Find all known ClinVar pathogenic variants in gene MITF associated with melanoma." → judge checks specific variant IDs are listed.
5. "Identify a CTCF candidate cCRE near TERT on chr5." → matches Antigravity Skills paper's reference task, gives direct cross-comparison.

Authoring rate: ~1–2 tasks per hour once template stabilizes. Total seeding: ~50 hours of work.

### 7.2 Tier 2 — external benchmarks

For each benchmark, a `benchmark-loader` reference + a `sciagent bench tier-2 --benchmark=<name>` invocation:
- **BixBench** — bioinformatics agent benchmark (referenced in Robin paper). Highest comparability with the field.
- **OpenProblems.bio batch-integration** — single-cell methods. Matches ERA's evaluation.
- **BioReason (VEP-Coding, VEP-Non-SNV)** — directly comparable with Antigravity Skills paper numbers.
- **PaperQA2 LitQA2** — for literature-search skills.

### 7.3 Tier 3 — end-to-end workflow

Pick 2–3 complete projects you've previously done by hand. For each:
- `bench/tier3/<project>/inputs/` — raw inputs only.
- `bench/tier3/<project>/expected_outputs/` — your published figures, tables.
- `bench/tier3/<project>/judge.md` — LLM-as-judge prompt scoring qualitative match.

Run: `sciagent bench tier-3 --project <name> --role <name>`. Measures wall-time, token cost, and qualitative reproduction quality. This is the *role-level* benchmark — not "does skill X work" but "does the role get from raw input to publishable output."

---

## 8. Roadmap

### 8.1 Four-week MVP

| Week | Goal | Deliverable |
|---|---|---|
| 1 | ADR-001/002/003/007 land | Skills migrated to namespaced metadata; `sciagent doctor` ships; no behavior change yet |
| 2 | `scorable` overlay v0 | `tree-search-controller` + `sandbox-execution` + `quality-score-template` skills shipped; `mutator`/`scorer`/`selector` agents shipped; `/optimize` works on a toy task (e.g., optimize a logistic regression on UCI dataset) |
| 3 | One real benchmark | Reproduce one ERA result on OpenProblems.bio batch-integration. Cost it on your hardware. This is the proof-of-life. |
| 4 | `sciagent bench tier-1` | Seed with 15–20 tier-1 tasks. Run against `base` and `base scorable`. Generate first comparison report. |

End of week 4 = blog post / arXiv preprint scope.

### 8.2 90-day strategic plan

**Days 1–30 (foundation)**:
- Land the four-week MVP above.
- Pick OrganoQuerry's target disease (ADR-008 disease pin).
- Identify one Chicago-area wet-lab collaborator for OrganoQuerry's eventual flow assay.

**Days 31–60 (loop closure)**:
- Ship ADR-005 (`sciagent record` + `--from-trace`).
- Author the first 5 OrganoQuerry skills (section 2.6).
- Tier-1 grows to 40–50 tasks.

**Days 61–90 (first results)**:
- Run OrganoQuerry on a small in-vitro pilot with the wet-lab collaborator.
- Bench tier-3 against one of your previous projects (you said you have a GSEA workflow — that's the obvious candidate).
- Draft a methods paper. Antigravity Skills paper template (intro → skill bundle → quality testing → results → limitations) is the cleanest model.

**Risk register**:
- *Cannot find a wet-lab collaborator*: OrganoQuerry becomes a comp-bio-only methods paper; smaller but still publishable.
- *Tree search doesn't beat baseline on first benchmark*: this happens; ERA's per-dataset solutions sometimes lose to foundation models. Pick a different benchmark with weaker baselines.
- *Spec collision on metadata keys*: ADR-001 mitigates; revisit if agentskills.io spec adds conflicting keys.

---

## 9. Open questions to decide before implementation

These are decisions that depend on you and aren't resolved by the design above.

1. **Migration strategy for ADR-001**. Big-bang codemod (one commit, all skills) or per-skill as touched? Recommend big-bang — keeps drift low.
2. **OrganoQuerry's disease commitment** (ADR-008). Atopic dermatitis (high pharma WTP, well-established cell models, your derm background) vs melanoma (organoid models more mature, harder competition from Recursion). My read: atopic dermatitis. Yours?
3. **Stack depth cap**. The architecture spec caps at 2. Stacking `scorable` + `organoquerry` requires the inject workaround (5.2 note). Acceptable, or worth raising the cap to 3 with explicit warnings?
4. **`uv` as hard dependency** (ADR-009). Anyone running SciAgent-toolkit needs `uv` installed. Acceptable for your target user base?
5. **External benchmark coverage in Tier 2**. BixBench is the obvious anchor — anyone scoring well on it is taken seriously. Is that enough for a first paper, or do you want OpenProblems.bio + BioReason for breadth?
6. **License for new skills**. Existing skills are MIT. Continue MIT, or switch to Apache 2.0 (matches agentskills.io spec license, gives patent grant)?
7. **Paper venue**. Bioinformatics? Nature Methods? Nature Communications? arXiv-only? Probably depends on Tier 3 results.

---

## 10. Things this spec deliberately doesn't address

- **MCP server authoring**. You have skills wrapping APIs; some might be cleaner as MCP servers. Out of scope for this spec; a separate ADR if you decide to take it on.
- **Multi-user / shared toolkit deployments**. Current design is per-user, per-project. If a Chicago lab wants to share roles across postdocs, that's a future ADR.
- **Cost tracking**. `/optimize` will burn tokens and cloud compute. A `sciagent budget` subcommand tracking spend would be useful but not foundational.
- **Cloud-lab integration** (ADR-010 explicitly defers).
- **Career strategy** (out of scope for a technical spec; see the conversation, not this file).

---

## 11. If you implement nothing else, do these three

In rough priority order:

1. **ADR-001 + ADR-007**: namespace + validator. ~1 day. Zero behavior change. Future-proofs everything. *Do this regardless.*
2. **ADR-006 Tier-1 only**: seed 20 hand-curated comp-bio tasks; build the LLM-judge harness. ~1 week. Gives you the evidence base every subsequent decision needs. *Do this before adding more skills.*
3. **ADR-004 `scorable` overlay**: the ERA-style optimizer. ~2 weeks. *This is the publishable contribution.* Either ship it or don't claim ERA-style search in any paper.

The remaining ADRs (002, 003, 005, 008, 009) are all valuable but compose well later. ADR-010 is a non-decision.

---

*End of spec. Sleep well.*
