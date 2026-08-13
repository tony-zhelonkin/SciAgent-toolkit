> **Updated** 2026-05-24 with grounded findings from local repos at `docs/.ref/pi/` and `docs/.ref/science-skills/`. Targeted edits to ADR-005, ADR-006, and cross-cutting observations.
> **Updated** 2026-05-24 (second pass) — deployment-shape and trace-necessity grilling deferred to file 11; pointers added in ADR-005, ADR-004, the decision queue, and cross-cutting observations.
> **Updated** 2026-05-24 (third pass) — ADR-004's sub-agent split bisected in file 11 §4a; pointer added below. For the architectural / engineering / taste taxonomy applied across all ADRs and a reading/decision order, see `docs/kickoff.md`.

# 10 — ADR-on-ADRs: a grilling matrix

> **Outcome note (superseded).** Superseded by the three-verb toolkit
> architecture documented in [`docs/architecture.md`](../../architecture.md).
> Retained as a dated research record; the original body below is unchanged.

For each of the 10 ADRs in `sciagent-extension-design-spec.md`, four sub-headers — decisional crispness, findings that change the picture, open questions for grilling, sequencing & risk. Then: recommended decision queue, cross-cutting observations.

---

## ADR-001 — Namespace custom metadata under `metadata.sciagent.*`

**Decisional crispness.** Crisp. Yes/no on whether to migrate. The *what* (namespace prefix) is well-specified.

**Findings that change the picture.** File 05 verifies the "reasonably unique key names" quote — ADR-001 is correctly motivated. **But:** the agentskills.io spec types `metadata` as `string → string`. Current sciagent skills (and ADR-001's proposal) put *structured* values (lists, nested maps) under `metadata`. ADR-001 inherits a mild spec violation. The skills-ref reference validator may flag it. Test before migrating. See file 05.

**Open questions for grilling.**
1. Has skills-ref been run against one current sciagent skill to confirm whether structured `metadata.sciagent.*` values pass or fail validation?
2. Is the migration big-bang (one commit, ~50 skills via codemod) or incremental (touch-as-you-go)? Spec §9.1 leans big-bang; what's the rollback plan if a skill silently breaks?
3. Why namespace `requires/scope/contraindications/complementary-skills/tier` under `sciagent` while keeping `version/last-reviewed/category/tier/tags/upstream-docs/skill-author` at `metadata.*`? The line between "sciagent-specific" and "generic ecosystem field" is fuzzy — see `tier` appearing in both lists in the current skills. Whose convention says `tier` is generic?
4. If agentskills.io spec evolves to standardise `requires` with different semantics, what's the migration path *out* of `metadata.sciagent.requires`?

**Sequencing & risk.** Foundational for ADR-002 (graph), ADR-003 (scope), ADR-007 (validator). Cheap to land first, low irreversibility — namespace renames are codemod-friendly. Should land before any new skills are written.

---

## ADR-002 — Formalize `requires` as a hard dependency graph

**Decisional crispness.** Design-space. Multiple defensible answers: soft vs hard, version-pinned vs name-only, in-file vs lockfile-external.

**Findings that change the picture.** File 09 shows cycle detection is cheap (Python `graphlib` or 30 lines of bash). The existing repo already has `lib/sciagent/skill_deps.sh`. Adding cycle detection there is incremental, not foundational.

**Open questions for grilling.**
1. Version pinning syntax: `[{name: harmonypy, version: ">=0.1.0"}]` is heavyweight. Most current skills have empty `requires:`. Is a less-typed form (`[harmonypy@>=0.1.0]`) preferable for the 95% case?
2. What does "resolved to latest" mean for an unversioned `requires: [harmonypy]`? Latest *in this checkout* (deterministic, but git-state-dependent) or latest *globally* (non-deterministic)?
3. Lockfile yes or no? (Not in ADR-002, but implied. File 09 sketches a `.sciagent/skills.lock` pattern. Worth grilling — is this premature?)
4. Soft vs hard for `complementary-skills`: spec keeps it soft. Confirm this is intentional and not asymmetry-for-its-own-sake.

**Sequencing & risk.** Depends on ADR-001 landing first (graph needs to live under a namespace). Foundational for ADR-007 (validator) and ADR-004 (ablation searches walk the graph). Reasonably cheap to land — the existing dep-resolution code is most of the way there. Hard to undo if you ship version constraints and later want to remove them — pinning becomes load-bearing.

---

## ADR-003 — Two-tier scope: `orchestrator` vs `atomic`

**Decisional crispness.** Crisp on the binary, but the spec also acknowledges "mid-tier gymnastics" cases. Worth confirming the rule for the boundary case is acceptable.

**Findings that change the picture.** Current skill `coresh-signature-search/SKILL.md` already has `scope: atomic` at `metadata.*` level. Migration is mostly renaming/relocating, not redesigning. Aligns with how anthropics/skills implicitly structures skills (file 05).

**Open questions for grilling.**
1. Why two tiers, not one (just `requires:` carrying the load)? Spec says "loses authoring guidance" — what guidance specifically, and does an enforced lint rule do the same work without a separate field?
2. The current `tier: standard` field across skills (see `coresh-signature-search`) overlaps with `scope`. ADR-003 says "drop `tier`" — confirm there are no consumers of `tier` elsewhere in the toolkit (`bin/sciagent`, agents, role YAML).
3. Are there genuinely mid-tier skills in the current 50+? Spot-check 2–3 and ask: orchestrator or atomic? If you find any where the answer is "really both," the binary is leaky.

**Sequencing & risk.** Lands with ADR-002 — they share migration. Low risk; ~50 skills × 1-line edit. Reversible.

---

## ADR-004 — Optimizer as overlay role (`scorable`), not baseline

**Decisional crispness.** Crisp on placement (overlay vs baseline). The contents of the overlay (skills 2.1–2.5, agents 3.1–3.3, commands 4.1–4.2) are a separate decision tree.

**Findings that change the picture.** File 06 verifies ERA exists, is open-source (Apache-2.0), and architecturally matches the spec's description. Naming nit: actual algorithm is "Flat UCB Tree Search (FUTS)," not "UCB1" — skill 2.1's body should reflect this. The ~20K-token-per-session estimate for the overlay (ADR-004 context) is unsourced; check before quoting in any external write-up.

**Open questions for grilling.**
1. The overlay assumes problems are scorable with a deterministic float. Single-cell integration via scIB *is* scorable; "interpret this DE result" is not. What fraction of sciagent use cases are genuinely scorable, and is that fraction big enough to justify a separate overlay?
2. Stack-depth-2 cap + `scorable` overlay forces an inject workaround when combining with `organoquerry` (spec §5.2 note). Is this acceptable, or does ADR-004 quietly motivate raising the cap (which spec explicitly resists)?
3. ERA repo is Jupyter-notebook-heavy (90%). The skill 2.1 plan ships `scripts/run_search.py` "derived from" ERA. How much code reuse vs reimplementation — is the Apache-2.0 attribution path clean?
4. Per-session cost: if a single `/optimize` run costs $50–500 in API calls (BixBench-scale), is interactive sessions the right interface, or should this run async / out-of-process? `sciagent bench --baseline=...` is async; `/optimize` in chat is not.

**Sequencing & risk.** Foundational for ADR-006 (the optimizer is what bench measures). Depends on ADR-002 (graph) for ablation-mode mutations. Depends on ADR-009 (uv) for sandboxed execution. Significant work (~2 weeks in spec roadmap). Hard to undo if you ship `/optimize` in chat and later need to make it async — that's a UX-breaking change.

**See file 11 §4** for: the deployment-shape decision applied to ADR-004 — the mutator/scorer/selector sub-agents assume Claude Code's Task tool, which Pi has no analog for (confirmed via four community Pi extensions in `docs/.ref/pi-extensions-examples/` that all simulate sub-agents via child Pi processes, not in-process primitives). Recommended shape: **out-of-process `sciagent optimize <objective>` CLI** rather than in-chat tool — matches Q4 above (per-session $50–500 is batch-shaped, not chat-shaped) and avoids the harness-coupling tax entirely. In-chat invocation can come later as a thin adapter if interactive UX matters.

**See file 11 §4a** for the *sub-agent bisection*: applying the agents-best-practices rule "subagents only when decomposition improves measured results" splits the three proposed sub-agents into one (`mutator`, [Arch] decomposition kept on ERA's borrowed evidence, mark for local re-A/B when ADR-006 tier-1 lands) plus two tools (`selector` — [Eng] drop sub-agent framing, "tool-use parity" doesn't pass the bar; `scorer` — [Taste] either sub-agent or tool, decided by retry/error preferences). Net effect: three sub-agents collapse to one + two tools. Less surface, less coordination tax, aligns naturally with the out-of-process shape.

---

## ADR-005 — `workflow-skill-creator` analog via session trace recording

**Decisional crispness.** Research question disguised as a decision. The spec proposes both (a) trace recording (`sciagent record`) and (b) distillation (`new skill --from-trace`). They are separable; the spec couples them.

**Findings that change the picture.** Substantial.
- File 04: there is **no Anthropic `workflow-skill-creator`**. The existing `skill-creator` interviews users, doesn't ingest traces. A `workflow-skill-creator` *does* exist in `docs/.ref/science-skills/skills/workflow_skill_creator/SKILL.md` (Google DeepMind), but it is **also interview-driven, not trace-driven** — uses agent memory of the just-completed conversation as seed, not a structured session log. ADR-005's trace-ingesting variant is genuinely net-new in both ecosystems. Reframe as net-new; cite the DeepMind skill as a precedent for the *phase model*, not for trace ingestion.
- File 01: Claude Code already writes JSONL transcripts; `transcript_path` is in every hook payload. Capture cost is near-zero — copying or symlinking the file.
- File 02 (revised): Pi coding agent writes JSONL with a **documented version-3 schema** (`docs/.ref/pi/packages/coding-agent/docs/session-format.md`), has a **fully documented hook surface** with ~30 events (`docs/.ref/pi/packages/coding-agent/docs/extensions.md:268-335`), and a **separate observability layer** for content-redacted spans (`docs/.ref/pi/packages/agent/docs/observability.md`). Cross-harness work is real but bounded — a documented-to-documented field mapping, ~150 LOC each direction.
- File 02 also: Pi has **no Task-tool / subagent primitive**. A Pi-side distiller must ship as `pi.registerCommand("sciagent-distill", …)` slash command, not an in-conversation sub-agent.
- File 03: do not adopt OTel GenAI or LangSmith as the primary format. Define a sciagent-native minimal JSONL schema with optional OpenInference exporter. Pi's content-vs-spans split suggests sciagent should also split the schema rather than collapse both into one file.

**Open questions for grilling.**
1. Does ADR-005 actually need to ship trace *recording*, or just trace *consumption*? Claude Code already records; Pi already records. sciagent could ship only `new skill --from-trace <path>` and consume existing files. Recording-as-hook is then a convenience, not a foundation.
2. The trace-distiller (sub-agent 3.4) is genuinely new. What's the success criterion — ship-ready SKILL.md, or draft for human edit? The Google DeepMind `workflow_skill_creator` enforces user approval gates between four phases (brainstorm → design → implement → validate, `docs/.ref/science-skills/skills/workflow_skill_creator/SKILL.md:23-270`). Should the sciagent distiller adopt the same gate model, or attempt one-shot autonomous distillation?
3. Cross-harness deployment is asymmetric: on Claude Code the distiller is a Task-tool sub-agent; on Pi it must be a slash command (no Task analog — confirmed at `docs/.ref/pi/packages/coding-agent/src/core/tools/`). Is "Claude-only sub-agent + Pi-only slash command + shared distiller logic in a package" the right factoring, or should sciagent ship one out-of-process distiller invoked the same way from both harnesses?
4. Should the sciagent trace schema split content (JSONL: prompts, tool calls, results) from spans (event log: timings, span tree, redacted metadata), mirroring Pi's session-vs-observability split? File 03's current single-schema proposal conflates them.
5. Privacy: traces contain everything the agent saw (file paths, prompt content, API outputs). Spec says "opt-in — never auto-on, never silent" — confirm this includes "never auto-pushed to a remote distiller." Pi's observability layer is content-redacted by default (`docs/.ref/pi/packages/agent/docs/observability.md:266-297`) — useful precedent for safe-by-default payloads.

**Sequencing & risk.** Decouples cleanly from ADRs 001–003. Depends on the trace-format decision (file 03). Medium risk — if the sciagent-native JSONL schema is wrong, every distiller and downstream tool inherits the wrongness. Spec underweights this risk. Pi's harness side is now lower-friction than initially scoped: recording is `pi.on("turn_end", …) + pi.appendEntry(…)`, ~50 LOC TypeScript; the unknowns are all on the distiller and cross-harness-deployment side.

**See file 11** for: (a) the deployment-shape decision matrix (out-of-process vs. core+adapters vs. adopt an existing Pi extension) and the recommended *core + thin per-harness adapters* shape, (b) the trace-necessity verdict: defer recording until ADR-004 or ADR-006 commits this cycle — ship ADR-005 as interview-driven with optional `--from-trace <path>` that consumes harness-native JSONL directly, no sciagent-native recorder yet. ADR-005 should be split into (distiller verb) + (recorder hook), with the second deferred.

---

## ADR-006 — First-class benchmarking via `sciagent bench`

**Decisional crispness.** Compound. The three tiers are three separate decisions. Tier-1 is crisp; tier-2 has open choices on which benchmarks; tier-3 is project-specific and the design is more guideline than spec.

**Findings that change the picture.** File 08 changes tier-2:
- BixBench: tier-2 anchor. Real cost (1–2 days wire-up, $200–800/run, 5–48h wallclock).
- OpenProblems batch-integration: second anchor, directly comparable to ERA. ~3–5 days wire-up.
- BioReason VEP: optional, lower priority.
- **GIFT-Eval: drop.** No bio relevance. Spec includes it because ERA used it; sciagent is bio-focused.
- PaperQA2 LitQA2: optional, depends on whether sciagent ships literature-search skills.

File 04 (revised): the "three-tier eval pattern" attributed to "Antigravity Skills paper" is unverifiable. Confirmed absent from `docs/.ref/science-skills/` itself — the Google DeepMind repo ships skills only, no eval harness, no `evals/` directory, no benchmark numbers. A linked tech report PDF (`docs/.ref/science-skills/README.md:62-64`) may contain the missing eval methodology and numbers, but is not in the local ref. Own the three-tier pattern as sciagent-original; if a citation surfaces in the DeepMind PDF, attribute then.

**Open questions for grilling.**
1. Tier-1 task authoring: 40–60 hand-curated tasks at ~1–2 tasks/hour = 50 hours of work. Who writes them, and against what? Each task needs `prompt.md`, `expected.yaml`, `judge.md`. Is the LLM-as-judge prompt rigorous enough to avoid grade inflation?
2. Tier-2 budget: BixBench is $200–800 per run × multiple runs (with baseline, across versions) = real money. Who pays?
3. Tier-3 is "previously-done projects reproduced from raw inputs." How many do you have ready? Spec says "2–3" — is one publication's worth enough?
4. `--baseline <commit-or-tag>` is a regression-detection feature. Are you set up to diff JSON across versions cleanly?
5. Cadence: `sciagent bench` per commit (CI gate) or on-demand? Tier-2 cost makes per-commit infeasible.

**Sequencing & risk.** Depends on ADR-005 for trace capture if benchmarks include "judge by inspecting trace." Depends on ADR-004 for `--role scorable` benchmarks. Foundational for any "ship more skills" claim — without benchmarks, skills are unaudited. High risk to defer; medium-high cost to land. Tier-1 only is the cheap entry.

---

## ADR-007 — `sciagent doctor` validator

**Decisional crispness.** Crisp on the binary (ship a validator) and the check list. Implementation choices have multiple defensible answers (file 09).

**Findings that change the picture.** File 09:
- skills-ref does **not** check cross-skill dependency graphs. ADR-007 must implement ADR-002/003 invariants natively.
- JSON Schema via `check-jsonschema` is recommended over Yamale — IDE integration is free, CI integration is one-liner.
- `--fix` mode should be narrowed to symlinks + manifest only; never auto-edit SKILL.md frontmatter.
- ~300 LOC bash + 50–150 LOC for schema check matches spec's own estimate.

**Open questions for grilling.**
1. JSON Schema vs Yamale vs Pydantic — pick a winner. Defaults aside, what does the rest of the toolkit code style favour?
2. `--fix` scope: narrow to safe ops (file 09 recommendation)? Or wider?
3. CI integration target: GitHub Actions only, or generic? `--json` mode covers the latter.
4. Does `doctor` validate the AGENTS.md managed-block hash? Architecture doc §4 makes this a first-class check.
5. Skills-ref delegation: per-skill linter or aggregated? Cost matters if `doctor` runs in CI.

**Sequencing & risk.** Depends on ADRs 001/002/003 to define what to validate. Can land as a stub before those — "validate YAML parses, validate symlinks" is useful immediately. Foundational for ADRs 005, 006 (both add files `doctor` should check). Low risk, decoupled from the optimizer/benchmark threads.

---

## ADR-008 — PathwayExplorer as overlay role, not separate repo

**Decisional crispness.** Two compound decisions packaged as one. (a) Support user-defined project overlays as first-class. (b) Ship a `lab-loop` template overlay. Both defensible separately.

**Findings that change the picture.** None from web research — this is internal architecture. Worth noting: the architecture doc's stack-depth-2 cap is restated in §5.2 with an inject workaround. ADR-008 + ADR-004 together create the case for raising the cap to 3, which the spec explicitly resists.

**Open questions for grilling.**
1. Is a *template* the right unit, or should it be a *generator*? Templates atrophy if not maintained; generators stay current with toolkit changes.
2. The five skill stubs (pathway-topology-query, organoid-protocol-selector, flow-fcs-analyzer, hypothesis-ranker, organoquerry-orchestrator) are extremely specific to one project shape (comp-bio + wet-lab loop closure). Are they the right primitive abstractions, or too narrowly scoped?
3. Spec §9.2 picks a disease (atopic dermatitis vs melanoma) as a follow-up decision. ADR-008 explicitly says "TBD." When does TBD become a blocker?
4. Multi-overlay stacking (scorable + organoquerry) requires the inject workaround. Live with it, or design overlays to compose natively?
5. If `lab-loop` is "deliberately incomplete," what's the first-run experience for a user activating it? Empty stubs that error helpfully, or pre-filled examples?

**Sequencing & risk.** Decoupled from the optimizer/benchmark threads but coupled to ADR-002 (`requires:` carries the orchestrator-to-atomic links). Premature commitment risk if you ship the template with too-strong opinions; toothless-spec risk if too weak. This is the only ADR where rolling-back commitments is genuinely expensive — once you ship a "this is how to structure a comp-bio project" template, it sets expectations.

---

## ADR-009 — uv-pinned execution environment

**Decisional crispness.** Crisp on the binary, but the user's own "Personal take" section flags the decision is unresolved.

**Findings that change the picture.** File 06 shows ERA's published material does not explicitly justify uv adoption (the abstract/blog/README don't mention it). ADR-009's appeal to "Antigravity Skills paper" reproducibility loss is unsourced. File 07 argues for a hybrid: optional in base role, mandatory in `scorable` overlay.

**Open questions for grilling.**
1. The user's container has 95% of packages pre-installed. For *interactive* analysis, does uv add anything they need? (File 07 says: a lockfile + per-project isolation, mostly nice-to-haves.)
2. For *LLM-generated code in tree search* (the ADR-004 use case), is there any path without uv? (File 07 says: no, the pre-installed container doesn't help when search generates novel imports.)
3. Hybrid (mandatory in `scorable`, optional in base) vs uniform (mandatory always). Which is less work to maintain and explain?
4. If skill 2.3 `sandbox-execution` is the only consumer of uv, can the uv dependency be local to that skill's `scripts/`, not a system-wide sciagent dep?

**Sequencing & risk.** Depends on whether ADR-004 lands. If ADR-004 ships first, ADR-009 is forced. If both deferred, no urgency. Low irreversibility — switching from "uv required" to "uv optional" is easier than the reverse.

---

## ADR-010 — Defer cloud-lab MCP integration

**Decisional crispness.** Crisp. Defer. The cost analysis ($250K/year for ECL) makes this a non-decision.

**Findings that change the picture.** None.

**Open questions for grilling.** Minimal.
1. The "cloud-lab-execution skill stub with documented but unimplemented contracts" — is even this stub worth the maintenance burden if nothing consumes it?
2. ADR-010 is two paragraphs of "we are not doing this." Why is it an ADR rather than a one-line note in §10 (Things this spec deliberately doesn't address)?

**Sequencing & risk.** No risk. No blockers. Could be a paragraph in §10 instead.

---

## Recommended decision queue

> File 11 reshuffles this. Grill-first now includes the deployment-shape question (out-of-process vs core+adapters vs adopt-existing-extension); ADR-005 splits into a distiller verb (do now, interview-driven) and a recorder hook (defer unless ADR-004 or ADR-006 ships same cycle).

### 1. Grill first — foundational + unresolved
- **Deployment shape** (file 11 §2). Decide before ADR-005 implementation: harness-agnostic core + thin per-harness adapters is recommended. Sets the substrate for every other in-conversation verb sciagent adds.
- **ADR-004** (`scorable` overlay). Foundational for ADR-006; the publishable contribution. Multiple sub-decisions need answers (per-session cost UX, interaction with stack-depth cap, code reuse from ERA repo). File 11 §4 recommends out-of-process `sciagent optimize` CLI shape.
- **ADR-006** (benchmark harness). Compound decision; tier-2 needs explicit pruning (drop GIFT-Eval, defer BioReason/PaperQA2); tier-1 authoring is 50+ hours of work that needs an owner.
- **ADR-005a** (distiller verb only — interview-driven, optionally `--from-trace <harness-native-path>`). Smaller decision than the original ADR-005. Independent of recorder.

### 2. Move to implementation — crisp, low-risk, decoupled
- **ADR-001** (namespace metadata) — but first test that structured values under `metadata.sciagent.*` pass skills-ref validation.
- **ADR-003** (orchestrator/atomic scope) — lands with ADR-001.
- **ADR-007** (`sciagent doctor`) — can land as a stub immediately, grow with ADR-001/002/003.

### 3. Defer or split
- **ADR-002** (hard dependency graph). Crisp in intent, but version-pinning syntax and lockfile question are deferrable. Land *cycle detection* + *resolution check* now; defer version semantics + lockfile to ADR-002.5.
- **ADR-005b** (sciagent-native JSONL recorder + cross-harness normalizer). Defer until ADR-004 or ADR-006 commits this cycle and provides a second consumer. Otherwise interview-driven distillation covers the common case; file 11 §3 decision rule.
- **ADR-008** (`lab-loop` template overlay). Coupled to project-disease decision which is also TBD. Split into (a) "first-class user-defined overlays" — already supported by current architecture, almost a non-ADR — and (b) "ship `lab-loop` template" — premature until ADR-002 lands.
- **ADR-009** (uv-pinned execution). Defer until ADR-004 ships and forces the decision. Soft-warn in `sciagent doctor` meanwhile.
- **ADR-010** (defer cloud-lab). Demote to a §10 paragraph; do not ship as an ADR.

---

## Cross-cutting observations

- **ADRs 004, 005, 006 all depend on a trace/result substrate that isn't itself an ADR.** Tree-search node histories (ADR-004), distillation traces (ADR-005), benchmark scoring logs (ADR-006) all want the same persistence layer. The sciagent-native JSONL schema in file 03 is the missing ADR. Consider promoting it: "ADR-011 — sciagent trace/result substrate."

- **The spec's empirical claims are systematically under-cited.** "49% → 93% internal reliability," "41% → 61% BioReason VEP," "ERA reports reproducibility loss with non-uv envs," "~20K tokens per session for the overlay" — none of these have URLs. Before any of this ships in a methods paper, audit every empirical claim for source.

- **The `.agents/` convention is more aspiration than implementation.** Confirmed at Pi source (`docs/.ref/pi/packages/coding-agent/src/core/package-manager.ts:2278-2339`): Pi natively walks `.agents/skills/` only — `agents/` and `commands/` siblings are never read. Claude Code doesn't use `.agents/` at all. The dual-track symlink topology (architecture §7) is sciagent-only — calling it "harness-agnostic" overstates what's portable. Pi *can* be pointed at Claude Code's `.claude/skills/` directly via settings.json (`docs/.ref/pi/packages/coding-agent/docs/skills.md:43-62`), making the dual-track partially redundant on Pi. Worth being honest about this in any external pitch.

- **Stack-depth-2 cap is rubber-banding under load.** ADR-004 (scorable overlay) + ADR-008 (project overlay like lab-loop) immediately want two overlays on top of base = depth 3. The spec workaround (`sciagent inject`) is a sleight-of-hand: it adds skills to an existing overlay rather than creating a new tier, but the user-visible behavior is "I now have three things stacked." Either accept depth 3 with explicit warnings or design overlay composition more carefully. This isn't an ADR yet but probably needs to be.

- **The "Antigravity Skills paper" referent is wrong, and now more precisely so.** Spec invokes it 3+ times. There is no Anthropic Antigravity paper. The local Google reference (`docs/.ref/science-skills/`) confirms: Google DeepMind ships a `workflow_skill_creator` skill (and 36 other science skills) under the Google Antigravity "Science" plugin, but the skill is **interview-driven, not trace-driven** — closest in shape to Anthropic's `skill-creator`, not to ADR-005's trace-distiller. The repo carries no eval harness or benchmark numbers; those would live (if anywhere) in the linked but not-locally-available DeepMind PDF. The attribution across the spec needs to: (a) drop "Anthropic Antigravity Skills paper," (b) optionally cite "Google DeepMind science-skills bundle" for the phase-model and CLI-script-default conventions if those are adopted, (c) source-hunt the eval numbers in the DeepMind PDF before quoting them anywhere external.

- **"Spec hygiene" thread (ADRs 001, 002, 003, 007) is cheap and decoupled; "loop closure" thread (005, 006) is expensive and entangled.** The cheap thread should ship first regardless of how the expensive thread resolves. Spec §11 already says this — re-read it.

- **Harness-agnostic-vs-coupled is itself a cross-cutting axis** (file 11 §1, §2). Every ADR that adds an in-conversation verb (`record`, `optimize`, `bench`, `new skill --from-trace`) implicitly takes a position. The four community Pi extensions (`pi-actors`, `pi-multiagent`, `pi-subagents`, `pi-packages`) converge on a tiny extension backbone — single tool registration, optional skills bundle, optional slash command, optional `session_shutdown` hook — small enough that the cost of writing a thin per-harness adapter against a harness-agnostic CLI core is bounded (~150 LOC TypeScript per harness, comparable to the ~50-LOC trace-recording figure from file 02). Adopting one of the four as a substrate is not viable (each carries 2K–20K LOC of opinions sciagent does not share). The user's standing harness-agnostic preference is satisfiable; the cost is the adapter layer, not a wholesale core rewrite.

- **Out-of-process is the right default for batch-shaped work.** ADR-004's optimizer ($50–500/run) and ADR-006's benchmarks (BixBench $200–800/run, 5–48h wallclock) are not chat-shaped. File 11 §4 makes this a recommendation for ADR-004; the same logic applies to ADR-006. In-conversation invocation can come later as a thin adapter — it should not be the v1 shape.
