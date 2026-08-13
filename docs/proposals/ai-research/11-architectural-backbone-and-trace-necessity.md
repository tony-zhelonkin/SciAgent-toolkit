> **Created** 2026-05-24 from local repos at `docs/.ref/pi-extensions-examples/` and `docs/.ref/agents-best-practices/`.
> **Updated** 2026-05-24 (second pass) — §4 expanded with the mutator/scorer/selector bisection derived from the agents-best-practices "subagents only when decomposition improves measured results" rule. See `docs/kickoff.md` for the architectural/engineering/taste taxonomy applied across ADRs.

# 11 — Architectural backbone and trace-recording necessity

> **Outcome note (superseded).** Superseded by the three-verb toolkit
> architecture documented in [`docs/architecture.md`](../../architecture.md).
> Retained as a dated research record; the original body below is unchanged.

Third grilling pass. Two halves: (a) what shape do real community Pi extensions converge on, and (b) does sciagent actually need a JSONL trace recorder right now? The user's two standing preferences — harness-agnostic where engineering cost allows, and minimal-but-no-less — are applied as constraints, not decoration.

---

## §1. Pi extension patterns — convergence and divergence

Four community extensions sit under `docs/.ref/pi-extensions-examples/`. Each one is solving "Claude-Code-style functionality on Pi" from a different angle.

**`pi-actors` (`@llblab/pi-actors` v0.20.1).** Single TypeScript extension (`./index.ts`, registered via `package.json:"pi".extensions`) exposing three durable verbs — `spawn`, `message`, `inspect` — plus a file-discovered registry at `~/.pi/agent/recipes/*.json` (`pi-actors/README.md:166-181`). The contract is operator-managed executable memory: filename = tool id, JSON template owns execution shape, run-actor owns lifecycle. Registers tools (`registerTool`), TUI widgets (`ctx.ui.setStatus`/`setWidget`, see `pi-actors/index.ts:64-100`), and slash commands (`/actors-inspector-toggle`). No out-of-process component for the recipe runner — actors are local processes spawned via `pi.run` templates. Skills bundled in-package (`pi-actors/skills/actors/SKILL.md`, `swarm/SKILL.md`) but pointed at by `package.json:"pi".skills`. Heaviest of the four; tries to be a full local-orchestration kernel.

**`pi-multiagent` (`pi-multiagent` v0.9.4).** One Pi tool: `agent_team` with seven sub-actions (`catalog`, `start`, `run_status`, `step_result`, `message`, `cancel`, `cleanup`) — `pi-multiagent/README.md:30-49`. Static-DAG graphs, child Pi processes launched as detached runs, evidence-first (artifacts not chat). Extension entry at `extensions/multiagent/index.ts:30-58`. Ships nine package agents (`agents/*.md`) and one skill (`skills/pi-multiagent/SKILL.md`). Hooks `session_shutdown` to cancel live runs (`extensions/multiagent/index.ts:38-44`), registers a flag (`registerFlag`), registers message renderer (`registerMessageRenderer`). Children inherit nothing by default — no parent transcript, no context files, no themes, no skills unless `--agent-team-subagent-skills enabled`. No out-of-process CLI; orchestration lives entirely inside the Pi session via the registered tool.

**`pi-subagents` (`pi-subagents` v0.25.0, originally `nicobailon/pi-subagents`).** A larger extension (`src/extension/index.ts`) registering one tool (`subagent`), eight built-in agents (`agents/*.md`), seven prompt templates (`prompts/*.md`), and slash commands (`/run`, `/chain`, `/parallel`, `/run-chain`, `/subagents-doctor`). Optional companion extension (`pi-intercom`) for child→parent comms. Has a thin out-of-process installer (`install.mjs`) that just `git clone`s into `~/.pi/agent/extensions/subagent/` — not a runtime CLI, only a bootstrap. Discovers agents from three scopes — builtin / user (`~/.pi/agent/agents/`) / project (`.pi/agents/` and legacy `.agents/`) — `pi-subagents/README.md:354-364`. Most Claude-Code-shaped of the four: explicit Task-tool-style delegation with foreground/background, chains, parallel, clarify UI.

**`pi-packages` (`@gotgenes/pi-packages`).** A pnpm monorepo of four siblings: `pi-autoformat`, `pi-github-tools`, `pi-permission-system`, `pi-subagents` (a fork of `tintinweb/pi-subagents`, not Bailon's). Each subpackage is an independent extension with the same shape — `package.json:"pi".extensions: ["./src/index.ts"]` — and bundled skills. Confirms the package-as-extension convention rather than introducing new primitives. Notable: `@gotgenes/pi-subagents` claims "Claude Code look & feel — same tool names, calling conventions, and UI patterns (`Agent`, `get_subagent_result`, `steer_subagent`)" (`pi-packages/packages/pi-subagents/README.md:18-19`) — an explicit deliberate Claude-port.

**Convergent backbone.** Across all four:

- Single TypeScript entry per extension, registered via `package.json:"pi".extensions: ["./path.ts"]`.
- One or a small number of `pi.registerTool(...)` calls. Sub-action shapes (`agent_team {action: ...}`, `subagent {action: ...}`) prefer one tool with discriminated-union actions over many narrow tools.
- Bundled skills via `package.json:"pi".skills`, each as `<dir>/SKILL.md`.
- Optional slash commands via `pi.registerCommand` for operator-facing workflows.
- Lifecycle hooks via `pi.on("session_start" | "session_shutdown" | "turn_end" | ...)` — almost always to cancel/cleanup detached state, not to mutate model behavior.
- Detached/long-running work persisted as JSON artifacts on disk, indexed by short run-ids; never as in-memory state alone.
- Child Pi processes for sub-agent behavior, never an in-process Task primitive — Pi has none (confirmed file 02, `pi-multiagent` and `pi-subagents` both shell out to `pi -p` children).

**Divergent design space.**

- *Trust model for children.* `pi-multiagent` starts children with nothing (no transcript, no skills, all-or-nothing skill propagation). `pi-subagents` defaults to forked context and selective skill injection. `pi-actors` treats actors as fully trusted local programs from the start.
- *Coordination surface.* `pi-actors` has rooms + mailboxes. `pi-multiagent` has static DAGs. `pi-subagents` has chains + parallel + saved `.chain.md` workflows. No convergence here yet.
- *State scope.* `pi-actors` writes recipes to `~/.pi/agent/recipes/` (user-global). `pi-subagents` reads agents from three scopes (builtin/user/project). `pi-multiagent` keeps run state in-process and tmp artifacts.

**Minimal architectural backbone for Pi-side Claude-Code-style functionality.** The smallest extension that delivers a Claude-Code-shaped sub-agent / Task-tool / slash-command surface on Pi is:

1. One `package.json` declaring `"pi".extensions: ["./index.ts"]` and optionally `"pi".skills: [...]`.
2. One `pi.registerTool({name, schema, run})` for any in-conversation delegate behavior.
3. Optional `pi.registerCommand(...)` for any out-of-flow operator action.
4. Optional `pi.on("session_shutdown", ...)` for cleanup of detached state.

That is the full convergent core. Everything else — actor rooms, static DAGs, clarify UIs, parallel groups, foregrounded TUI widgets — is choice, not necessity. None of the four extensions can be cleanly *depended on* by sciagent as a substrate: they have overlapping ambitions but incompatible coordination models, and each carries 2K–20K LOC of TUI/lifecycle/policy that sciagent does not need.

---

## §2. Three deployment shapes for sciagent

sciagent today is a bash CLI (`bin/sciagent`) that wires symlinks into `.claude/` and `.agents/`. The ADR spec proposes extending it with `record`, `optimize`, `bench`, `doctor`, and `new skill --from-trace` verbs. The shape question: how does the extension reach the harnesses' in-conversation surface?

**A. Harness-agnostic out-of-process core only.** sciagent stays a bash/Python CLI invoked identically from both harnesses through a slash command (`/sciagent record`, `/sciagent doctor`) or a registered tool wrapper. Both harnesses already support shelling out: Claude Code via Bash tool, Pi via `pi.registerTool` calling `bash`. *Engineering cost:* lowest — no per-harness code, no TypeScript dependency, current bash/Python continues to work. *Adoption cost:* friction inside conversations — operators must remember to invoke `/sciagent X` explicitly; nothing surfaces in the TUI lifecycle; no `tool_use` block carries a structured sciagent result. *Where it fails:* the trace-distiller use case (ADR-005) needs to run right after the conversation ends; an out-of-process CLI can be triggered by the user but cannot react to `turn_end` or `session_shutdown` events.

**B. Harness-agnostic core + thin per-harness adapter.** Same standalone core (bash + Python), plus a ~150-LOC Claude Code sub-agent shim (markdown + frontmatter, calls the core via Bash tool) and a ~150-LOC Pi extension shim (TypeScript, registers one `sciagent` tool plus optional `pi.on("session_shutdown", ...)` hook that calls the core). The shims are *thin* — they translate harness events into core CLI calls and translate core JSON output into harness-native UX. The core never imports harness types. *Engineering cost:* moderate — two adapters to maintain, each tracking their harness's API. Bounded because the surface is small (file 02 showed Pi recording is ~50 LOC; the inverse — invoking the core — is similar). *Adoption cost:* low — each harness gets native UX. *Where it fails:* requires sciagent to keep the core's CLI surface stable enough for adapters to wrap; if the core's CLI grows N verbs, both adapters need N updates.

**C. Adopt an existing Pi extension as base.** Pick one of the four reviewed extensions and depend on it.
- `pi-actors` is closest in spirit (file-discovered tools, recipe registry) but heavyweight (~20KB `index.ts`, full actor kernel) and its abstraction — verbs over addresses — is orthogonal to sciagent's role/skill/symlink model. **Disqualified by size.**
- `pi-multiagent` exposes a static-DAG delegation tool; sciagent has no DAG. **Disqualified by mismatch.**
- `pi-subagents` is the closest Claude-Code-port; it ships `/run`, `/chain`, `/parallel` and an `Agent`-shaped tool. Could plausibly host sciagent's distiller as one of its registered agents. **Closest fit, but** it carries strong opinions on agent scopes (builtin/user/project), chain files, intercom bridge, clarify UI — adopting it would import all of that. And sciagent's existing `agents/*.md` already follow a Claude-Code-style frontmatter scheme that overlaps but doesn't exactly match `pi-subagents`'s agent format. **Verdict: none qualifies as a clean dependency.** The convergent backbone (§1) is small enough that rolling it is cheaper than absorbing any single extension's opinions.

**Recommended shape: B.** Rationale, against the two standing preferences:

- *Harness-agnostic preference:* B keeps the load-bearing logic in a harness-neutral core. Both adapters are thin and dispensable — a third harness (Codex, an `agentskills.io`-compatible runtime, etc.) gets a third adapter, not a fork of the core. A is *more* harness-agnostic but loses the in-conversation seam that ADR-005 needs. C optimizes for Pi at the cost of Claude Code and any future harness.
- *Minimal-but-no-less preference:* the adapter is the minimum needed to satisfy the one concrete use case (post-conversation distillation hook) that A cannot satisfy. ~150 LOC TypeScript per harness is the floor identified by file 02 and confirmed by the convergent backbone in §1 (single-tool registrations, optional hook). It is the smallest seam that gives sciagent ride-along access to harness lifecycle without coupling the core.
- *Evidence from §1:* the four extensions converge on a tiny registration surface. Building a thin shim against that surface is a few-hundred-LOC effort, not a multi-week project. The cost is real but defensible.

A is the right answer if the user decides ADR-005 should ship as "invoke `sciagent record` manually after a session" rather than as a `turn_end` hook. That decision belongs to §3.

---

## §3. The JSONL trace-recording necessity question

The user's open question: *"Should we even design and architect a record command that would record and trace the JSON log of conversation in some particular manner?"*

### 3a. Where JSONL trace recording brings unique value

Claude Code already writes JSONL transcripts (file 01); Pi writes JSONL sessions with documented v3 schema (file 02). The question is not *will the trace exist* — it does — but *will sciagent consume it programmatically*?

The interview pattern (Google DeepMind's `workflow_skill_creator`, Anthropic's `skill-creator`) substitutes the agent's own recollection of the just-finished conversation for the structured trace. It works well when the conversation is recent and the agent can be trusted to articulate intent. It fails specifically where:

1. **Cross-session distillation.** "Turn the last *three* analyses I did into a skill." The agent does not remember three sessions ago. Only the on-disk traces do. Interview cannot recover what the agent never had in context.
2. **Tool-call-pattern extraction.** "Which Bash commands did I run that mattered, in what order, with what intermediate failures?" Interview produces a sanitized recollection; the trace has the actual command sequence including the dead ends. For sciagent's domain — turning real analysis sessions into reusable atomic skills — the dead ends are load-bearing (they become the "Gotchas" section that file 01 of the agents-best-practices skill insists on, `references/skills-and-connectors.md:88-117`).
3. **Ablation / regression replay (ADR-004, ADR-006).** Tree-search mutation needs node histories; benchmark scoring needs per-tool-call timing and result diffing. Interview cannot reconstruct either. The cross-cutting observation in file 10 ("ADRs 004, 005, 006 all depend on a trace/result substrate that isn't itself an ADR") is the real motivation — trace recording is a *substrate* whose first consumer happens to be distillation but whose other consumers are the optimizer and the benchmark harness.
4. **Pi observability split.** Pi separates content-bearing session JSONL from content-redacted observability spans (file 02). Interview pattern cannot use the spans; trace consumption can. Span-level timing data is exactly what a benchmark harness needs.

### 3b. Where the interview pattern is strictly better

- **Single-session "save this as a skill"** — the most common case. The agent saw the whole conversation; the four-phase interview gates (file 04: brainstorm → design → implement → validate) produce a cleaner SKILL.md than autonomous trace ingestion. Distillation of a one-shot trace is *more* work than the interview, not less.
- **Skill creation by non-technical users.** Interview asks "what did we just do?" — easy. Trace ingestion asks "where is your `~/.claude/projects/.../session_id.jsonl`?" — already friction.
- **When the trace is contaminated.** Long sessions where the agent went down a wrong path and was course-corrected. Articulated intent cleanly excludes the dead branch; trace distillation has to decide what to drop, and the decision is hard.
- **Tax of maintaining a recorder.** A `sciagent record` hook is small (~50 LOC TypeScript per harness per file 02 + §2 above). The schema underneath is *not* small — file 03 already noted this. Every schema change ripples to distiller, optimizer, benchmark. The schema is the load-bearing artifact, and the cost is ongoing.

### 3c. Decision rule

**Ship trace recording now if and only if ADR-004 (optimizer) or ADR-006 (benchmark) is also being implemented this cycle.** ADR-005 alone does not justify it — the interview pattern covers the common case at lower cost. The substrate cost is only amortized when more than one consumer reads from it.

Concretely:
- If the user is implementing ADRs 001/003/007 this cycle (the "cheap thread" from file 10), defer the recorder. Ship ADR-005 as an interview-driven `sciagent new skill` analog with optional `--from-trace <path>` for the cross-session and ablation cases — *consume* the harness-native JSONL files directly, do not *record* into a sciagent-native schema yet.
- If the user commits to ADR-004 (`scorable` overlay) or ADR-006 (`bench`) this cycle, ship the recorder alongside, because both will need the schema regardless. In that case the recorder is paying for itself across three consumers.

**Re-evaluation trigger.** Revisit if any of: (a) the interview-driven distiller starts producing skills that miss dead-end gotchas users actually hit, (b) more than two sciagent verbs need to consume harness output, or (c) a third harness (e.g., Codex) joins and the cost of writing per-harness consumers exceeds the cost of one normalized schema.

---

## §4. What this means for ADRs 005 and 004

**ADR-005 should be reframed as two separable sub-decisions, not one.** Per §3:

1. *Distillation entrypoint* — the `sciagent new skill --from-X` verb. This is the user-facing decision and is largely independent of the recorder question. The interview gate model from Google DeepMind (file 04) is the right shape; adopt it. The verb takes `--from-conversation` (interview), `--from-trace <path>` (consume harness-native JSONL), or `--from-traces <dir>` (cross-session). All three feed the same distiller.
2. *Recorder* — `sciagent record` as a per-harness hook that writes a sciagent-native schema. Defer this until either ADR-004 or ADR-006 lands. Until then, the `--from-trace` path consumes Claude Code or Pi JSONL directly (the normalizer is bounded work per file 02; both schemas are documented).

The deployment-shape decision from §2 (recommend B) applies to both halves: the interview-driven distiller is core CLI logic; the per-harness hook is a thin adapter that calls the core only after the conversation ends.

ADR-005 grilling in file 10 already raised the deployment-asymmetry question (Q3: "Claude-only sub-agent + Pi-only slash command + shared distiller logic?"). §2 answers it: shared core + per-harness adapter, with the adapter being a slash command on Pi (no Task analog) and a sub-agent on Claude Code. Both adapters call the same `sciagent distill` CLI verb.

**ADR-004 (`scorable` overlay) inherits the same deployment matrix.** The spec's `mutator`/`scorer`/`selector` sub-agents assume Claude Code's Task tool. Per §1, Pi has no Task primitive — child Pi processes are how `pi-multiagent` and `pi-subagents` simulate one. Three options, briefly:

- A (out-of-process only) — the tree search runs as `sciagent optimize <objective>` in a terminal, not in a conversation. *Best fit for §3's "defer recording" path*: the optimizer writes its own trace files on its own schedule, no harness recorder needed. The per-session cost concern raised in file 10 (Q4: "$50–500 per `/optimize` run — should this run async / out-of-process?") points the same direction.
- B (core + adapters) — the tree search exposes a thin tool on each harness that delegates to the out-of-process optimizer. Useful if interactive `/optimize` UX matters; not necessary if not.
- C (adopt `pi-subagents` on Pi) — would let the mutator/scorer/selector ride the `subagent` tool. Same disqualification as §2: imports too many opinions for too little win.

**Recommended for ADR-004: A**, with B as an upgrade path. The optimizer is the use case where the out-of-process shape is *most* defensible — sessions are long, costs are real, and the work is genuinely batch-shaped. Putting `/optimize` in chat is the wrong UX even before the harness-coupling question. File 10's Q4 already framed this; §2 confirms it.

### 4a. Inside the optimizer: the sub-agent bisection

The agents-best-practices rule "*subagents only when decomposition improves measured results*" (`docs/.ref/agents-best-practices/SKILL.md:172`, `references/checklists.md:164`) is sharp enough to bisect the three proposed sub-agents in ADR-004 individually. Each fails or passes the bar for a different reason:

- **`selector` — fails the bar; use a tool, not a sub-agent.** It is deterministic Python running Flat UCB Tree Search math. The spec's justification is "*invoked as a sub-agent for tool-use parity*" (`sciagent-extension-design-spec.md` §3.3). "Tool-use parity" is architectural prettiness, not measured benefit. Wrapping a pure function in a sub-agent boundary adds prompt overhead, latency, and a place for context to leak — for zero accuracy gain. Drop the sub-agent framing; expose `selector` as a plain tool the orchestrator calls.

- **`scorer` — passes the bar by structural argument, not performance.** The scorer runs LLM-generated code in a sandbox and returns a float. Sandbox isolation is a *structural* requirement, not a "decomposition improves results" claim. Sub-agent or tool both work; this is a [Taste] call about retry semantics and error surfaces, not an [Arch] one.

- **`mutator` — keep as sub-agent on borrowed evidence; mark for local revalidation.** The mutator's prompt structure (creative code rewriting given a research idea) is genuinely different from the orchestrator's (decide what to expand next). ERA's paper measures the decomposed shape against monolithic baselines and the decomposed shape wins. Borrowed evidence is honest provided it is *cited as borrowed* — re-A/B locally when ADR-006 tier-1 lands. This is the only place in ADR-004 where the rule's measurement bar genuinely bites.

Net effect: three sub-agents collapse to one sub-agent (`mutator`) plus two tools (`selector`, `scorer`). Strictly less surface, less prompt overhead, less coordination tax. Aligns naturally with the §4 out-of-process `sciagent optimize` shape — out-of-process tools are easier than out-of-process sub-agents anyway.

### 4b. The rule's bootstrap dependency on ADR-006

Applying "*decomposition only when measured*" requires a measurement apparatus — which is ADR-006. The rule creates a soft coupling: ADR-004's decomposition decisions inherit prior evidence (ERA's published numbers) until ADR-006 tier-1 lands and produces local A/Bs. This is fine — the rule is a *retrospective accountability check*, not a *prospective shipping gate*. Ship the decomposed `mutator` on borrowed evidence; mark a re-eval trigger ("revisit when tier-1 produces ≥5 task A/Bs"). The cross-cutting observation in file 10 ("ADRs 004, 005, 006 share a substrate that isn't an ADR") was already this; the bisection makes it concrete.

---

## Summary

The convergent Pi extension backbone is small (one tool registration, optional skills bundle, optional slash command, optional `session_shutdown` hook). sciagent should ship as a harness-agnostic core (bash + Python) plus thin per-harness adapters (§2-B), not as a per-harness rewrite. JSONL trace recording should not ship as part of ADR-005 alone; defer it until ADR-004 or ADR-006 commits this cycle. ADR-005 splits cleanly into (interview-driven distiller now) + (recorder when the second consumer arrives). ADR-004 belongs out-of-process as `sciagent optimize`, not in chat.
