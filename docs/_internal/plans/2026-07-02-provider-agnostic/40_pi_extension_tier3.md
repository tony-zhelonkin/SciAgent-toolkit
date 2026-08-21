# Tier 3 — Hard: the pi role-routing extension

**Effort:** ~2–4 weeks. **Risk:** medium–high (new subsystem; pins to pi's extension API).
**Depends on:** Tier 2 pi adapter (context+skills+settings). **Delivers:** SciAgent *roles* routed into
pi as **parallel or chained isolated sub-agents**, auto-selected from the role menu or named
explicitly — "orchestrate the way I would," not just "inject context."

This is the tier that dissolves the old blocker. Prior-art (`ai-research/02`) said *"pi has no
Task/sub-agent primitive."* **That is now wrong**: the installed pi package
(`@earendil-works/pi-coding-agent@0.79.4`) ships a complete working sub-agent example
(`examples/extensions/subagent/`), and `@aliou/pi-guardrails@0.13.1` is a canonical guardrail-extension
package. Both are on this machine and are the reference implementations.

---

## T3.1 — Ship SciAgent as a pi *package* (not a loose extension)

Model on `@aliou/pi-guardrails` (its `package.json` declares `"pi":{"extensions":[...]}`, ships
`extensions/`, `src/`, a JSON-schema-backed config, and `commands/` submenus). The SciAgent pi package:

```
sciagent-pi/                      # vendored into the toolkit; published as a pi package
  package.json                    # "pi": { "extensions": ["./dist/index.js"], "skills": [...] }
  src/index.ts                    # export default (pi: ExtensionAPI) => { ... }
  src/roles.ts                    # role discovery + role→pi-agent-md generation
  src/fanout.ts                   # vendored/adapted from examples/extensions/subagent
  src/guardrails.ts               # pi.on("tool_call") CRAFT/style enforcement
  schema.json                     # extension config (default trust, concurrency, role dir)
```

**Vendor, don't depend** (per `ai-research/11` "don't depend on any single pi extension"). Pin the
tested pi version (0.79.4). The sub-agent fan-out is an extension *convention*, not a core guarantee —
own the code.

## T3.2 — Map a SciAgent role → a pi agent

A SciAgent **role** (`roles/*.yaml` = `skills[] + agents[] + commands[] + output_style`) generates a pi
agent markdown (pi's `agents/*.md` convention: frontmatter `name/description/tools/model`, body =
context):

- **frontmatter** `name` = role name; `description` = the role's one-liner (this is what auto-selection
  matches on); `model`/`tools` from the role/overlay.
- **body** = the role's `system-prompt` (from `system-prompts/`, already provider-agnostic) + the CRAFT
  guardrail block + a manifest of the role's skills (pi loads `.agents/skills/` natively, so the skills
  ride along without copying).
- Generated into `.pi/agents/<role>.md` (project) and/or `~/.pi/agent/agents/<role>.md` (user), with
  pi's `agentScope` gating. **Note:** pi's sub-agent example discovers `~/.pi/agent/agents/` and
  `.pi/agents/` — *not* SciAgent's `.agents/agents/`. The extension must point discovery at whichever
  tree we standardize on (trivial one-liner in the vendored `discoverAgents`).

## T3.3 — Fan-out modes (both already demonstrated in the pi example)

The vendored fan-out registers one tool with three modes (from `examples/extensions/subagent/`):

- **Single** `{role, task}` — one isolated `pi -p --mode json` child playing the role.
- **Parallel** `{tasks:[{role,task}...]}` — up to 8 tasks, 4 concurrent, isolated contexts, streamed,
  usage-tracked, Ctrl-C-propagating. *This is the "true parallel sub-agents just the way I want them"
  the user asked for.*
- **Chain** `{chain:[{role,task}...]}` — sequential agent-over-agent with a `{previous}` placeholder
  threading each output into the next. *This is the "sequential agent over agent" behavior.*

Each child inherits the SciAgent package (so guardrails + skills apply to sub-agents too, not just the
parent). Outputs cap at ~50 KB to the parent; full result in tool details.

## T3.4 — Role selection: explicit AND automatic (both native to pi)

- **Explicit** — a `pi.registerCommand("role", ...)` slash command (`/role:science-architect <task>`)
  or naming the role in the prompt. Deterministic, user-driven.
- **Automatic / unsupervised** — register the fan-out as a **tool** whose `description` enumerates the
  available roles and their one-liners. pi's progressive-disclosure system-prompt injection lets the
  model pick a role by matching task→description — the same mechanism skills use. A `session_start`
  hook injects the active-role menu so the model knows the options up front.

This satisfies the user's stated design: *route a role either manually via prompt or automatically in
an unsupervised way where a fanned-out agent pulls up the necessary role context on the fly and plays
that role.*

## T3.5 — Guardrails on parent and children

Use the `@aliou/pi-guardrails` pattern: `pi.on("tool_call", handler)` can `return {block:true, reason}`
to enforce SciAgent coding-style / safety rules (no-ephemeral, figure-save discipline, provenance)
*before* a tool runs — on the parent and, because children load the same package, on the sub-agents.
This is layer (c) "unskippable" realized natively in pi, complementing the `validate` floor.

## T3.6 — Trace substrate (pulls ADR-011 forward)

pi writes a **documented `version: 3` JSONL session schema** (`~/.pi/agent/sessions/`) plus a separate
content-redacted observability event stream. Claude writes an undocumented JSONL. The deferred
cross-harness plan (`ai-research/03`) proposed *"sciagent-native JSONL on disk + optional OpenInference
exporter."* If Tier 3 ships sub-agent orchestration, we now have a **second trace consumer** — which is
exactly the documented trigger to un-defer the recorder (`ai-research/11:96`: *"ship trace recording iff
ADR-004 or ADR-006 also ships"* — generalize: iff a second consumer appears). **Promote ADR-011 —
sciagent trace/result substrate** (the "missing ADR" flagged in `ai-research/10:220`) and decide it
before Tier 3 fan-out lands, so sub-agent runs are captured in one normalized format.

---

## T3.7 — Feasibility verdict

**Highly feasible.** Every primitive the design needs is public and demonstrated in the *installed* pi
package: `export default (pi)=>`, `pi.registerTool`, `pi.registerCommand`, `pi.on(event)`, child
`pi -p --mode json` spawning, `agents/*.md` discovery with scope gating, and a blocking `tool_call`
guardrail hook. The work is integration + role→agent generation + vendoring the fan-out, not research.

### Caveats to design around
- Fan-out is an extension convention → pin pi 0.79.4, vendor the example, add a version check.
- Project-local agents are trust-gated → headless/CI needs `defaultProjectTrust:"always"` or `-a`.
- Sub-agent output is capped (50 KB/task) → design role tasks to return structured summaries.
- Agents re-discovered per invocation (hot-editable) → good for authoring, but cache role→md generation.

## Tier-3 acceptance checks

- [ ] `pi` in a SciAgent project loads the package; `/role:science-architect "…"` runs an isolated
      role sub-agent.
- [ ] Parallel mode fans 3 roles concurrently; chain mode threads output through 3 roles sequentially.
- [ ] Auto-selection: with no explicit role, pi picks the right role from task→description matching.
- [ ] `tool_call` guardrail blocks an ephemeral-file write in both parent and a sub-agent.
- [ ] Sub-agent runs are captured in the ADR-011 trace format (if ADR-011 shipped).
- [ ] Package is vendored + version-pinned; no runtime dependency on a third-party pi extension.
