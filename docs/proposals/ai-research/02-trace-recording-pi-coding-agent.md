> **Updated** 2026-05-24 with grounded findings from local repo at `docs/.ref/pi/`.
> **See also** file 11 §3 for the *necessity* decision: this file describes how to record (~50 LOC TypeScript hook) and confirms it's tractable, but file 11 argues recording should be deferred unless ADR-004 or ADR-006 ships in the same cycle, so a second consumer amortizes the schema-maintenance cost. The Pi extension shape itself (single `pi.registerTool` + `pi.on("session_shutdown", ...)`) is the convergent backbone identified across `docs/.ref/pi-extensions-examples/` — file 11 §1 documents this; no change to the recording mechanics described here.

# 02 — Trace recording in Pi coding agent

## What "Pi coding agent" actually is

Identifiable. Pi is an open-source terminal coding harness developed by earendil-works (Mario Zechner), published as `@earendil-works/pi-coding-agent`. The `pi-mono` monorepo contains four packages: `pi-ai` (provider-agnostic LLM API), `pi-agent-core` (runtime), `pi-coding-agent` (CLI), `pi-tui` (TUI library). License MIT, harness in TypeScript. The earlier `@mariozechner/pi-coding-agent` npm name and `badlogic/pi-mono` GitHub path are pre-rename references.

- repo root: `docs/.ref/pi/README.md`
- packages: `docs/.ref/pi/packages/{agent,ai,coding-agent,tui}/`

## Hook / extension surface — fully documented, richer than first thought

Pi exposes an in-process, typed event-driven extension API documented end-to-end in `docs/.ref/pi/packages/coding-agent/docs/extensions.md`. The prior "only `tool_call` shown publicly" hedging was wrong — the full event taxonomy is enumerated, with payload shapes, mutation semantics, and ordering guarantees. The hooks-runtime design lives at `docs/.ref/pi/packages/agent/docs/hooks.md`.

Complete event surface (`docs/.ref/pi/packages/coding-agent/docs/extensions.md:268-335`):

| category | events |
|---|---|
| session | `session_start`, `session_before_switch`, `session_before_fork`, `session_before_compact`, `session_compact`, `session_before_tree`, `session_tree`, `session_shutdown` |
| resources | `resources_discover` |
| agent | `before_agent_start`, `agent_start`, `agent_end`, `turn_start`, `turn_end` |
| messages | `message_start`, `message_update`, `message_end` |
| provider | `before_provider_request`, `after_provider_response`, `context` |
| tool | `tool_execution_start`, `tool_execution_update`, `tool_execution_end`, `tool_call` (blocking), `tool_result` (mutating) |
| input | `input` (intercept/transform/handle), `user_bash` |
| model | `model_select`, `thinking_level_select` |

Extensions are TypeScript modules auto-discovered from `~/.pi/agent/extensions/`, `.pi/extensions/`, npm packages, or `--extension <path>`. They receive a typed `ExtensionAPI` (`pi.on(...)`, `pi.registerTool(...)`, `pi.registerCommand(...)`, `pi.appendEntry(...)`, `pi.sendMessage(...)`) and an `ExtensionContext` with `ctx.ui`, `ctx.sessionManager`, `ctx.signal`. Mutation semantics are formally specified per event in `docs/.ref/pi/packages/agent/docs/hooks.md:156-271` — observers vs handlers, sequential transform chains, early-exit on `cancel`/`block`, etc.

The dedicated `hookMessage`/`custom` entry type (`docs/.ref/pi/packages/coding-agent/docs/session-format.md:25-26`) lets extensions persist arbitrary state into the session JSONL itself via `pi.appendEntry()` — survives restart, never enters LLM context.

In short: Pi's hook story is more capable than Claude Code's out-of-process hook scripts. Recording traces from Pi is a `pi.on("turn_end", …)` + `pi.appendEntry()` exercise.

## Session persistence — documented schema

Sessions are JSONL at `~/.pi/agent/sessions/--<path>--/<timestamp>_<uuid>.jsonl` (`docs/.ref/pi/packages/coding-agent/docs/session-format.md:7-9`). Currently `version: 3`; auto-migrated on load. Tree-structured via `id`/`parentId` per entry, enables `/tree` navigation, `/fork`, `/clone`, `/resume`, `/export [file]` (HTML), `/share` (private gist).

Entry types are formally defined (`docs/.ref/pi/packages/coding-agent/docs/session-format.md:39-169`): `SessionHeader`, `UserMessage`, `AssistantMessage`, `ToolResultMessage`, `BashExecutionMessage`, `CustomMessage`, `BranchSummaryMessage`, `CompactionSummaryMessage`. Each carries `id`, `parentId`, `timestamp`. The TypeScript type union (`AgentMessage`) lives in `docs/.ref/pi/packages/agent/src/types.ts`.

This is a stable, documented schema — adequate to consume without source-diving. The prior file's "format is undocumented" hedge is rescinded. The schema differs in shape from Claude Code's JSONL (different field names, different role taxonomy, `bashExecution` and `branchSummary` are Pi-specific), so a normalizer is still required for ADR-005, but it's a known-to-known translation.

## Observability layer — separate from hooks

Pi has a dedicated observability subsystem distinct from extension hooks (`docs/.ref/pi/packages/agent/docs/observability.md`). It defines runtime-agnostic `traceId`/`spanId` records with `runWithPiContext` + `traceOperation` primitives, and a stable event name list: `pi.agent.prompt`, `pi.agent.turn`, `pi.agent.tool_call`, `pi.ai.provider.request`, `pi.agent.session.append_entry`, etc. Designed to bridge to OTel/Sentry via adapters without binding Pi to either. Default payloads explicitly redact prompts, completions, tool args, and tool results (safe-by-default — `docs/.ref/pi/packages/agent/docs/observability.md:266-297`).

For ADR-005 this matters: Pi already separates "what the user said and the tool returned" (session JSONL, content-bearing) from "spans and timings" (observability events, content-redacted). A sciagent trace-recorder consuming Pi has both surfaces available.

## Skills loader — `.agents/skills/` only, not `.agents/agents/` or `.agents/commands/`

Confirmed at source. `docs/.ref/pi/packages/coding-agent/src/core/package-manager.ts:2278-2339`:

```text
// Project skills from .agents/ (each with its own baseDir)
for (const agentsSkillsDir of projectAgentsSkillDirs) { … addResources("skills", …) }
// User skills from ~/.agents/ (with its own baseDir)
addResources("skills", collectAutoSkillEntries(userAgentsSkillsDir, "agents"), …)
```

Only the `skills/` child of `.agents/` is walked. `agents/` and `commands/` siblings are not consumed. Discovery rules (`docs/.ref/pi/packages/coding-agent/docs/skills.md:36-39`): "In `~/.agents/skills/` and project `.agents/skills/`, root `.md` files are ignored" — directories with `SKILL.md` only. There is no `pi.agents` or `pi.commands` analog in Pi's resource model.

Pi does, however, allow Claude-side skills directories via `settings.json` (`docs/.ref/pi/packages/coding-agent/docs/skills.md:43-62`):

```json
{ "skills": ["~/.claude/skills", "../.claude/skills"] }
```

So sciagent's dual-track symlink topology (skills under `.claude/skills/` *and* `.agents/skills/`) is partially redundant on the Pi side — Pi can read directly from `.claude/skills/` if pointed there. The dual-track is genuinely useful only as a harness-portability convention for *other* harnesses that read `.agents/`. Worth re-stating that this is sciagent's own convention, not Pi's.

Name validation: Pi deliberately diverges from the agentskills.io spec — it does **not** require skill `name` to match parent directory, citing shared skill directories used across harnesses as the rationale (`docs/.ref/pi/packages/coding-agent/docs/skills.md:138-157`). This is a downstream-friendly choice that helps sciagent (and matches Antigravity science-skills, where directory `alphagenome_single_variant_analysis` carries name `alphagenome-single-variant-analysis`).

## Subagent / Task-tool mechanism — none

No `Task` tool. `docs/.ref/pi/packages/coding-agent/src/core/tools/` ships seven built-ins: `bash`, `edit`, `edit-diff`, `find`, `grep`, `ls`, `read`, `write`. No subagent-spawning primitive. Sub-task patterns are user-built via `pi.sendUserMessage(…, { deliverAs: "followUp" })`, `ctx.newSession({ withSession })`, `ctx.fork(entryId, { withSession })`. These are session-replacement primitives, not isolated short-lived sub-agents. Extensions can implement a "summarizer subagent" or similar — see `examples/extensions/summarize.ts` reference — but it's a manual pattern, not a tool the LLM calls.

For ADR-005 this implies: a Pi-side trace recorder cannot use a Task-tool-driven distillation handoff. Distillation has to live as a `pi.registerCommand("distill", …)` slash command or a post-session pipeline. Claude Code's pattern (Task tool spawning the distiller) is not portable to Pi.

## Implications for ADR-005

1. **Pi trace recording is now a tractable extension, not research.** The needed surface is `pi.on("session_start", …)` + `pi.on("turn_end", …)` + `pi.appendEntry("sciagent-trace", …)` or a sibling JSONL file. ~50 lines of TypeScript. The earlier "research gap" framing is overstated.

2. **Cross-harness JSONL portability remains zero by default.** Pi `version 3` JSONL and Claude Code JSONL share the parent-pointer DAG idea but disagree on field names (`role` taxonomy, `customType`, `toolCall.id` vs `tool_use.id`), entry types (Pi has `bashExecution`, `branchSummary`, `compactionSummary` natively; Claude Code does not), and timestamp format (Pi: Unix ms in `timestamp`; Claude Code: ISO 8601). A normalizer is real work — but bounded, because both schemas are now fully documented.

3. **Two persistence surfaces, not one.** Session JSONL (content-bearing, in-context) vs observability events (content-redacted, span-tree). sciagent's trace-recorder design (file 03's minimal JSONL schema) currently conflates these. A clean Pi mapping is: session JSONL → sciagent `user_prompt`/`tool_call`/`tool_result`/`assistant_text`; observability events → sciagent `subagent_start`/`subagent_stop`/timing metadata. Worth surfacing in ADR-005.

4. **The `.agents/` convention is `.agents/skills/`-only across both harnesses.** Pi confirms; Claude Code never adopted `.agents/`. The toolkit's `.agents/agents/` and `.agents/commands/` directories are pure sciagent convention, consumed by nothing else. Architecture doc §3 should drop the "Pi extension reading `.agents/agents/`" hand-wave entirely.

5. **No Task-tool analog forces a different distiller deployment.** On Pi, the distiller is a slash command or an out-of-process pipeline, not an in-conversation sub-agent. ADR-005's "sub-agent 3.4 trace-distiller" framing is Claude-Code-specific. Document explicitly that Pi gets a `pi.registerCommand("sciagent-distill", …)` shipped by a sciagent Pi extension.

## Negative findings retained, with refinement

- "Pi has no documented hook system." **Rescinded.** Pi has the most thoroughly documented hook system of any harness sciagent considers, with formal mutation semantics (`docs/.ref/pi/packages/agent/docs/hooks.md`).
- "Pi JSONL schema is undocumented." **Rescinded.** `docs/.ref/pi/packages/coding-agent/docs/session-format.md` is the spec, with TypeScript types in `docs/.ref/pi/packages/agent/src/types.ts`.
- "Pi `.agents/agents/` and `.agents/commands/` walking is an open question." **Resolved: no.** Pi walks `.agents/skills/` only. The dual-track convention is sciagent-only.
- "Pi has no Task-tool analog." **Confirmed.** Sub-agent / Task-tool patterns must be Claude-Code-only or built as Pi slash commands per harness.

Sources:
- `docs/.ref/pi/README.md`
- `docs/.ref/pi/packages/coding-agent/docs/extensions.md`
- `docs/.ref/pi/packages/coding-agent/docs/skills.md`
- `docs/.ref/pi/packages/coding-agent/docs/session-format.md`
- `docs/.ref/pi/packages/coding-agent/src/core/package-manager.ts:2278-2339`
- `docs/.ref/pi/packages/coding-agent/src/core/tools/`
- `docs/.ref/pi/packages/agent/docs/hooks.md`
- `docs/.ref/pi/packages/agent/docs/observability.md`
- `docs/.ref/pi/packages/agent/src/types.ts`
