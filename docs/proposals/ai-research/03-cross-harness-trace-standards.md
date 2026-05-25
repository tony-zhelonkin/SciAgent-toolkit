> **Updated** 2026-05-24 with grounded findings from local repo at `docs/.ref/pi/`.

# 03 — Cross-harness trace standards

The single most consequential research finding for ADR-005. The question: should sciagent define its own JSONL schema, adopt OpenTelemetry GenAI, adopt OpenInference, or hybrid?

---

## OpenTelemetry GenAI semantic conventions

Status: **Development / experimental**, not stable. Migration to stable is gated on an `OTEL_SEMCONV_STABILITY_OPT_IN` env var with values like `gen_ai_latest_experimental`. Source: [OpenTelemetry semconv — gen-ai](https://opentelemetry.io/docs/specs/semconv/gen-ai/).

Three span types relevant to sciagent:

- **Generative AI client spans** — one LLM call. Attributes: `gen_ai.system`, `gen_ai.operation.name`, `gen_ai.request.model`, `gen_ai.response.finish_reasons`, `gen_ai.input.messages`, `gen_ai.output.messages`. Tool calls captured inside `gen_ai.output.messages` and via `gen_ai.tool.definitions`. [gen-ai-spans](https://opentelemetry.io/docs/specs/semconv/gen-ai/gen-ai-spans/)
- **Agent spans** — `create_agent`, `invoke_agent`, `invoke_workflow`. Attributes: `gen_ai.agent.name`, `gen_ai.agent.id`, `gen_ai.agent.version`, `gen_ai.agent.description`. Sub-agent invocation = nested `invoke_agent` spans under an `invoke_workflow` parent. [gen-ai-agent-spans](https://opentelemetry.io/docs/specs/semconv/gen-ai/gen-ai-agent-spans/)
- **Tool execution spans** — `execute_tool` with tool-name + arguments attributes.

Provider-specific extensions exist for Anthropic, OpenAI, AWS Bedrock, Azure AI Inference.

Field map to sciagent needs:

| sciagent need | OTel GenAI attribute |
|---|---|
| user prompt | `gen_ai.input.messages[].content` (role=user) |
| agent response | `gen_ai.output.messages[].content` (role=assistant) |
| tool call name | `tool_call.function.name` inside output messages |
| tool call args | `tool_call.function.arguments` |
| tool result | next event in message history, role=tool |
| sub-agent invocation | child `invoke_agent` span under parent `invoke_workflow` |
| timestamp | span start/end times (native OTel) |

Covers everything sciagent needs. Cost: it's spans, not files. Requires an exporter — OTLP gRPC/HTTP — and a collector. Not a JSONL.

## OpenInference (Arize)

Status: active, mature spec. Built on top of OpenTelemetry — i.e., it's an OpenTelemetry-attribute convention, not a competing transport. Source: [Arize-ai/openinference — semantic conventions](https://arize-ai.github.io/openinference/spec/semantic_conventions.html).

**Span kinds**: `LLM`, `EMBEDDING`, `CHAIN`, `RETRIEVER`, `RERANKER`, `TOOL`, `AGENT`, `GUARDRAIL`, `EVALUATOR`, `PROMPT`. Richer typology than OTel GenAI's four agent operations.

**Attributes relevant to sciagent**:
- input: `input.value`, `llm.input_messages`, `message.content`, `message.role`
- output: `llm.output_messages`, `llm.finish_reason`, `llm.token_count.*`
- tool: `tool.name`, `tool.description`, `tool.json_schema`, `tool_call.id`, `tool_call.function.name`, `tool_call.function.arguments`
- agent hierarchy: `agent.name`, parent-child via `graph.node.parent_id`

Natively supported by [arize-phoenix](https://github.com/Arize-ai/phoenix) (local, open-source observability UI), works with any OTel backend.

## LangSmith trace format

Run-based (analogous to OTel spans). Top-level fields: `id`, `trace_id`, `parent_run_id`, `run_type`, `start_time`, `end_time`, `inputs`, `outputs`, `error`, `name`, `tags`, `total_tokens`, `total_cost`. `run_type` values: `llm`, `chain`, `tool`, `retriever`, `embedding`, `prompt`, `parser`. Hierarchy encoded in a `dotted_order` field — sortable key of the form `<timestamp>Z<uuid>` separated by dots. Source: [Run (span) data format — LangChain docs](https://docs.langchain.com/langsmith/run-data-format).

Proprietary platform; format is documented but is a LangSmith-API contract first, a portable schema second. Adoption outside LangChain ecosystem is limited.

## Langfuse trace format

Open-source, OTel-built. Top-level: `trace` (single request/operation, container). Within a trace: `observations` typed as `event` | `span` | `generation` | `agent` | `tool` | `chain` | `retriever` | `evaluator`. Source: [Langfuse data model](https://langfuse.com/docs/observability/data-model).

Better-typed than LangSmith for agent workflows. Has a self-hostable platform.

## Phoenix / Helicone

Both consume OpenInference (Phoenix is Arize's open-source UI; Helicone is a separate proxy-based observability tool). Mentioning for completeness — neither defines a competing schema; they're sinks.

---

## Recommendation: hybrid — sciagent-native JSONL on disk + optional OpenInference exporter

### Why not pure OTel GenAI

- Still in development. Tying ADR-005 to an evolving spec invites breakage. `OTEL_SEMCONV_STABILITY_OPT_IN` exists precisely because the gen_ai surface is moving.
- Requires a collector/exporter pipeline. sciagent is a per-project bash CLI. Standing up an OTLP collector to record a session is an order of magnitude more infrastructure than `cat >> .sciagent/traces/<ts>.jsonl`.
- The trace-distiller sub-agent (ADR-005 §3.4) operates on a file, not a stream. It needs random access by line, not a span graph.

### Why not pure OpenInference

- Same infrastructure cost as OTel GenAI (it *is* OTel-based).
- Vendor-leaning even if the spec is open.
- Distillation needs message-by-message ordering, which OpenInference encodes via OTel span timestamps — fine in principle, awkward for a bash + Python distiller.

### Why not pure LangSmith / Langfuse

- Both are tied to product platforms. Even where self-hostable (Langfuse), they're a heavier dependency than sciagent should take on.

### Why sciagent-native JSONL on disk

- Claude Code and Pi both already write JSONL. A sciagent-native schema can be a **subset/normalisation** of those — no new persistence mechanism, just a fixed key layout for the distiller to consume.
- Distillation is the primary consumer. Distiller reads file → identifies invariants → emits SKILL.md draft. No queries, no aggregations, no UI.
- ADR-002's `requires:` graph and ADR-006's bench reports already live as files in the repo. JSONL traces fit the same idiom.

### Update — Pi schema is documented, not undocumented

Earlier framing of the Pi side as "format is undocumented, blocks cross-harness work" is wrong. Pi ships a formal `version: 3` JSONL spec (`docs/.ref/pi/packages/coding-agent/docs/session-format.md`) with typed entries (`UserMessage`, `AssistantMessage`, `ToolResultMessage`, `BashExecutionMessage`, `CustomMessage`, `BranchSummaryMessage`, `CompactionSummaryMessage`), 8-char `id` + `parentId` tree linkage, and TypeScript source-of-truth at `docs/.ref/pi/packages/agent/src/types.ts`. The Claude-Code-vs-Pi normalizer becomes a documented field-mapping exercise, not a reverse-engineering project. Both ends of the cross-harness pipeline are now known.

Pi additionally exposes a separate **observability** event stream (`docs/.ref/pi/packages/agent/docs/observability.md`) for content-redacted timings/spans — distinct from session JSONL. The sciagent-native schema above conflates content events (`user_prompt`, `tool_call`, `tool_result`) and structural events (`session_start`, `subagent_start`, `subagent_stop`). Pi's split suggests a cleaner sciagent design: one JSONL for content (consumed by the distiller), one event log for spans (consumed by the bench scorer and optional OTel/OpenInference exporter). This is a tractable refinement, not a redesign.

### Why also ship an optional OpenInference exporter

- The day someone wants Phoenix / Arize / Langfuse UI over sciagent traces, the conversion is mechanical. Building this exporter is a 1-2 day job once the native JSONL is stable.
- Future-proofs the architecture without committing to it upfront. If OTel GenAI becomes stable and ubiquitous in 2027, swap the exporter target — the native format does not change.
- This is the same pattern Langfuse adopted: native model + OTel ingest.

### Minimal sciagent JSONL schema (proposed)

```jsonl
{"v":1,"event":"session_start","session_id":"...","harness":"claude-code","timestamp":"...","cwd":"..."}
{"v":1,"event":"user_prompt","session_id":"...","timestamp":"...","prompt":"..."}
{"v":1,"event":"tool_call","session_id":"...","timestamp":"...","tool":"Bash","input":{...},"call_id":"..."}
{"v":1,"event":"tool_result","session_id":"...","timestamp":"...","call_id":"...","output":"...","error":null}
{"v":1,"event":"subagent_start","session_id":"...","timestamp":"...","agent":"trace-distiller","agent_id":"..."}
{"v":1,"event":"subagent_stop","session_id":"...","timestamp":"...","agent_id":"...","result":"..."}
{"v":1,"event":"assistant_text","session_id":"...","timestamp":"...","text":"..."}
{"v":1,"event":"session_end","session_id":"...","timestamp":"...","reason":"..."}
```

Versioned (`v`), flat, one event per line, no nesting. Claude Code → sciagent and Pi → sciagent normalizers fit in ~150 lines of Python each. The trace-distiller reads `events == "tool_call"` and `tool` to infer invariant operations, `events == "user_prompt"` for the parameter surface.

## Tradeoff table

| dimension | sciagent-native JSONL | OTel GenAI | OpenInference | LangSmith | Langfuse |
|---|---|---|---|---|---|
| schema stability | controlled by us | development | mature, on top of dev OTel | proprietary, stable | OSS, stable |
| ecosystem | none yet | growing | Phoenix + many | LangChain | self-host + Langfuse Cloud |
| infrastructure to record | zero | OTLP collector | OTLP collector | platform account | platform/self-host |
| distillation simplicity | trivial (line by line) | medium (span tree walk) | medium (span tree walk) | medium | medium |
| cross-harness | yes via normalizer | yes if both export | yes if both export | no | no |
| vendor lock | none | OTel ecosystem | Arize-leaning | LangChain | Langfuse |

**Position:** sciagent-native JSONL is the load-bearing surface. OpenInference exporter is an opt-in convenience. Do not adopt OTel GenAI or LangSmith as the primary format.

Sources:
- [OpenTelemetry GenAI semantic conventions](https://opentelemetry.io/docs/specs/semconv/gen-ai/)
- [GenAI client spans](https://opentelemetry.io/docs/specs/semconv/gen-ai/gen-ai-spans/)
- [GenAI agent spans](https://opentelemetry.io/docs/specs/semconv/gen-ai/gen-ai-agent-spans/)
- [OpenInference semantic conventions](https://arize-ai.github.io/openinference/spec/semantic_conventions.html)
- [Arize-ai/openinference repo](https://github.com/Arize-ai/openinference)
- [LangSmith run data format](https://docs.langchain.com/langsmith/run-data-format)
- [Langfuse data model](https://langfuse.com/docs/observability/data-model)
- [Inside the LLM Call — OpenTelemetry blog 2026](https://opentelemetry.io/blog/2026/genai-observability/)
