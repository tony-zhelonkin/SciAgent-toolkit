# 01 — Trace recording in Claude Code

> **Outcome note (superseded).** Superseded by the three-verb toolkit
> architecture documented in [`docs/architecture.md`](../../architecture.md).
> Retained as a dated research record; the original body below is unchanged.

Use case for sciagent: ADR-005 requires recording sessions so a `trace-distiller` sub-agent can convert ad-hoc workflows into shipping skills. This file answers: *what trace surface does Claude Code already expose?*

---

## On-disk JSONL transcript — already comprehensive

Claude Code persists every session as an append-only JSONL file. Path pattern:

```
~/.claude/projects/<url-encoded-project-path>/<session-uuid>.jsonl
```

The project-path encoding swaps `/` for `-` (`/home/user/myapp` → `-home-user-myapp`). Each line is one JSON object — one event in the session DAG.

**Top-level envelope** (every line):

| field | meaning |
|---|---|
| `type` | `user`, `assistant`, `tool_result`, `system`, `summary`, `result`, `file-history-snapshot` |
| `uuid` | this entry's id |
| `parentUuid` | the entry this one descends from — enables branching DAG, not a flat list |
| `timestamp` | ISO-8601 UTC |
| `sessionId` | session-wide UUID |
| `cwd` | working directory when entry was written |
| `gitBranch`, `version` | runtime context |
| `message` | the actual content (shape depends on `type`) |

**Concrete excerpt** — assistant turn with thinking + tool use, then tool result:

```jsonl
{"type":"assistant","uuid":"3fa8...","parentUuid":"1a2b...","timestamp":"2026-05-22T09:14:32.441Z","sessionId":"abc","cwd":"/home/u/p","message":{"role":"assistant","model":"claude-opus-4-7","content":[{"type":"thinking","thinking":"I need to..."},{"type":"text","text":"Reading the file."},{"type":"tool_use","id":"toolu_01abc","name":"Read","input":{"file_path":"/x/y.py"}}],"usage":{"input_tokens":12840,"output_tokens":631,"cache_read_input_tokens":8200}}}
{"type":"tool_result","uuid":"7d2e...","parentUuid":"3fa8...","timestamp":"2026-05-22T09:14:33.102Z","sessionId":"abc","toolUseResult":{"tool_use_id":"toolu_01abc","content":"...","is_error":false}}
```

Tool calls correlate by `tool_use_id`. Content blocks inside assistant messages can be `thinking`, `text`, `tool_use`. Token usage is per-turn.

This format is already sufficient for a distiller — it has user prompts, assistant text, tool names + args + results, sub-agent boundaries (via `SubagentStart` / `SubagentStop` hook events written into the transcript). There is no need to "subscribe" if the consumer can read the file after the fact.

**Documentation status:** unofficial. The format is reverse-engineered in community write-ups [databunny.medium.com](https://databunny.medium.com/inside-claude-code-the-session-file-format-and-how-to-inspect-it-b9998e66d56b), [claude-dev.tools](https://claude-dev.tools/docs/jsonl-format). Anthropic does not publish a schema reference — it is observed-from-the-wild. This is a stability risk for ADR-005.

---

## Hook surface — for in-flight subscription

Hook events fire at specific lifecycle moments and receive JSON payloads on **stdin**. Every payload includes `session_id`, `transcript_path` (path to the JSONL above), `cwd`, `hook_event_name`. Source: [Claude Code Hooks reference](https://code.claude.com/docs/en/hooks).

Events relevant to trace capture:

| event | when | payload-specific fields |
|---|---|---|
| `SessionStart` | startup / resume / clear / compact | `source`, `model` |
| `UserPromptSubmit` | before Claude processes a user prompt | `prompt` |
| `UserPromptExpansion` | slash command → prompt | `command_name`, `command_args`, `prompt` |
| `PreToolUse` | tool args computed, before execution | `tool_name`, `tool_use_id`, `tool_input` |
| `PostToolUse` | tool succeeded | `tool_name`, `tool_input`, `tool_result` |
| `PostToolUseFailure` | tool failed | `tool_name`, `tool_input`, `error` |
| `PostToolBatch` | full parallel batch resolved | `tool_calls[]` |
| `SubagentStart` / `SubagentStop` | sub-agent boundaries | `agent_type`, `agent_id`, `result` (on stop) |
| `Stop` | assistant turn ends | `response` |
| `SessionEnd` | session terminates | `reason` |

Example `PreToolUse` payload for a Bash call:

```json
{
  "session_id": "abc123",
  "transcript_path": "/home/u/.claude/projects/-home-u-p/abc123.jsonl",
  "cwd": "/home/u/p",
  "permission_mode": "default",
  "hook_event_name": "PreToolUse",
  "tool_name": "Bash",
  "tool_use_id": "toolu_01...",
  "tool_input": {"command": "npm test", "description": "Run test suite"}
}
```

**Practical implication for ADR-005:** sciagent does not need to invent capture. It can register a `SessionStart` hook that opens an append-only file under `.sciagent/traces/<timestamp>.jsonl`, then a `Stop` / `SessionEnd` hook that closes it — or it can simply copy/symlink `transcript_path` after the session ends. The latter is parasitic and cheap.

---

## `/resume` — the transcript is the recovery format

Claude Code's `/resume` reads the JSONL transcript and reconstructs context. This means the transcript is already "complete enough" to drive a distiller — if Claude itself can reconstruct a session from it, a sub-agent can read it to identify invariant operations (constant commands, API endpoints) vs variant ones (user inputs, paths).

---

## OpenTelemetry / official export

No native OTel exporter in Claude Code as of May 2026. The JSONL is the only persistence surface. Anthropic does not publish a stable schema. Any external trace standard (OTel GenAI, OpenInference) has to be reached by writing a normalizer over the JSONL — see file 03.

---

## Bottom line for ADR-005

- **Capture cost:** near-zero. The JSONL exists whether sciagent records or not. A `sciagent record` command is mostly a convenience wrapper that pins `transcript_path` to a known location and tags it with a user-supplied label.
- **Schema risk:** the format is undocumented and could change without notice. A distiller written against today's shape might break on a future Claude Code release.
- **Cross-harness portability:** zero. This is Claude-specific. The Pi story is in file 02.

Sources:
- [Hooks reference — Claude Code docs](https://code.claude.com/docs/en/hooks)
- [Inside Claude Code: The Session File Format — Yi Huang, Medium](https://databunny.medium.com/inside-claude-code-the-session-file-format-and-how-to-inspect-it-b9998e66d56b)
- [Claude Code JSONL transcript format explained — claude-dev.tools](https://claude-dev.tools/docs/jsonl-format)
- [withLinda/claude-JSONL-browser — community reader](https://github.com/withLinda/claude-JSONL-browser)
