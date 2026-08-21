# CLI capability matrix (source-grounded, 2026-07-02)

Ground truth = actually-installed binaries/packages on this machine, cross-checked against official
docs. Where installed source and web disagree, source wins.

Installed: pi `@earendil-works/pi-coding-agent@0.79.4` (npm global, with `docs/`+`examples/`),
`@aliou/pi-guardrails@0.13.1`, `context-mode@1.0.151`; agy `~/.local/bin/agy` (183 MB Go binary);
codex config at `~/.codex/` (binary off PATH but config real); claude `~/.local/bin/claude` v2.1.198.
opencode not installed (docs-only).

## Context-file support (the portable substrate)

| CLI | `AGENTS.md` native | Also reads | Shim for canonical AGENTS.md |
|-----|:---:|---|---|
| codex | ✅ canonical | configurable fallbacks (`project_doc_fallback_filenames`) | none |
| pi | ✅ (`AGENTS.md` **or** `CLAUDE.md`) | `CLAUDE.md`; `.pi/SYSTEM.md` overrides system prompt | none |
| agy | ✅ (`AGENTS.md` **and** `GEMINI.md`) | `GEMINI.md`; writes rules back into AGENTS.md | none |
| opencode | ✅ | `instructions[]` globs | none |
| claude | ❌ (reads `CLAUDE.md`) | `@path` imports | `CLAUDE.md` = `@AGENTS.md` (already used) |

**Conclusion: `AGENTS.md` + `.agents/skills/` is the portable substrate. Only Claude needs a one-line
shim it already has.**

## Config location & format

| CLI | Global config | Project config | Format | Merge/precedence |
|-----|---|---|---|---|
| codex | `~/.codex/config.toml` (`$CODEX_HOME`) | `.codex/config.toml` (trusted only) | TOML | CLI > `-c` > profile > project > user > `/etc/codex` |
| pi | `~/.pi/agent/settings.json` | `.pi/settings.json` | JSON (deep-merge) | project overrides global |
| agy | `~/.gemini/antigravity-cli/` + `mcp_config.json` | workspace `AGENTS.md` | JSON | Global vs Workspace Customizations Root; honors `XDG_CONFIG_HOME` |
| opencode | `~/.config/opencode/opencode.json` | `opencode.json` / `.opencode/` | JSON/JSONC | project > `OPENCODE_CONFIG` > global (merged) |
| claude | `~/.claude/settings.json` | `.claude/settings.json` + `.local.json` | JSON | user < project < local < managed |

## Sub-agents / skills / hooks / MCP

| CLI | Sub-agents | Skills convention | Hooks | MCP |
|-----|---|---|---|---|
| codex | none first-class (parallel via multiple `codex exec`); `exec fork` requested | `~/.codex/skills/` (agentskills.io) | execpolicy `.rules` | ✅ `[mcp_servers.*]` |
| pi | **example extension** (single/parallel≤8,4-conc/chain via `pi -p --mode json` children) | `.agents/skills/`, reads `~/.claude/skills`,`~/.codex/skills` if listed | ✅ `pi.on("tool_call")` blocking + full event bus | ✅ |
| agy | RPC `define_subagent`/`invoke_subagent` | `.agents/skills/<n>/SKILL.md` + `~/.gemini/antigravity-cli/skills/` | ✅ PreToolUse/json-hooks | ✅ `mcp_config.json` |
| opencode | **first-class** primary/subagent, Task tool + `@mention` | `.opencode/skills/` | `.opencode/plugins/` | ✅ `mcp` |
| claude | Task tool + `.claude/agents/*.md` | `.claude/skills/`, `~/.claude/skills/` | ✅ settings.json hooks | ✅ |

## Headless invocation

- **codex:** `codex exec [-m MODEL] -s SANDBOX [-o out] [--json] [--skip-git-repo-check] "PROMPT"`
- **pi:** `pi -p "PROMPT"` (merges stdin); `--mode json`/`rpc`; `--model provider/id:thinking`
- **agy:** `agy -p "PROMPT"` (**no `--model`** — auto-selects); `-c` continue; `--sandbox`
- **opencode:** `opencode run "PROMPT" -m provider/model [--agent] [--continue] [--session]`
- **claude:** `claude -p "PROMPT" --output-format json|stream-json --model ...`

## pi extension API (from installed `docs/extensions.md` + `examples/` + `@aliou/pi-guardrails`)

- Load: `~/.pi/agent/extensions/`, `.pi/extensions/`, package `"pi":{"extensions":[...]}`, `-e`.
- Entry: `export default function (pi: ExtensionAPI) { ... }` (async ok).
- API: `pi.registerTool`, `pi.registerCommand`, `pi.on(event,h)`, `pi.events.on/emit`,
  `pi.registerProvider`; helpers `defineTool`, `getAgentDir`, `parseFrontmatter`, `withFileMutationQueue`.
- Events: `session_start/shutdown`, `before_agent_start`, `agent_start/end`, `turn_start/end`,
  `message_*`, `before_provider_request/after_provider_response`, `tool_call` (blockable), `tool_result`,
  `input`, `user_bash`, `model_select`.
- Sub-agent example: `examples/extensions/subagent/` — modes Single/Parallel/Chain; agents from
  `~/.pi/agent/agents/*.md` + `.pi/agents/*.md`, frontmatter `name/description/tools/model`,
  `agentScope: user|project|both`.

## Corrections to prior-art

- `ai-research/02` "pi has no Task/sub-agent primitive" — **outdated**; the installed package ships a
  full sub-agent example extension.
- `ai-research/02` "pi never adopted AGENTS.md" — **wrong**; pi reads `AGENTS.md` or `CLAUDE.md` natively.
- There is no Anthropic "Antigravity Skills paper" (confirmed prior-art `10:228`); the real reference is
  Google DeepMind's `workflow_skill_creator`.
