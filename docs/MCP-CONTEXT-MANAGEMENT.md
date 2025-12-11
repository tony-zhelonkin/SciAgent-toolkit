# MCP Context Management Guide

**Last Updated:** 2025-12-10
**Status:** Best practices for managing MCP tool context consumption across Claude, Gemini, and Codex.

---

## The Problem: Context Window Bloat

MCP tools consume context tokens **before your conversation even starts**. Each tool definition includes its name, description, parameters, and schema - this can add up quickly.

**Claude Opus 4.5 has a 200k token window.** Loading all ToolUniverse tools (600+) exceeds this limit immediately.

| MCP Server | Approx. Context Cost | Use Case |
|------------|---------------------|----------|
| Sequential Thinking | ~1.6k tokens | Structured reasoning |
| Context7 | ~1.8k tokens | Library documentation |
| PAL (22 tools) | ~22k tokens | Multi-model collaboration |
| Serena | ~50k tokens | Semantic code search |
| ToolUniverse (ALL 600 tools) | ~500k+ tokens | Scientific research |
| ToolUniverse (14 curated) | ~20k tokens | Targeted research |

---

## Solution: Profile-Based MCP Management

Our unified profile switcher configures **Claude**, **Gemini**, and **Codex** simultaneously, ensuring a consistent environment or specialized hybrid workflows.

### Profile Location

Profiles are stored in the SciAgent-toolkit repository:

```
toolkits/SciAgent-toolkit/
├── templates/mcp-profiles/
│   ├── minimal.mcp.json         # ~3k tokens   - Default coding
│   ├── coding.mcp.json          # ~25k tokens  - + PAL
│   ├── codebase.mcp.json        # ~75k tokens  - + PAL + Serena
│   ├── research-lite.mcp.json   # ~30k tokens  - + ToolUniverse (6 tools)
│   ├── research-full.mcp.json   # ~50k tokens  - + ToolUniverse (14 tools)
│   ├── hybrid-research.mcp.json # Specialized  - Claude=Lite, Gemini=Heavy
│   └── full.mcp.json            # ~100k tokens - PAL + ToolUniverse + Serena
└── scripts/
    └── switch-mcp-profile.sh    # Unified Profile Switcher
```

### Using the Profile Switcher

```bash
# Switch to a profile (Apply to all CLIs)
./toolkits/SciAgent-toolkit/scripts/switch-mcp-profile.sh minimal
./toolkits/SciAgent-toolkit/scripts/switch-mcp-profile.sh coding
./toolkits/SciAgent-toolkit/scripts/switch-mcp-profile.sh hybrid-research

# Restart your agent (Claude, Gemini, or Codex) to apply changes.
```

---

## Profile Descriptions & Workflow Mapping

The switcher intelligently maps the selected profile to the specific needs of each CLI agent.

### 1. `minimal` (~3k tokens)
**Best for:** Quick questions, small scripts, maximum context availability.
- **Claude/Codex:** Sequential Thinking, Context7.
- **Gemini:** Same (via "coding" template).

### 2. `coding` (~25k tokens) - **Default**
**Best for:** General software engineering, debugging, code review.
- **Claude/Codex:** Sequential Thinking, Context7, **PAL** (for multi-model consensus).
- **Gemini:** Same.
- **Workflow:** Use PAL to ask other models for opinions (`Use pal to ask gemini-pro...`).

### 3. `codebase` (~75k tokens)
**Best for:** Understanding large, unfamiliar codebases.
- **Claude/Codex:** Adds **Serena** (Semantic Code Search).
- **Gemini:** Adds Serena.
- **Workflow:** `Use serena to find symbol definitions...`

### 4. `research-lite` / `research-full` (~30k - 50k tokens)
**Best for:** Literature review, drug discovery, protein analysis.
- **Claude/Codex:** Adds **ToolUniverse** (Curated set of 6 or 14 tools).
- **Gemini:** Adds ToolUniverse (Curated set + Google Search).
- **Note:** Gemini's 2M context window handles these tools easily.

### 5. `hybrid-research` (Specialized)
**Best for:** **The "Orchestrator" Workflow.**
Keep Claude lightweight for coding while leveraging Gemini for heavy research.
- **Claude:** `coding` profile (Lightweight, PAL enabled).
- **Gemini:** `research` profile (ToolUniverse enabled).
- **Codex:** `coding` profile (Lightweight).
- **Workflow:**
    1.  Work in Claude for coding.
    2.  Need research? Use PAL to spawn Gemini CLI.
    3.  `Use pal to run gemini "Search for recent papers on CRISPR..."`
    4.  Gemini (running in background) uses its ToolUniverse tools and returns the summary to Claude.

---

## CLI-Specific Details

### Claude Code (`.mcp.json`)
- Uses stdio transport for all servers.
- PAL is configured to allow spawning sub-agents.

### Gemini CLI (`.gemini/settings.json`)
- Automatically configured by the switcher.
- Uses **Compact Mode** for ToolUniverse to save tokens (exposing only 4-5 core tools like `search_tools` and `execute_tool`).
- API Keys (`GEMINI_API_KEY`, `OPENAI_API_KEY`) are injected safely from your environment.

### Codex CLI (`~/.codex/config.toml`)
- Configured to match Claude's profile.
- ToolUniverse is dynamically enabled/disabled based on the profile.

---

## ToolUniverse Configuration Best Practices

**NEVER load all 600 tools.** Always specify what you need via `--include-tools`.

The profile switcher handles this automatically:

**Literature Research Set:**
`EuropePMC_search_articles`, `openalex_literature_search`, `pubmed_search`, `semantic_scholar_search`

**Drug Discovery Set:**
`ChEMBL_search_similar_molecules`, `ChEMBL_get_molecule`, `DrugBank_search`, `FDA_adverse_events`

---

## Troubleshooting

### PAL "Failed to connect"
Ensure you have API keys in `.devcontainer/.env`:
```bash
GEMINI_API_KEY=your-key
OPENAI_API_KEY=your-key
```
Then re-run the profile switcher.

### Gemini CLI not showing tools?
Verify the configuration:
```bash
cat .gemini/settings.json
```
Ensure `tooluniverse` is present in `mcpServers`.

### "Thinking not supported" Error in Gemini
Some Gemini models (like `flash-thinking`) do not support the `sequential-thinking` tool because they have native thinking capabilities. Use `gemini-2.0-flash-001` or similar standard models when using MCP thinking tools.

---

## References

- [Anthropic: Advanced Tool Use](https://www.anthropic.com/engineering/advanced-tool-use)
- [ToolUniverse Documentation](https://github.com/greedyai/tool-universe)
- [PAL MCP Server](https://github.com/BeehiveInnovations/pal-mcp-server)