# MCP Context Management Guide

**Last Updated:** 2025-12-10
**Status:** Best practices for managing MCP tool context consumption

---

## The Problem: Context Window Bloat

MCP tools consume context tokens **before your conversation even starts**. Each tool definition includes its name, description, parameters, and schema - this can add up quickly:

| MCP Server | Approx. Context Cost | Use Case |
|------------|---------------------|----------|
| Sequential Thinking | ~1.6k tokens | Structured reasoning |
| Context7 | ~1.8k tokens | Library documentation |
| PAL (22 tools) | ~22k tokens | Multi-model collaboration |
| Serena | ~50k tokens | Semantic code search |
| ToolUniverse (ALL 600 tools) | ~500k+ tokens | Scientific research |
| ToolUniverse (14 curated) | ~20k tokens | Targeted research |

**Claude Opus 4.5 has a 200k token window.** Loading all ToolUniverse tools alone exceeds this limit!

**CRITICAL: Never load all 600 ToolUniverse tools.** Always use `--include-tools` to specify only what you need.

---

## Solution: Profile-Based MCP Management

### Profile Location

Profiles are stored in the SciAgent-toolkit repository:

```
toolkits/SciAgent-toolkit/
├── templates/mcp-profiles/
│   ├── minimal.mcp.json       # ~3k tokens   - Default coding
│   ├── coding.mcp.json        # ~25k tokens  - + PAL
│   ├── codebase.mcp.json      # ~75k tokens  - + PAL + Serena (large codebase work)
│   ├── research-lite.mcp.json # ~30k tokens  - + ToolUniverse (6 tools)
│   ├── research-full.mcp.json # ~50k tokens  - + ToolUniverse (14 tools)
│   └── full.mcp.json          # ~100k tokens - PAL + ToolUniverse + Serena
└── scripts/
    └── switch-mcp-profile.sh  # Profile switcher
```

### Using the Profile Switcher

```bash
# Show available profiles
./toolkits/SciAgent-toolkit/scripts/switch-mcp-profile.sh

# Switch to a profile
./toolkits/SciAgent-toolkit/scripts/switch-mcp-profile.sh minimal
./toolkits/SciAgent-toolkit/scripts/switch-mcp-profile.sh codebase
./toolkits/SciAgent-toolkit/scripts/switch-mcp-profile.sh research-lite

# Restart Claude Code to apply
exit && claude
```

---

## Profile Descriptions

### `minimal` (~3k tokens) - Default
**Servers:** Sequential Thinking, Context7
**Use for:**
- Standard coding sessions
- Quick questions and debugging
- Maximum free context for conversation

### `coding` (~25k tokens)
**Servers:** Sequential Thinking, Context7, PAL
**Use for:**
- Code review with multi-model consensus
- Planning sessions with external AI collaboration
- Debugging with multiple perspectives (Gemini, GPT, etc.)

### `codebase` (~75k tokens)
**Servers:** Sequential Thinking, Context7, PAL, Serena
**Use for:**
- Large codebase exploration and navigation
- Complex refactoring across many files
- Semantic code search and symbol analysis
- Understanding unfamiliar codebases

### `research-lite` (~30k tokens)
**Servers:** Sequential Thinking, Context7, ToolUniverse (6 tools)
**Tools:** EuropePMC, OpenAlex, ChEMBL (molecule/search), UniProt, ClinicalTrials
**Use for:**
- Quick literature searches
- Drug/protein lookups
- When you know exactly which tools you need

### `research-full` (~50k tokens)
**Servers:** Sequential Thinking, Context7, ToolUniverse (14 tools + SummarizationHook)
**Tools:** Full research set with summarization for long outputs
**Use for:**
- Broader research sessions
- When you might need multiple databases
- Long-running research tasks

### `full` (~100k tokens)
**Servers:** Sequential Thinking, Context7, PAL, ToolUniverse (14 tools), Serena
**Use for:**
- Combined code + research sessions
- Maximum capability (but uses 50% of context)
- Use sparingly - less room for conversation

---

## ToolUniverse Tool Filtering

### The Key Feature: `--include-tools`

**NEVER load all 600 tools.** Always specify what you need:

```bash
# BAD: All tools (~500k+ tokens) - NEVER DO THIS
uv run tooluniverse-smcp-stdio

# GOOD: Curated tools (~20k tokens)
uv run tooluniverse-smcp-stdio \
  --include-tools EuropePMC_search_articles,ChEMBL_search_similar_molecules,search_clinical_trials
```

### Recommended Tool Sets by Workflow

**Literature Research:**
```
--include-tools EuropePMC_search_articles,openalex_literature_search,pubmed_search,semantic_scholar_search
```

**Drug Discovery:**
```
--include-tools ChEMBL_search_similar_molecules,ChEMBL_get_molecule,DrugBank_search,FDA_adverse_events
```

**Genomics/Proteomics:**
```
--include-tools UniProt_search_proteins,UniProt_get_protein,protein_interaction_search
```

**Clinical Research:**
```
--include-tools search_clinical_trials,get_clinical_trial_details,FDA_drug_labels
```

---

## Context Monitoring

Always check context usage:

```bash
# Inside Claude Code:
/context
```

### Healthy Usage by Profile

| Profile | Expected MCP Tokens | Free Space |
|---------|--------------------:|------------|
| minimal | ~3k (1.5%) | ~150k (75%) |
| coding | ~25k (12.5%) | ~130k (65%) |
| codebase | ~75k (37.5%) | ~80k (40%) |
| research-lite | ~30k (15%) | ~125k (62%) |
| research-full | ~50k (25%) | ~105k (52%) |
| full | ~100k (50%) | ~55k (27%) |

### Warning Signs

- MCP tools > 100k tokens: Switch to lighter profile
- MCP tools > 150k tokens: Conversation will be severely limited
- MCP tools > 200k tokens: **Unusable** - context overflow

---

## Guardrails and Best Practices

### DO:
- Start sessions with `minimal` profile by default
- Use `/context` to monitor token usage
- **Always** use `--include-tools` with ToolUniverse
- Use `codebase` profile for large codebase exploration
- Switch profiles based on task needs

### DON'T:
- ❌ Load all 600 ToolUniverse tools
- ❌ Use `full` profile for simple coding tasks
- ❌ Ignore context warnings
- ❌ Keep heavy profiles active when not needed

### Profile Selection Flowchart

```
Starting new session?
    │
    ├── Standard coding/debugging? → minimal
    │
    ├── Need multi-model AI perspectives? → coding
    │
    ├── Large codebase exploration? → codebase
    │
    ├── Literature/drug/protein research?
    │       │
    │       ├── Know which tools? → research-lite
    │       │
    │       └── Exploratory? → research-full
    │
    └── Combined code + research? → full (use sparingly)
```

---

## Future Improvements

### Anthropic's Tool Search Tool (Coming to Claude Code)
- Will allow `defer_loading: true` on tools
- Tools discoverable but not loaded until needed
- 85% token reduction possible
- Currently API-only beta

### Subagent-Only MCP Tools (Feature Request #6915)
- Most-requested feature (144+ upvotes)
- Would allow `subagentOnly: true` flag
- Tools available to subagents, not main agent
- Vote at: https://github.com/anthropics/claude-code/issues/6915

---

## Troubleshooting

### PAL shows "Failed to connect"
PAL requires at least one AI provider API key. Add to `.devcontainer/.env`:
```bash
GEMINI_API_KEY=your-key    # https://aistudio.google.com/apikey
# OR
OPENAI_API_KEY=your-key    # https://platform.openai.com/api-keys
```
Then: `./toolkits/SciAgent-toolkit/scripts/configure_mcp_servers.sh --force`

### Context exceeds 200k immediately
1. Run `/context` to see breakdown
2. Identify which MCP servers consume most
3. Switch to `minimal` profile
4. Use `--include-tools` for ToolUniverse

### Profile switcher can't find ToolUniverse
The switcher searches for `tooluniverse-env/` in several locations. If not found:
```bash
# Run AI setup first
./toolkits/SciAgent-toolkit/scripts/setup-ai.sh
```

---

## References

- [Anthropic: Advanced Tool Use](https://www.anthropic.com/engineering/advanced-tool-use)
- [Optimising MCP Context Usage](https://scottspence.com/posts/optimising-mcp-server-context-usage-in-claude-code)
- [MCP Context Bloat Problem](https://medium.com/@yashfortunate/understanding-mcp-the-context-bloat-problem-in-ai-agents-and-the-risks-of-mcp-7a5046600a4b)
- [McPick CLI Tool](https://github.com/spences10/mcpick)
- [Feature Request: Subagent-Only MCP Tools](https://github.com/anthropics/claude-code/issues/6915)
- [ToolUniverse Documentation](https://github.com/greedyai/tool-universe)
