# Hybrid Research Workflow: Claude + Gemini

This workflow demonstrates how to use the "Hybrid Research" profile, where Claude Code handles the coding/orchestration while delegating heavy research tasks to Gemini CLI (leveraging its 2M token context window).

## Prerequisites

1. **Install SciAgent-toolkit**:
   ```bash
   ./scripts/setup-ai.sh
   ```

2. **Configure API Keys**:
   Ensure `.env` has both `GEMINI_API_KEY` and `OPENAI_API_KEY`.

3. **Switch to Hybrid Profile**:
   ```bash
   ./scripts/switch-mcp-profile.sh hybrid-research
   ```

## Workflow Steps

### 1. Start in Claude Code

Open your terminal and start Claude:

```bash
claude
```

Claude is now configured in "Lightweight" mode (fast, low context usage) but has access to the **PAL** MCP server, which can talk to Gemini.

### 2. The Research Task

Let's say you need to write a Python script that analyzes the latest literature on a specific gene, but you don't want to clutter Claude's context with hundreds of search results.

**Prompt to Claude:**
> "I need to write a Python script to analyze recent papers on the BRCA1 gene.
> First, use PAL to ask Gemini to search PubMed and EuropePMC for papers from 2024-2025.
> Ask Gemini to summarize the key methodologies used in those papers.
> Then, use that summary to write the Python analysis script."

### 3. What Happens Behind the Scenes

1. **Claude** receives your request.
2. **Claude** calls `pal.chat` (or `pal.thinkdeep`) targeting `gemini-pro`.
3. **PAL** forwards the request ("Search PubMed...") to **Gemini**.
4. **Gemini** (running in the background with the "Research" profile) executes:
   - `ToolUniverse.pubmed_search(...)`
   - `ToolUniverse.EuropePMC_search_articles(...)`
5. **Gemini** processes the raw search results (potentially hundreds of lines) within its large context window.
6. **Gemini** generates a concise summary of methodologies.
7. **PAL** returns just this summary to **Claude**.
8. **Claude** uses the summary to write the Python script.

### 4. Benefits

- **Cost Efficient:** Claude doesn't process the raw search tokens.
- **Context Efficient:** Claude's context window stays clean for code generation.
- **Best Tool for the Job:** Gemini's long context is perfect for reading papers; Claude's reasoning is great for coding.

## Example: Drug Discovery

**Prompt:**
> "Use PAL to have Gemini find FDA-approved drugs for Alzheimer's that target Beta-secretase 1.
> Have Gemini list their chemical structures (SMILES).
> Once you get that list, write a RDKit script to visualize them."

## Troubleshooting

- **"PAL failed to connect"**: Check `GEMINI_API_KEY` in `.env`.
- **"Gemini cannot find tools"**: Ensure `switch-mcp-profile.sh hybrid-research` was run successfully. Check `.gemini/settings.json`.