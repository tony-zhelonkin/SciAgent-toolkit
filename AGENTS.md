# AI Agents in SciAgent-toolkit

This toolkit integrates several powerful AI agents and MCP servers to assist with scientific research.

## PAL (Personal AI Lead)

PAL is a primary orchestrator agent designed to help with collaboration, planning, and code analysis. It provides a suite of tools that act as specialized sub-agents.

### Capabilities

- **Collaboration & Planning**:
    - `chat`: Brainstorm ideas, get second opinions, and validate approaches.
    - `thinkdeep`: Extended reasoning and analysis for complex problems.
    - `planner`: Break down complex projects into structured, actionable plans.
    - `consensus`: Get expert opinions from multiple AI models (if configured).

- **Code Analysis & Quality**:
    - `debug`: Systematic investigation and root cause analysis.
    - `codereview`: Professional reviews with actionable feedback.
    - `analyze`: Understand architecture and patterns (disabled by default).

### Usage

PAL tools are available directly in Claude Code after installation.

**Example Prompts:**
- "Use pal planner to create a roadmap for this analysis."
- "Ask pal thinkdeep to evaluate the potential flaws in this experimental design."
- "Use pal codereview to check this R script for best practices."

## Domain-Specific Agents

The toolkit also includes specialized agents defined in `agents/`:

### Bioinformatics Research Librarian
- **Expertise**: Finding bioinformatics tools, documentation, and resources.
- **Tools**: Web search, documentation lookups.
- **Use Case**: "Find the best tool for scATAC-seq peak calling."

### RNA-seq Methods Writer
- **Expertise**: Generating publication-ready Methods sections from code.
- **Tools**: Code analysis, file reading.
- **Use Case**: "Write a methods section based on the analysis in `./scripts`."

## Configuration

PAL is configured automatically by `./scripts/setup_mcp_infrastructure.sh`.
Key settings can be adjusted in the generated `.mcp.json` file or by modifying the setup script.
