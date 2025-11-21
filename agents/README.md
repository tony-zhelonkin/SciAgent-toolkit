# Claude Agents for Bioinformatics

Pre-configured Claude agents optimized for bioinformatics and scientific research workflows.

## Available Agents

### Bioinformatics Research Librarian
**File:** `bioinf-librarian.md`

Expert agent for finding bioinformatics tools, documentation, and resources. Use when you need to:
- Research specific tools (Seurat, ArchR, MACS2)
- Find parameter documentation
- Locate GitHub repositories
- Compare tool versions
- Investigate best practices

**Example usage:**
```
"I need documentation for MACS2 peak calling with ATAC-seq data"
```

### RNA-seq Methods Writer
**File:** `rnaseq-methods-writer.md`

Automatically generates publication-ready Methods sections from RNA-seq analysis code. Use when you need to:
- Document completed RNA-seq analyses
- Generate reproducible methods descriptions
- Extract statistical models from code
- Create formal scientific writing from scripts

**Example usage:**
```
"Generate a methods section from my RNA-seq analysis in ./analysis"
```

## Using These Agents

1. Agents are automatically loaded by Claude Code from this directory
2. Reference agents in your prompts when needed
3. Agents have specialized tools and knowledge for their domains

## Creating Custom Agents

See the existing agent files for template structure. Key components:
- Frontmatter with name, description, model, color
- Agent identity and expertise
- Methodology and protocols
- Example usage patterns

## Contributing New Agents

Have an idea for a useful bioinformatics agent? Please:
1. Create agent definition following existing format
2. Test thoroughly
3. Document use cases
4. Submit a pull request

See [CONTRIBUTING.md](../CONTRIBUTING.md) for details.

## Agent File Location

These agent files are symlinked to `.claude/agents/` so Claude Code can automatically discover them. The canonical location is this `agents/` directory in the repository root.
