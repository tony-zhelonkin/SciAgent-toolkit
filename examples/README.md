# SciAgent Toolkit Examples

This directory contains practical examples and workflows for using the SciAgent Toolkit.

## Available Examples

- [Drug Discovery Workflow](drug-discovery-workflow.md) - Find and evaluate potential drug candidates
- [Genomics Research Workflow](genomics-research-workflow.md) - Analyze proteins, genes, and pathways
- [Literature Review Workflow](literature-review-workflow.md) - Conduct comprehensive literature searches

## How to Use These Examples

Each example includes:
- **Scientific context** and use case description
- **Prerequisites** needed (tools, knowledge, API keys)
- **Step-by-step instructions** with specific queries
- **Expected outputs** to validate your results
- **Common variations** for different scenarios
- **Troubleshooting tips** for common issues

### Getting Started

1. Ensure you have completed the [installation](../docs/INSTALLATION.md)
2. Verify MCP servers are loaded: `/mcp` in Claude Code or Codex CLI
3. Choose an example workflow that matches your research needs
4. Follow the step-by-step instructions
5. Adapt the queries for your specific research question

### Tips for Success

- **Start with the examples as written** to understand the workflow
- **Modify queries** to fit your specific research needs
- **Combine workflows** for complex research questions
- **Save useful queries** for reuse in your own research
- **Check API limits** for external services (PubMed, ChEMBL, etc.)

## Contributing Examples

Have a useful workflow to share? We welcome contributions!

See [CONTRIBUTING.md](../CONTRIBUTING.md) for guidelines on adding new examples.

Good examples to contribute:
- Domain-specific workflows (oncology, neuroscience, etc.)
- Multi-tool integration patterns
- Advanced analysis pipelines
- Real-world research case studies

## Example Template

When creating new examples, follow this structure:

```markdown
# Workflow Title

## Overview
Brief description of what this workflow accomplishes

## Prerequisites
- Required MCP servers
- External accounts/API keys needed
- Background knowledge assumed

## Use Case
Specific research scenario this addresses

## Steps

### Step 1: [Action Name]
**Query:**
```
"Your query here"
```

**Expected Output:**
Description of what you should see

**Notes:**
- Important considerations
- Variations
- Common issues

### Step 2: [Action Name]
...

## Variations

How to adapt this workflow for:
- Different research questions
- Alternative tools
- Specific domains

## Troubleshooting

Common issues and solutions

## Further Reading

Related resources and documentation
```

## Additional Resources

- [Installation Guide](../docs/INSTALLATION.md) - Setup instructions
- [Configuration Guide](../docs/CONFIGURATION.md) - Advanced configuration
- [FAQ](../docs/FAQ.md) - Frequently asked questions
- [Troubleshooting](../docs/TROUBLESHOOTING.md) - Common issues
