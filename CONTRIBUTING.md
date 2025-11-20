# Contributing to SciAgent Toolkit

Thank you for your interest in contributing! This toolkit thrives on community contributions.

## How to Contribute

### Reporting Issues
- Use GitHub Issues to report bugs
- Include reproduction steps, environment details, and error messages
- Check existing issues first to avoid duplicates

### Suggesting Enhancements
- Open an issue with the "enhancement" label
- Describe the use case and expected behavior
- Consider creating a proof-of-concept

### Contributing Code

#### Getting Started
1. Fork the repository
2. Clone your fork: `git clone [your-fork-url]`
3. Create a branch: `git checkout -b feature/your-feature-name`

#### Development Guidelines
- Follow existing code style (bash scripts with proper error handling)
- Add comments for complex logic
- Test scripts on multiple platforms (macOS, Linux)
- Update documentation for new features

#### Adding New MCP Servers
1. Create setup script in `scripts/mcp_servers/setup_newserver.sh`
2. Follow the existing script structure (error handling, logging, verification)
3. Add to main orchestrator in `scripts/setup_mcp_infrastructure.sh`
4. Document in `docs/INSTALLATION.md`

Example script structure:
```bash
#!/usr/bin/env bash
#
# Script Name: setup_newserver.sh
# Description: Install and configure NewServer MCP
# Author: Your Name
# Version: 1.0
#

set -e  # Exit on error

echo "Installing NewServer MCP..."

# Installation logic here
# - Check prerequisites
# - Install dependencies
# - Configure server
# - Test installation

echo "NewServer MCP installed successfully"
```

#### Adding New Agents
1. Create agent definition in `agents/your-agent-name.md`
2. Follow the existing agent format
3. Document in `agents/README.md`

Agent file structure:
```markdown
---
name: "Your Agent Name"
description: "Brief description"
model: "sonnet"
color: "blue"
---

# Agent Identity
Your agent's role and capabilities...

# Methodology
How the agent works...

# Examples
Example use cases...
```

#### Pull Request Process
1. Update relevant documentation
2. Test your changes thoroughly
3. Commit with clear, descriptive messages
4. Push to your fork
5. Open a Pull Request with:
   - Clear description of changes
   - Related issue numbers
   - Testing performed
   - Screenshots (if UI changes)

### Documentation Contributions
- Fix typos, improve clarity, add examples
- Documentation is as valuable as code!
- Update docs when adding features

### Community Guidelines
- Be respectful and inclusive
- Follow the [Code of Conduct](CODE_OF_CONDUCT.md)
- Help others in issues and discussions

## Development Setup

```bash
# Clone and test
git clone [repository-url]
cd SciAgent-toolkit
./scripts/setup_mcp_infrastructure.sh --help
```

## Testing Your Changes

### Test Installation Scripts

```bash
# Test individual components
./scripts/install_claude.sh
./scripts/mcp_servers/setup_serena.sh

# Test full installation
./scripts/setup_mcp_infrastructure.sh
```

### Test Configuration Changes

```bash
# Validate JSON
python3 -m json.tool .mcp.json

# Test in Claude Code
claude
/mcp

# Run diagnostics
claude doctor
```

### Test Documentation

- Check markdown rendering
- Verify all links work
- Test example commands
- Ensure consistency across docs

## Contributing Examples

Example workflows in `examples/` are highly valuable! When adding examples:

1. Create a markdown file in `examples/`
2. Include:
   - Use case description
   - Step-by-step instructions
   - Expected outputs
   - Common variations
3. Test the workflow thoroughly
4. Update `examples/README.md`

Example structure:
```markdown
# Workflow Title

## Overview
What this workflow accomplishes...

## Prerequisites
- Required tools
- API keys needed
- Knowledge assumed

## Steps

### Step 1: [Action]
```
Query or command here
```

Expected output:
```
Output example
```

### Step 2: [Action]
...

## Common Variations
- Alternative approaches
- Different use cases

## Troubleshooting
- Common issues
- Solutions
```

## Commit Message Guidelines

Use clear, descriptive commit messages:

```bash
# Good examples
git commit -m "Add support for custom ToolUniverse filters"
git commit -m "Fix PATH configuration in install_claude.sh"
git commit -m "Update INSTALLATION.md with macOS M1 instructions"

# Bad examples
git commit -m "Fix bug"
git commit -m "Update docs"
git commit -m "Changes"
```

Format:
- Use imperative mood ("Add feature" not "Added feature")
- Keep first line under 72 characters
- Add detailed description if needed

## Code Style

### Bash Scripts

```bash
#!/usr/bin/env bash
#
# Always include descriptive header
#

set -e  # Exit on error

# Use functions for reusability
function install_dependency() {
    local dep_name=$1
    echo "Installing ${dep_name}..."
    # Installation logic
}

# Clear error messages
if ! command -v python3 &> /dev/null; then
    echo "Error: Python 3 is required but not installed"
    exit 1
fi

# Proper logging
echo "Step 1: Checking prerequisites..."
echo "✓ Prerequisites met"
```

### Documentation

- Use clear headings
- Include code examples
- Provide context and rationale
- Link to related docs
- Test all commands

## Review Process

1. **Automated checks**: All PRs run automated validation
2. **Code review**: Maintainers review changes
3. **Testing**: Changes are tested on multiple platforms
4. **Documentation**: Ensure docs are updated
5. **Merge**: Approved PRs are merged

## Questions?

- Open a discussion on GitHub
- Check the [FAQ](docs/FAQ.md)
- Review [existing issues](link-to-issues)

## Recognition

Contributors are recognized in:
- GitHub contributors list
- Release notes
- Project documentation (for significant contributions)

Thank you for contributing to SciAgent Toolkit!
