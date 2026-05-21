# SciAgent Toolkit

A modular framework for building AI-powered scientific research agents. Provides a role-based system for configuring agents and skills per project, with a curated library of bioinformatics agents and single-cell analysis skills.

## Features

- **Role-Based Configuration** via declarative YAML role definitions
- **Agent Library** — pre-configured Claude agents for bioinformatics workflows
- **Skills Library** — 50+ single-cell and multi-omics analysis skills
- **Template System** for consistent AI context files across projects

## Quick Start

```bash
# Clone the repository
git clone https://github.com/tony-zhelonkin/SciAgent-toolkit
cd SciAgent-toolkit

# Setup project (installs templates, activates base role)
./scripts/setup-ai.sh

# Or just activate a role directly
./scripts/activate-role.sh base --project-dir /path/to/project
```

## Role System & Custom Agents

```bash
# Activate a role (symlinks agents/skills into .claude/)
./scripts/activate-role.sh base --project-dir /path/to/project

# List available roles
ls roles/*.yaml
```

> **Scope:** one `.claude/` per invocation. No cascade into nested repos — re-run with `--project-dir` for each target. See `docs/workflows/architect/01-architecture.md` for details.

### Pre-configured Agents (Base Role)

**Research & Documentation**
- **Bioinformatics Research Librarian** — Find tools, docs, and resources via web research
- **Bio-Research Visualizer** — Deep biological mechanism research + visualization recommendations

**Data Exploration & Analysis**
- **RNA-seq Insight Explorer** — Explore RNAseq results with scientific skepticism

**Publication & Documentation**
- **RNA-seq Methods Writer** — Auto-generate publication methods sections from code
- **Figure Caption Generator** — Publication-quality captions for figures and tables (fire-and-forget)
- **Repo Doc Curator** — Audit and consolidate repository documentation

**Code Review & Quality**
- **Refactor Stage Reviewer** — Peer review of refactored analysis stages

**Session Management**
- **Handoff** — Timestamped session handoff documentation

### Creating Custom Roles

Create `roles/my-role.yaml`:
```yaml
name: my-role
description: Custom workflow role
agents:
  - bioinf-librarian
  - my-custom-agent
skills: []
```

See [agents/README.md](agents/README.md) for details.

## Docker/Container Deployment

The Docker test images in `docker/test/` are for CI/CD validation only.

For production containerized deployments, use [scbio-docker](https://github.com/tony-zhelonkin/scbio-docker).

## License

MIT License — see [LICENSE](LICENSE) for details.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.
