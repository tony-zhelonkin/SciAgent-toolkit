# Claude Agents for Bioinformatics

Pre-configured Claude agents optimized for bioinformatics and scientific research workflows.

## Role-Based Activation

Agents are activated via the **role system**. When you activate a role, the specified agents are symlinked to `.claude/agents/`:

```bash
# Activate the base role (includes all 8 agents)
./scripts/activate-role.sh base --project-dir /path/to/project

# View available roles
ls roles/*.yaml
```

### Base Role (`roles/base.yaml`)

The base role includes all 8 agents organized by function:

```yaml
name: base
description: Default bioinformatics analysis role with full agent suite
mcp_profile: coding

agents:
  # Research & Documentation
  - bioinf-librarian            # Tool/documentation research
  - bio-research-visualizer     # Biological mechanism research + visualization

  # Data Exploration & Analysis
  - rnaseq-insight-explorer     # RNAseq data exploration with skepticism

  # Publication & Documentation
  - rnaseq-methods-writer       # Methods section generation from code
  - figure-caption-generator    # Publication-quality figure captions
  - repo-doc-curator            # Repository documentation cleanup

  # Code Review & Quality
  - refactor-stage-reviewer     # Code refactoring peer review

  # Session Management
  - handoff                     # Session handoff documentation

skills: []
```

## Available Agents (8 total)

### Research & Documentation

#### Bioinformatics Research Librarian
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

#### Bio-Research Visualizer
**File:** `bio-research-visualizer.md`

Deep web research into biological mechanisms underlying bioinformatics findings, followed by visualization recommendations. Use when you need to:
- Explain biological mechanisms behind gene expression patterns
- Research pathway interactions and crosstalk
- Identify cell types from marker expression
- Get visualization recommendations for biological insights

**Example usage:**
```
"I found upregulation of ISG15, MX1, OAS1. What's the biological mechanism and how should I visualize this?"
```

### Data Exploration & Analysis

#### RNAseq Insight Explorer
**File:** `rnaseq-insight-explorer.md`

Explore and interpret RNAseq analysis results with scientific skepticism. Use when you need to:
- Explore DESeq2/edgeR results stored in RDS/CSV files
- Validate user interpretations against actual data
- Identify unexpected patterns or overlooked findings
- Get visualization recommendations based on data characteristics

**Key Feature:** Applies critical thinking to avoid confirmation bias - verifies claims against actual data.

**Example usage:**
```
"I have DESeq2 results in results.rds, can you help me understand what's interesting here?"
```

### Publication & Documentation

#### RNA-seq Methods Writer
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

#### Figure Caption Generator
**File:** `figure-caption-generator.md`

Generates publication-quality scientific captions for bioinformatics outputs. **Fire-and-forget design** - writes directly to README.md files.

Use when you need to:
- Document figures, tables, and data files
- Trace outputs back to generating scripts
- Create Cell/Nature/Science quality figure legends
- Document output directories with comprehensive READMEs

**Usage Pattern (IMPORTANT):**
```
# Launch with run_in_background: true
# DO NOT call TaskOutput - wastes context
# Verify completion: ls <target_dir>/README.md
```

**Example usage:**
```
"Document the plots in 03_results/plots/QC/"
```

#### Repository Documentation Curator
**File:** `repo-doc-curator.md`

Audit, consolidate, and improve repository documentation. Use when you need to:
- Clean up messy documentation after refactoring
- Validate documentation accuracy against code
- Consolidate scattered docs into unified structure
- Prepare repository for external sharing

**Example usage:**
```
"The repo has gotten messy with outdated docs. Can you help organize it?"
```

### Code Review & Quality

#### Refactor Stage Reviewer
**File:** `refactor-stage-reviewer.md`

Independent peer review of refactored analysis stages. Use when you need to:
- Compare refactored code against original
- Verify compliance with project plan
- Assess reproducibility
- Review scientific integrity of code changes

**Example usage:**
```
"I've finished refactoring the preprocessing script. Can you review it against the original?"
```

### Session Management

#### Handoff
**File:** `handoff.md`

Create timestamped handoff documentation for session continuity. Use when you need to:
- Document completed work before ending a session
- Create checkpoint documentation
- Enable seamless session resumption

**Features:**
- Timestamped filenames (`handoff_YYYYMMDD_HHMMSS.md`)
- Auto-archives previous handoffs to `.handoff_archive/`
- Focused on immediate next steps

**Example usage:**
```
"Let's wrap up for today. We got the integration working."
```

## Using These Agents

1. **Activate a role** to symlink agents to `.claude/agents/`
2. Agents are automatically discovered by Claude Code
3. Reference agents in your prompts when needed
4. Agents have specialized tools and knowledge for their domains

## Creating Custom Agents

See the existing agent files for template structure. Key components:
- Frontmatter with name, description, model, color
- Agent identity and expertise
- Methodology and protocols
- Example usage patterns

### Adding to a Role

After creating an agent, add it to a role definition:

```yaml
# roles/my-custom-role.yaml
name: my-custom-role
description: Custom workflow role
agents:
  - bioinf-librarian
  - my-new-agent  # Your new agent
skills: []
```

## Directory Structure

```
SciAgent-toolkit/
├── agents/                           # Canonical agent definitions (8 total)
│   ├── bioinf-librarian.md           # Tool/documentation research
│   ├── bio-research-visualizer.md    # Biological mechanism research
│   ├── rnaseq-insight-explorer.md    # RNAseq data exploration
│   ├── rnaseq-methods-writer.md      # Methods section generation
│   ├── figure-caption-generator.md   # Publication figure captions
│   ├── repo-doc-curator.md           # Repository documentation
│   ├── refactor-stage-reviewer.md    # Code refactoring review
│   └── handoff.md                    # Session handoff
├── skills/                           # Custom skills (empty by default)
├── roles/
│   └── base.yaml                     # Role that references all agents
└── scripts/
    └── activate-role.sh              # Symlinks agents to .claude/agents/
```

When a role is activated, symlinks are created:
```
project/
└── .claude/
    └── agents/
        ├── bioinf-librarian.md -> /path/to/SciAgent-toolkit/agents/bioinf-librarian.md
        ├── bio-research-visualizer.md -> ...
        ├── rnaseq-insight-explorer.md -> ...
        ├── rnaseq-methods-writer.md -> ...
        ├── figure-caption-generator.md -> ...
        ├── repo-doc-curator.md -> ...
        ├── refactor-stage-reviewer.md -> ...
        └── handoff.md -> ...
```

## Contributing New Agents

Have an idea for a useful bioinformatics agent? Please:
1. Create agent definition following existing format
2. Add to appropriate role(s) in `roles/`
3. Test thoroughly
4. Document use cases
5. Submit a pull request

See [CONTRIBUTING.md](../CONTRIBUTING.md) for details.
