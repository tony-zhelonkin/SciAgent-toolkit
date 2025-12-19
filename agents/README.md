# Claude Agents for Bioinformatics

Pre-configured Claude agents for bioinformatics and scientific research workflows.

## Quick Reference

| Agent | Purpose | Trigger |
|-------|---------|---------|
| `docs-librarian` | Find tool/package documentation | "What are the parameters for...?" |
| `bio-interpreter` | Research biological mechanisms | "What's the mechanism behind...?" |
| `insight-explorer` | Explore data files with skepticism | "What's interesting in this data?" |
| `captions` | Generate figure legends | "Document the plots in..." |
| `doc-curator` | Clean up repo documentation | "Help organize the docs" |
| `code-reviewer` | Review refactored code | "Review my refactoring" |
| `handoff` | Create session handoff docs | "Let's wrap up for today" |

## Agent Separation Logic

```
┌─────────────────────────────────────────────────────────────────┐
│                 docs-librarian                                   │
│  Find DOCUMENTATION for tools/packages                          │
│  Input: Tool/package name                                        │
│  Output: web_notes.md with docs, parameters, best practices     │
│  Tools: Context7, WebSearch, WebFetch                           │
└─────────────────────────────────────────────────────────────────┘
                              ↓ user has findings
┌─────────────────────────────────────────────────────────────────┐
│                    WORKFLOW HANDOFF                              │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  insight-explorer                    bio-interpreter             │
│  ─────────────────                   ──────────────             │
│  INPUT: Data files (RDS/CSV)         INPUT: Gene/pathway names  │
│  ACTION: Statistical exploration     ACTION: Literature research │
│  OUTPUT: Data patterns + viz recs    OUTPUT: web_notes.md       │
│                                                                  │
│  "What patterns are in this data?"   "What does this mean       │
│                                       biologically?"             │
│                                                                  │
│         ↓ discovers patterns                                     │
│         ↓ identifies genes/pathways                             │
│         └──────────────────────────→ researches mechanisms      │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

## Role-Based Activation

Agents are activated via the role system:

```bash
# Activate base role (includes all 7 agents)
./scripts/activate-role.sh base --project-dir /path/to/project

# View available roles
ls roles/*.yaml
```

### Base Role (`roles/base.yaml`)

```yaml
name: base
description: Default bioinformatics analysis role
mcp_profile: hybrid-research

agents:
  # Research & Documentation
  - docs-librarian      # Tool/package documentation
  - bio-interpreter     # Biological mechanism research
  # Data Exploration
  - insight-explorer    # Data exploration with skepticism
  # Publication
  - captions            # Figure caption generation
  - doc-curator         # Repository documentation
  # Code Quality
  - code-reviewer       # Refactoring peer review
  # Session Management
  - handoff             # Session handoff documentation

skills:
  - annotate-rnaseq-data
```

## Agent Details

### docs-librarian
**Purpose:** Find documentation for bioinformatics tools, libraries, packages, workflows.

**Tools:** Context7 (for packages), WebSearch, WebFetch

**Use when:**
- Need current parameters for a tool
- Looking for GitHub issues/solutions
- Researching best practices
- Comparing tool versions

**Output:** Findings recorded in `web_notes.md`

---

### bio-interpreter
**Purpose:** Research biological mechanisms underlying findings via literature.

**Use when:**
- Gene expression patterns need biological explanation
- Pathway enrichment needs mechanistic context
- Cell types need functional characterization

**Output:** Mechanism explanation in `web_notes.md`

---

### insight-explorer
**Purpose:** Explore data files (RDS, CSV) with scientific skepticism.

**Use when:**
- Have results files that need exploration
- User makes claims that should be validated
- Need visualization recommendations based on data

**Output:** Structured report with data patterns and viz recommendations

---

### captions
**Purpose:** Generate publication-quality figure captions. Fire-and-forget design.

**Use when:**
- Need to document plots/figures
- Creating README.md for output directories

**Output:** Writes directly to README.md (don't call TaskOutput)

---

### doc-curator
**Purpose:** Audit and consolidate repository documentation.

**Use when:**
- Documentation has become messy/scattered
- Preparing repo for sharing
- Need to validate docs against code

**Output:** Cleaned, consolidated documentation

---

### code-reviewer
**Purpose:** Rigorous peer review of refactored code.

**Use when:**
- Completed a refactoring
- Need to verify functionality preserved
- Want reproducibility assessment

**Output:** Structured review report with verdict

---

### handoff
**Purpose:** Create timestamped session handoff documentation.

**Use when:**
- Completing a work session
- Reaching a checkpoint
- Need to document progress for continuity

**Output:** `handoff_YYYYMMDD_HHMMSS.md` with archive management

## Directory Structure

```
agents/
├── docs-librarian.md     # Tool/package documentation
├── bio-interpreter.md    # Biological mechanism research
├── insight-explorer.md   # Data exploration
├── captions.md           # Figure caption generation
├── doc-curator.md        # Repository documentation
├── code-reviewer.md      # Refactoring review
├── handoff.md            # Session handoff
├── README.md             # This file
└── .deprecated/
    └── rnaseq-methods-writer.md  # Deprecated
```

## Creating Custom Agents

1. Create `agents/new-agent.md` with YAML frontmatter:
   ```yaml
   ---
   name: new-agent
   description: |
     When to use this agent...
     <example>...</example>
   tools: Tool1, Tool2, ...
   model: sonnet
   color: blue
   ---

   Agent instructions here...
   ```

2. Add to role in `roles/base.yaml`:
   ```yaml
   agents:
     - new-agent
   ```

3. Activate role:
   ```bash
   ./scripts/activate-role.sh base --project-dir /path/to/project
   ```

## Contributing

See [CONTRIBUTING.md](../CONTRIBUTING.md) for guidelines.
