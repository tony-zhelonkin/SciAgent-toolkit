# SciAgent-toolkit v3 Architecture Design

> **Date:** 2026-03-07
> **Status:** Design proposal
> **Prerequisite reading:** `../research/RESEARCH_20260307.md`
> **Scope:** Complete architecture redesign for provider-agnostic, human-in-the-loop scientific workflows

---

## Table of Contents

1. [Design Philosophy](#1-design-philosophy)
2. [Revised Tier System](#2-revised-tier-system)
3. [Component Architecture](#3-component-architecture)
4. [Workflow System Design](#4-workflow-system-design)
5. [Verification Gate Pattern](#5-verification-gate-pattern)
6. [Knowledge Persistence Design](#6-knowledge-persistence-design)
7. [Provider Agnosticism](#7-provider-agnosticism)
8. [Migration Path](#8-migration-path)

---

## 1. Design Philosophy

### 1.1 Core Principles

**P1: Human-in-the-loop by design, not afterthought.**
The AI is a mentor and coding partner. It proposes, the human verifies. No step
proceeds without the human confirming the output makes sense. This is not a
safety guardrail bolted on -- it is the fundamental interaction model.

**P2: Progressive complexity -- basics before production.**
Every workflow starts small: subset the data, test the method, understand the
output, then scale up. The system enforces this progression structurally. You
cannot skip the "test on 5000 cells" step to go straight to 1M cells.

**P3: Reproducibility through structure, not rigidity.**
The framework provides scaffolding (phases, gates, checkpoints) but does not
prescribe specific tools or parameters. The scientist chooses; the framework
records and enforces verification.

**P4: Knowledge accumulates across sessions.**
Every session builds on what came before. Decisions, parameter choices, failed
approaches, and validated findings persist in a structured knowledge base.
Sessions do not start from zero.

**P5: Provider-agnostic at the workflow level.**
Workflows, skills, and knowledge are defined in plain Markdown and YAML.
Provider-specific machinery (MCP servers, agent definitions, output styles)
lives in a thin adapter layer. Switching from Claude Code to Gemini CLI changes
the adapter, not the workflow.

**P6: Small surface area, deep capability.**
Fewer concepts that compose well beat many concepts that overlap. The v3 tier
system has 6 tiers (up from 4) but each tier has a single, clear responsibility
with no overlap.

### 1.2 Anti-Patterns This Design Prevents

| Anti-Pattern | How v3 Prevents It |
|---|---|
| AI runs ahead without verification | Verification gates block phase transitions |
| Session drift (over-engineering, tangents) | Workflow phases constrain scope per session |
| "Works on my machine" non-reproducibility | Knowledge base records environment, parameters, decisions |
| Starting over each session | Checkpoints + knowledge base provide structured resume |
| Black-box scripts that "just work" | Progressive complexity requires understanding before scaling |
| Provider lock-in | Workflow/skill layer is pure Markdown; adapter layer is thin |

---

## 2. Revised Tier System

### 2.1 Current System (v2)

```
Tier 1: Roles          -- what to activate
Tier 2: Agents/Skills  -- what's available
Tier 3: MCP Profiles   -- what tools are connected
Tier 4: Addons         -- what persists across switches
```

**Gaps in v2:**
- No concept of workflow phases or progression
- No verification enforcement
- No knowledge persistence
- Output styles are disconnected from roles
- Skills are flat files with no supporting assets
- Addons only cover MCP servers, not skills or styles

### 2.2 Proposed System (v3)

```
Tier 0: Knowledge Base       -- what we know (persists across everything)
Tier 1: Workflows            -- what phase we're in (NEW)
Tier 2: Roles                -- what persona is active (evolved)
Tier 3: Skills               -- what domain knowledge is loaded (evolved)
Tier 4: Agents               -- what sub-agents are available (stable)
Tier 5: Profiles + Addons    -- what tools are connected (merged)
```

### 2.3 Why This Ordering

The tiers are ordered from most persistent (top) to most volatile (bottom):

```
+---------------------------------------------------+
|  Tier 0: KNOWLEDGE BASE                           |
|  Persists: across all sessions, all workflows     |
|  Changes: only by appending findings/decisions     |
+---------------------------------------------------+
|  Tier 1: WORKFLOWS                                |
|  Persists: across sessions within a project        |
|  Changes: when moving to a new analysis phase      |
+---------------------------------------------------+
|  Tier 2: ROLES                                    |
|  Persists: within a workflow phase                 |
|  Changes: when switching interaction mode           |
+---------------------------------------------------+
|  Tier 3: SKILLS                                   |
|  Persists: within a role                           |
|  Changes: when role switches or addons toggle       |
+---------------------------------------------------+
|  Tier 4: AGENTS                                   |
|  Persists: within a role                           |
|  Changes: when role switches                        |
+---------------------------------------------------+
|  Tier 5: PROFILES + ADDONS                        |
|  Persists: until explicitly switched               |
|  Changes: when role switches or manual override     |
+---------------------------------------------------+
```

### 2.4 Key Differences from v2

| Aspect | v2 | v3 |
|---|---|---|
| Top-level concept | Role (persona) | Workflow (progression) |
| Knowledge | Ad-hoc (MEMORY.md, handoff.md) | Structured knowledge base (Tier 0) |
| Phase tracking | None | Workflow phases with gates |
| Verification | Manual, inconsistent | Encoded in workflow, logged |
| Skills format | Single .md file | Directory: SKILL.md + references/ + scripts/ + assets/ |
| Output style | Manual, global | Role-aware, per-role default with override |
| Addons | MCP servers only | MCP servers + skills + output styles |
| Profiles + Addons | Separate tiers | Merged into single "tooling" tier |

---

## 3. Component Architecture

### 3.1 Directory Structure (Toolkit)

```
SciAgent-toolkit/
    workflows/                        # NEW: Phase-based workflow definitions
        scrna-atlas/
            workflow.yaml             #   Phase definitions, gates, prerequisites
            phases/
                01-explore.yaml       #   Per-phase: skills, role, verification
                02-integrate.yaml
                03-annotate.yaml
                04-analyse.yaml
            knowledge-schema.yaml     #   What knowledge this workflow produces
        multiome-grn/
            workflow.yaml
            phases/
                ...

    roles/                            # EVOLVED: Now include output_style
        base.yaml
        scrna-atlas.yaml
        mentor.yaml                   #   NEW: Pedagogical mode
        production.yaml               #   NEW: Efficient execution mode

    skills/                           # EVOLVED: Directory-based format
        anndata/
            SKILL.md                  #   Metadata + navigation + triggers
            references/               #   Detailed topic guides
                memory-efficiency.md
                concat-strategies.md
            scripts/                  #   Utility scripts
            assets/                   #   Templates, configs
        scanpy/
            SKILL.md
            references/
                ...
        scvi-basic.md                 #   BACKWARD COMPAT: flat files still work

    agents/                           # STABLE: Agent definitions
        docs-librarian.md
        bio-interpreter.md
        ...

    templates/
        mcp-profiles/                 # STABLE: MCP server configs
        mcp-addons/                   # STABLE: Persistent MCP configs
        skill-addons/                 # NEW: Persistent skill configs
        output-styles/                # NEW: Role-aware output style templates
            mentor.md
            production.md
            socratic.md
        vendor/                       # STABLE: Project templates
        knowledge/                    # NEW: Knowledge base templates
            project-knowledge.yaml.template
            session-log.yaml.template

    scripts/
        activate-role.sh              # EVOLVED: Also sets output style
        activate-workflow.sh          # NEW: Phase-based workflow activation
        switch-mcp-profile.sh         # STABLE
        manage-addon.sh               # EVOLVED: Skills + MCP + styles
        verify-gate.sh                # NEW: Verification gate checker
        knowledge-update.sh           # NEW: Append to knowledge base

    knowledge/                        # NEW: Toolkit-level knowledge
        patterns/                     #   Reusable analysis patterns
            pseudobulk-de.md
            batch-aware-integration.md
        pitfalls/                     #   Common mistakes and fixes
            cluster-level-annotation.md
            naive-cross-batch-de.md
```

### 3.2 Directory Structure (Project)

```
project/
    .claude/
        agents/           -> symlinks to toolkit agents/
        skills/           -> symlinks to toolkit skills/ (dirs or files)
        output-styles/    -> symlinks to toolkit templates/output-styles/
        active-addons.json
        settings.local.json

    .sciagent/                        # NEW: Project-level SciAgent state
        active-workflow.yaml          #   Current workflow + phase pointer
        knowledge/
            findings.yaml             #   Accumulated findings
            decisions.yaml            #   Parameter choices with rationale
            failed-approaches.yaml    #   What didn't work and why
            checkpoints.yaml          #   Data checkpoint registry
            session-log/              #   Per-session logs
                2026-03-07_001.yaml
                2026-03-07_002.yaml
        verification/
            gates-passed.yaml         #   Which verification gates passed
            evidence/                 #   Verification outputs (plots, stats)
```

### 3.3 Component Relationships

```
                    +------------------+
                    |  WORKFLOW        |
                    |  (scrna-atlas)   |
                    +--------+---------+
                             |
                    selects phase
                             |
                             v
            +----------------+----------------+
            |                |                |
     +------v------+  +-----v------+  +------v------+
     | PHASE:      |  | PHASE:     |  | PHASE:      |
     | 01-explore  |  | 02-integ   |  | 03-annotate |
     +------+------+  +-----+------+  +------+------+
            |                |                |
            | declares       | declares       | declares
            v                v                v
     +------+------+  +-----+------+  +------+------+
     | ROLE:       |  | ROLE:      |  | ROLE:       |
     | mentor      |  | mentor     |  | production  |
     +------+------+  +-----+------+  +------+------+
            |                |                |
            | activates      | activates      | activates
            v                v                v
     +------+------+  +-----+------+  +------+------+
     | Skills:     |  | Skills:    |  | Skills:     |
     | anndata     |  | scvi-basic |  | scvi-scanvi |
     | scanpy      |  | scanpy     |  | cellxgene.. |
     | sc-rna-qc   |  | anndata    |  | anndata     |
     +------+------+  +-----+------+  +------+------+
            |                |                |
            | loads          | loads          | loads
            v                v                v
     +------+------+  +-----+------+  +------+------+
     | Agents:     |  | Agents:    |  | Agents:     |
     | bio-interp  |  | bio-interp |  | bio-interp  |
     | insight-exp |  | insight-exp|  | code-review |
     +------+------+  +-----+------+  +------+------+
            |                |                |
            | connects       | connects       | connects
            v                v                v
     +------+------+  +-----+------+  +------+------+
     | Profile:    |  | Profile:   |  | Profile:    |
     | coding      |  | coding     |  | coding      |
     | + addons    |  | + addons   |  | + addons    |
     +-------------+  +------------+  +-------------+

            |                |                |
            | reads/writes   | reads/writes   | reads/writes
            v                v                v
     +----------------------------------------------+
     |              KNOWLEDGE BASE                   |
     |  findings | decisions | checkpoints | gates   |
     +----------------------------------------------+
```

### 3.4 Component Specifications

#### 3.4.1 Workflows (NEW)

**Purpose:** Define phase-based progressions for multi-session analysis projects.

**File:** `workflows/<name>/workflow.yaml`

```yaml
name: scrna-atlas
description: scRNA-seq atlas construction from raw counts to annotated object
version: 1.0.0

# Phases are sequential. Each phase has prerequisites (gates from prior phase).
phases:
  - id: explore
    file: phases/01-explore.yaml
    description: Explore each dataset independently
    prerequisites: []      # First phase has no prerequisites

  - id: integrate
    file: phases/02-integrate.yaml
    description: Batch integration with scVI
    prerequisites:
      - gate: explore/datasets-characterized
      - gate: explore/batch-structure-understood

  - id: annotate
    file: phases/03-annotate.yaml
    description: Cell type annotation
    prerequisites:
      - gate: integrate/integration-validated
      - gate: integrate/batch-mixing-confirmed

  - id: analyse
    file: phases/04-analyse.yaml
    description: Differential expression, pathway analysis, visualization
    prerequisites:
      - gate: annotate/labels-validated
      - gate: annotate/marker-consistency-confirmed

# What knowledge this workflow produces
knowledge_outputs:
  - type: finding
    description: Cell type composition across conditions
  - type: decision
    description: Integration parameters and batch key choices
  - type: checkpoint
    description: Annotated AnnData objects
```

#### 3.4.2 Phase Definitions (NEW)

**File:** `workflows/<name>/phases/<nn>-<phase>.yaml`

```yaml
id: explore
name: Dataset Exploration
description: |
  Explore each dataset independently before combining.
  Goal: understand what you have before deciding how to integrate.

# What role to use during this phase
default_role: mentor

# What skills this phase needs (merged with role skills)
phase_skills:
  - anndata
  - scanpy
  - single-cell-rna-qc

# Small-sample directive: test before full scale
small_sample:
  enabled: true
  strategy: random-subset
  size: 5000         # cells
  purpose: "Verify QC thresholds and basic structure before processing full dataset"

# Verification gates that must pass before moving to next phase
gates:
  - id: datasets-characterized
    description: "Each dataset explored: shape, genes, batches, basic QC metrics"
    evidence_type: summary     # summary | plot | statistic | checklist
    verification: |
      Confirm for each dataset:
      - [ ] Number of cells and genes reported
      - [ ] Batch structure documented (how many batches, cells per batch)
      - [ ] Basic QC metrics computed (genes/cell, counts/cell, mito%)
      - [ ] QC thresholds chosen with rationale

  - id: batch-structure-understood
    description: "Batch composition documented, integration strategy decided"
    evidence_type: checklist
    verification: |
      Confirm:
      - [ ] Batch key identified and documented
      - [ ] Confounding variables identified (tissue, donor, technology)
      - [ ] Decision recorded: which covariates to model in integration

# Suggested session structure
session_guide:
  - title: "Session 1: Load and inspect first dataset"
    goal: "Understand shape, metadata, available layers"
    checkpoint: "01_dataset1_raw.h5ad"

  - title: "Session 2: QC on first dataset"
    goal: "Compute QC metrics, choose thresholds, filter"
    checkpoint: "01_dataset1_qc.h5ad"

  - title: "Session 3: Repeat for remaining datasets"
    goal: "Apply consistent QC across all datasets"
    checkpoint: "01_all_datasets_qc/"

# Human decision points -- AI must stop and ask
decision_points:
  - trigger: "QC threshold selection"
    instruction: |
      Present the QC metric distributions (violin plots or histograms).
      Propose thresholds with rationale. DO NOT apply filters until the
      human confirms the thresholds look reasonable for their biology.

  - trigger: "Integration strategy"
    instruction: |
      Present the batch structure summary. Propose an integration approach
      (scVI with specific covariates). Explain trade-offs. Wait for
      human approval before proceeding.
```

#### 3.4.3 Roles (EVOLVED)

Roles gain two new fields: `output_style` and `verification_mode`.

**File:** `roles/<name>.yaml`

```yaml
name: mentor
description: >
  Pedagogical interaction mode. AI explains reasoning, proposes one step at a
  time, verifies understanding before proceeding. Default for exploration and
  learning phases.
mcp_profile: coding

# NEW: Output style to activate with this role
output_style: socratic

# NEW: How strictly to enforce verification
verification_mode: strict    # strict | normal | relaxed
# strict:  AI must ask for confirmation after every data transformation
# normal:  AI asks for confirmation at gate boundaries
# relaxed: AI proceeds but logs what it would have verified

agents:
  - docs-librarian
  - bio-interpreter
  - insight-explorer
  - captions
  - doc-curator
  - code-reviewer
  - handoff

skills:
  - anndata
  - scanpy
  - single-cell-rna-qc
  # Phase-specific skills are merged on top of role skills
```

```yaml
name: production
description: >
  Efficient execution mode. AI provides complete code blocks, focuses on
  correctness over explanation. For phases where the approach is already
  understood and the goal is execution.
mcp_profile: coding

output_style: production
verification_mode: normal

agents:
  - docs-librarian
  - bio-interpreter
  - code-reviewer
  - handoff

skills: []   # Populated entirely by workflow phase
```

#### 3.4.4 Skills (EVOLVED)

Skills support both the current flat-file format (backward compatible) and the
new directory format from claude-scientific-skills.

**Flat format (v2, still supported):**
```
skills/
    scvi-basic.md                    # Single file, works as before
```

**Directory format (v3, preferred for new/complex skills):**
```
skills/
    anndata/
        SKILL.md                     # Required: metadata, triggers, navigation
        references/                  # Optional: detailed guides
            memory-efficiency.md     #   When to use backed mode, chunked ops
            concat-strategies.md     #   outer vs inner join, batch handling
            subsetting-patterns.md   #   View vs copy, boolean vs fancy indexing
        scripts/                     # Optional: utility scripts
            validate_anndata.py      #   Check h5ad integrity
        assets/                      # Optional: templates, schemas
            anndata_checklist.md     #   Post-load verification checklist
```

**SKILL.md frontmatter (v3):**
```yaml
---
name: anndata
description: "AnnData object manipulation, I/O, subsetting, concatenation"
version: 1.0.0
tags: ["scverse", "data-structure", "foundation"]
complementary_skills:
  - scanpy
  - scvi-framework
contraindications:
  - "For Seurat objects, use anndatar-seurat-scanpy-conversion instead"
  - "For multimodal data, use multimodal-anndata-mudata instead"
references:
  - memory-efficiency.md
  - concat-strategies.md
  - subsetting-patterns.md
---

# AnnData: Annotated Data Containers

## When to Use
...

## When NOT to Use
...

## Core Concepts
...

## Quick Reference
...

## See Also
- `scanpy` skill for analysis operations on AnnData objects
- `scvi-framework` skill for model setup with AnnData
```

**activate-role.sh evolution:** The script detects whether a skill name
resolves to `skills/<name>.md` (flat) or `skills/<name>/` (directory).
For directories, it symlinks the entire directory. Claude Code reads
`.claude/skills/<name>/SKILL.md` as the entry point and can access
`references/`, `scripts/`, `assets/` when it needs deeper information.

#### 3.4.5 Agents (STABLE)

No structural changes. Agents remain `.md` files with YAML frontmatter in
`agents/`. The current 8 agents continue to work.

One addition: agents can optionally declare `workflow_affinity` to indicate
which workflow phases they are most useful in:

```yaml
---
name: insight-explorer
description: |
  Explore data with scientific skepticism...
model: sonnet
workflow_affinity:
  - explore
  - integrate
---
```

This is informational only -- it helps roles and workflows select appropriate
agent sets, but does not restrict agent usage.

#### 3.4.6 Profiles + Addons (MERGED)

Profiles and addons merge into a single "tooling" tier. The conceptual
distinction remains (profiles are base configurations, addons persist across
switches), but they are managed together.

**New addon types:**

| Addon Type | v2 | v3 |
|---|---|---|
| MCP servers | Yes | Yes |
| Skills | No | Yes (skill-addons) |
| Output styles | No | Yes (style-addons) |

**Skill addons** allow skills to persist across role switches (like MCP addons
persist across profile switches). Use case: the `anndata` skill is useful in
every phase. Rather than listing it in every role, make it a skill addon.

**File:** `templates/skill-addons/anndata.addon.yaml`
```yaml
_meta:
  name: anndata
  description: "AnnData skill persists across all role switches"
  type: skill
skill: anndata
```

**active-addons.json evolution:**
```json
{
  "addons": {
    "notebook-tools": { "enabled": true, "type": "mcp" },
    "anndata": { "enabled": true, "type": "skill" },
    "socratic": { "enabled": false, "type": "style" }
  }
}
```

#### 3.4.7 Output Styles (EVOLVED)

Output styles move from project-level manual config to toolkit-level templates
that roles can select.

**File:** `templates/output-styles/socratic.md`

This is the same format as the current `cs101_v0.1.md` -- a Markdown file with
YAML frontmatter and interaction instructions. The change is that roles
declare their default style, and `activate-role.sh` symlinks the appropriate
style into `.claude/output-styles/`.

**Role-style mapping:**
```
mentor role      -> socratic style (explain, verify, one-cell-at-a-time)
production role  -> production style (complete code, minimal explanation)
planning role    -> structured style (outlines, checklists, decision matrices)
```

Users can always override: `activate-role.sh mentor --output-style production`

#### 3.4.8 Knowledge Base (NEW)

The knowledge base is the most significant new component. It replaces ad-hoc
MEMORY.md entries, handoff documents, and implicit session context with a
structured, queryable record of what the project knows.

**Location:** `.sciagent/knowledge/` in the project directory.

**Components:**

```
.sciagent/knowledge/
    findings.yaml            # Validated scientific findings
    decisions.yaml           # Parameter and method choices with rationale
    failed-approaches.yaml   # What was tried and why it failed
    checkpoints.yaml         # Data checkpoint registry
    session-log/             # Per-session structured logs
        2026-03-07_001.yaml
```

**findings.yaml format:**
```yaml
findings:
  - id: F001
    date: 2026-03-07
    phase: explore
    gate: datasets-characterized
    summary: "PanCancer tumor dataset has 4.1M cells, 29 cancer types, 6 tissue categories"
    details: |
      - cancerType: 29 unique values
      - tissue: Tumor, Adjacent, Blood, Metastasis, PreLesion, PleuralFluids
      - ann1 has 605K Myeloid cells
      - majorCluster has 637K Myeloid cells
      - 21,812 genes
    evidence: "02_analysis/notebooks/L01_explore_pancancer.ipynb, cell 5"
    confidence: high

  - id: F002
    date: 2026-03-07
    phase: integrate
    gate: integration-validated
    summary: "scVI integration with batch_key='batch' produces good mixing across 115 batches"
    details: |
      - n_latent=30, n_hidden=128
      - batch_key: combined dataset+donor identifier
      - 1,051,354 cells after integration
    evidence: "03_results/checkpoints/myeloid_integrated_v2.h5ad"
    confidence: high
```

**decisions.yaml format:**
```yaml
decisions:
  - id: D001
    date: 2026-03-07
    phase: integrate
    question: "What batch key to use for scVI integration?"
    choice: "Combined dataset+donor identifier as 'batch'"
    rationale: |
      Each donor within each dataset is a biological replicate.
      Using dataset alone would under-correct; using sample_id
      might over-correct. The combined key balances batch
      correction strength.
    alternatives_considered:
      - "dataset only: too coarse, leaves donor effects"
      - "sample_id: too fine, risks removing biological signal"
    reversible: true
    impact: "All downstream integration results depend on this"
```

**failed-approaches.yaml format:**
```yaml
failed_approaches:
  - id: X001
    date: 2026-03-07
    phase: annotate
    what: "Cluster-level annotation via max(cluster_mean_score)"
    why_failed: |
      Assigns one label per Leiden cluster. Mixed clusters get wrong label.
      cDC1/cDC2 markers swapped in dotplot validation. DC3 became 31% catch-all.
    lesson: "Always annotate at cell level. Always validate with dotplot."
    reference: "02_analysis/scripts/.archive_v2/03c_annotate_dc.py"
```

**checkpoints.yaml format:**
```yaml
checkpoints:
  - id: CP001
    date: 2026-03-07
    phase: explore
    file: "03_results/checkpoints/dcverse_myeloid.h5ad"
    description: "DC-VERSE myeloid subset, 34,479 x 93,195"
    size_mb: 2100
    depends_on: []

  - id: CP002
    date: 2026-03-07
    phase: integrate
    file: "03_results/checkpoints/myeloid_integrated_v2.h5ad"
    description: "scVI-integrated myeloid cells, 1,051,354 x 21,812"
    size_mb: 49000
    depends_on: [CP001]
    parameters:
      n_latent: 30
      n_hidden: 128
      batch_key: batch
```

**Session log format:**
```yaml
session:
  id: "2026-03-07_001"
  start: "2026-03-07T14:00:00"
  workflow: scrna-atlas
  phase: explore
  role: mentor

  objectives:
    - "Explore PanCancer tumor metadata structure"
    - "Compute QC metrics for myeloid subset"

  outcomes:
    - type: finding
      ref: F001
    - type: checkpoint
      ref: CP001

  gates_attempted:
    - gate: explore/datasets-characterized
      status: partial
      remaining: ["QC thresholds not yet chosen"]

  next_steps:
    - "Choose QC thresholds for myeloid subset"
    - "Explore DC-VERSE healthy dataset with same approach"
```

### 3.5 Relationship to Existing MEMORY.md and Handoff Patterns

The knowledge base does NOT replace MEMORY.md. They serve different purposes:

| Component | Scope | Format | Managed by |
|---|---|---|---|
| `MEMORY.md` | User's personal notes, preferences, gotchas | Free-form Markdown | Claude Code auto-memory |
| `.sciagent/knowledge/` | Project analysis state, findings, decisions | Structured YAML | SciAgent scripts + AI |
| `handoff_*.md` | Session-to-session continuity | Free-form Markdown | Handoff agent |

**MEMORY.md** is personal ("I prefer scaffold-style code", "anndataR API
gotcha: pass layer names not 'X'"). The knowledge base is scientific ("scVI
integration used n_latent=30, batch_key='batch'").

Handoff documents become a lightweight summary pointing into the knowledge base
rather than containing all context themselves:

```markdown
# Session Handoff 2026-03-07

## Workflow: scrna-atlas, Phase: explore
## Knowledge updates: F001, D001, CP001
## Gates: explore/datasets-characterized PARTIAL
## Next: Choose QC thresholds (see decisions.yaml pending)
```

---

## 4. Workflow System Design

### 4.1 Lifecycle

```
activate-workflow.sh scrna-atlas --project-dir .
    |
    v
[1] Read workflows/scrna-atlas/workflow.yaml
    |
    v
[2] Create .sciagent/ directory structure
    |
    v
[3] Write .sciagent/active-workflow.yaml
    (points to phase: explore, the first phase)
    |
    v
[4] Activate the phase's default role
    (calls activate-role.sh mentor)
    |
    v
[5] Load phase-specific skills on top of role skills
    |
    v
[6] AI reads .sciagent/active-workflow.yaml at session start
    to know current phase, gates, and knowledge state
```

### 4.2 Phase Transitions

Phase transitions require all gates to pass:

```
Phase N (current)
    |
    | work happens, checkpoints saved, knowledge accumulated
    |
    v
Human says: "I think we're ready for the next phase"
    |
    v
AI runs: verify-gate.sh explore/datasets-characterized
    |
    +--> PASS: gate logged in .sciagent/verification/gates-passed.yaml
    |          evidence saved to .sciagent/verification/evidence/
    |
    +--> FAIL: AI reports what's missing
    |          suggests specific steps to satisfy the gate
    |          does NOT proceed
    |
    v (all gates pass)
Human confirms: "Advance to integrate phase"
    |
    v
activate-workflow.sh scrna-atlas --advance
    |
    v
[1] Verify all gates for current phase pass
[2] Update .sciagent/active-workflow.yaml (phase: integrate)
[3] Activate new phase's default role
[4] Load new phase's skills
```

### 4.3 Small-Sample Testing

Each phase can declare a `small_sample` directive. When present, the AI's
instructions include:

> Before running this on the full dataset, test on a {size}-cell random
> subset. Verify the output makes sense. Only then proceed to full scale.

This is enforced through the output style, not through code. The `socratic`
output style includes the instruction: "When a phase has small_sample
enabled, always propose a subset test first." The `production` style includes:
"When a phase has small_sample enabled, run subset test automatically and
report results before proceeding."

### 4.4 Decision Points

Decision points are moments where the AI MUST stop and present options to the
human. They are defined in the phase YAML and injected into the AI's context.

```
AI encounters a decision point
    |
    v
[1] Present the relevant data (plots, statistics, summaries)
[2] Propose options with trade-offs
[3] STOP. Wait for human input.
[4] Record the decision in decisions.yaml
[5] Proceed with the chosen approach
```

The AI does not skip decision points even in `production` mode. Decision points
are about scientific judgment, not about speed.

### 4.5 Session Continuity

At the start of every session, the AI reads:

1. `.sciagent/active-workflow.yaml` -- current phase
2. `.sciagent/knowledge/findings.yaml` -- what we know
3. `.sciagent/knowledge/decisions.yaml` -- what we decided
4. `.sciagent/knowledge/checkpoints.yaml` -- what data exists
5. Latest session log -- what happened last time

This replaces the "read the handoff doc" pattern with structured state. The AI
can construct a precise summary:

> We are in the **explore** phase of the **scrna-atlas** workflow.
> Last session we characterized the PanCancer tumor dataset (F001).
> The gate `datasets-characterized` is partially satisfied.
> Remaining: QC thresholds not yet chosen.
> Available checkpoints: dcverse_myeloid.h5ad (CP001).

### 4.6 Workflow Example: scrna-atlas

```
PHASE 1: EXPLORE
    Goal: Understand each dataset independently
    Role: mentor
    Skills: anndata, scanpy, single-cell-rna-qc
    Gates:
        [x] datasets-characterized
        [x] batch-structure-understood
    Knowledge produced: F001-F005, D001, CP001-CP003
        |
        | (all gates pass, human confirms)
        v
PHASE 2: INTEGRATE
    Goal: Batch integration with scVI
    Role: mentor
    Skills: anndata, scanpy, scvi-basic, scvi-framework
    Small sample: 5000 cells first
    Gates:
        [x] integration-validated
        [x] batch-mixing-confirmed
    Decision points:
        - batch_key selection
        - n_latent, n_hidden choices
        - covariate modeling strategy
    Knowledge produced: F006-F010, D002-D005, CP004-CP005
        |
        v
PHASE 3: ANNOTATE
    Goal: Cell type annotation
    Role: mentor -> production (once approach is validated)
    Skills: scvi-scanvi, cellxgene-census-annotation, scanpy
    Small sample: 10000 cells first
    Gates:
        [x] labels-validated
        [x] marker-consistency-confirmed
    Decision points:
        - Annotation strategy (scANVI vs manual vs hybrid)
        - Reference dataset selection
        - Marker panel validation
    Knowledge produced: F011-F015, D006-D008, CP006-CP007
        |
        v
PHASE 4: ANALYSE
    Goal: DE, pathway analysis, visualization
    Role: production
    Skills: scanpy, genenmf-metaprogram-discovery, scvi-mrvi
    Gates:
        [x] de-results-validated
        [x] figures-generated
    Knowledge produced: F016+, tables, figures
```

---

## 5. Verification Gate Pattern

### 5.1 Gate Anatomy

A verification gate has four components:

```yaml
gate:
  id: datasets-characterized           # Unique within workflow
  description: "Each dataset explored"  # Human-readable purpose
  evidence_type: checklist              # What constitutes proof
  verification: |                       # What to check
    - [ ] Number of cells and genes reported
    - [ ] Batch structure documented
    - [ ] Basic QC metrics computed
```

### 5.2 Evidence Types

| Type | What the AI produces | Storage |
|---|---|---|
| `summary` | Structured text summary of findings | `evidence/<gate_id>_summary.md` |
| `plot` | Specific plot(s) that demonstrate the finding | `evidence/<gate_id>_*.png` |
| `statistic` | Quantitative metric(s) with thresholds | `evidence/<gate_id>_stats.yaml` |
| `checklist` | All items checked with references | `evidence/<gate_id>_checklist.yaml` |

### 5.3 Gate Verification Flow

```
verify-gate.sh <gate-id>
    |
    v
[1] Read gate definition from workflow phase YAML
    |
    v
[2] Check .sciagent/verification/gates-passed.yaml
    for existing passing evidence
    |
    +--> Already passed? Report and exit.
    |
    v (not yet passed)
[3] For each checklist item or evidence requirement:
    - Check if finding/checkpoint exists in knowledge base
    - Check if evidence file exists
    |
    v
[4] Report status:
    PASS: All requirements met, log to gates-passed.yaml
    PARTIAL: Some requirements met, list remaining
    FAIL: No evidence found, list all requirements
```

### 5.4 How Verification Prevents AI From Skipping Steps

The verification system works at three levels:

**Level 1: Workflow YAML (structural).**
Phase transitions require all gates from the current phase to pass. The
`activate-workflow.sh --advance` command checks this programmatically.

**Level 2: Output style (behavioral).**
The `socratic` output style includes explicit instructions:

> After any data transformation, propose a verification step before
> proceeding. Do not write the next cell until the human has confirmed
> the current output looks correct.

**Level 3: Session context (contextual).**
The AI reads the active workflow state at session start. The phase definition
lists decision points and gates. The AI's system prompt includes:

> You are in phase: explore. This phase has 2 verification gates that
> must pass before advancing. Gate `datasets-characterized` is PARTIAL.
> Remaining items: QC thresholds not yet chosen.

This makes skipping verification unnatural -- the AI would have to actively
ignore its own context to skip a gate.

### 5.5 Verification Log Format

```yaml
# .sciagent/verification/gates-passed.yaml
gates:
  - id: explore/datasets-characterized
    passed_date: 2026-03-07
    passed_by: human     # human | automated
    evidence:
      - file: evidence/datasets-characterized_checklist.yaml
      - finding_refs: [F001, F002, F003]
    notes: "All 3 datasets characterized. QC thresholds chosen after violin plot review."

  - id: explore/batch-structure-understood
    passed_date: 2026-03-07
    passed_by: human
    evidence:
      - file: evidence/batch-structure_summary.md
      - finding_refs: [F004]
    notes: "115 batches total. batch_key = dataset + donor combination."
```

---

## 6. Knowledge Persistence Design

### 6.1 Two-Tier Knowledge

```
+----------------------------------+
|  TOOLKIT KNOWLEDGE               |
|  Location: SciAgent-toolkit/     |
|            knowledge/            |
|  Scope: Reusable across projects |
|  Content: Patterns, pitfalls,    |
|           best practices         |
|  Example: "Cluster-level         |
|   annotation fails when          |
|   clusters are mixed"            |
+----------------------------------+
               |
               | informs
               v
+----------------------------------+
|  PROJECT KNOWLEDGE               |
|  Location: .sciagent/knowledge/  |
|  Scope: This project only        |
|  Content: Findings, decisions,   |
|           checkpoints, sessions  |
|  Example: "scVI n_latent=30      |
|   chosen because 115 batches     |
|   need sufficient capacity"      |
+----------------------------------+
```

### 6.2 Toolkit Knowledge

Toolkit knowledge is curated, version-controlled content that represents
hard-won lessons applicable across projects.

**Patterns** (`knowledge/patterns/`):
```markdown
# Pseudobulk DE for Cross-Batch Comparisons

## Problem
Naive cell-level DE between conditions (e.g., tumor vs healthy) in integrated
scRNA-seq data detects batch effects, not biology. scVI corrects the latent
space but does not modify raw counts.

## Solution
Aggregate counts to pseudobulk (sum per donor per cell type), then use
standard bulk DE tools (DESeq2, edgeR) with batch as a covariate.

## When to Use
- Comparing conditions across batches
- Any cross-batch DE after scVI integration

## When NOT to Use
- Within-batch comparisons (cell-level DE is fine)
- Very few replicates per condition (pseudobulk needs n >= 3 per group)

## Code Template
...
```

**Pitfalls** (`knowledge/pitfalls/`):
```markdown
# Cluster-Level Annotation Failure Mode

## Symptom
Marker gene dotplots show swapped or inconsistent markers across annotated
cell types. Large "other" or "catch-all" categories.

## Root Cause
Assigning one label per Leiden cluster via max(cluster_mean_score). Mixed
clusters (containing multiple cell types) get the wrong label.

## Fix
Annotate at cell level using scANVI or cell-level marker scoring.
Always validate with dotplot before trusting annotations.

## Reference
DC_hum_verse v2 annotation failure (2026-03-06). See
.sciagent/knowledge/failed-approaches.yaml X001.
```

### 6.3 Knowledge Update Flow

```
Session work
    |
    v
AI proposes knowledge update:
  "I'd suggest recording this as a finding: ..."
    |
    v
Human confirms or edits
    |
    v
knowledge-update.sh --type finding --phase explore \
    --summary "PanCancer has 29 cancer types" \
    --evidence "notebook L01, cell 5"
    |
    v
Appends to .sciagent/knowledge/findings.yaml
Auto-increments ID (F001, F002, ...)
Timestamps automatically
```

### 6.4 Knowledge at Session Start

When a session begins, the AI's context includes a knowledge summary. For
large knowledge bases, this is a condensed view:

```
=== PROJECT KNOWLEDGE SUMMARY ===
Workflow: scrna-atlas | Phase: integrate | Gate status: 1/2 passed

Recent findings (last 5):
  F008: scVI converges in 200 epochs on 5K subset
  F009: batch mixing score = 0.87 (good) on 5K subset
  F010: UMAP shows expected separation of myeloid subtypes
  ...

Active decisions:
  D002: batch_key = "batch" (dataset+donor)
  D003: n_latent = 30, n_hidden = 128
  D004: covariate: tissue as categorical

Latest checkpoint: myeloid_integrated_v2.h5ad (CP005)

Failed approaches:
  X001: Cluster-level annotation (swapped markers)
  X002: n_latent=10 (insufficient for 115 batches)
```

This replaces reading MEMORY.md for project state. MEMORY.md continues to
store personal preferences and environment-specific notes.

---

## 7. Provider Agnosticism

### 7.1 Abstraction Layers

```
+--------------------------------------------------+
|  PROVIDER-AGNOSTIC LAYER                         |
|                                                  |
|  workflows/     - YAML phase definitions         |
|  skills/        - Markdown + references          |
|  knowledge/     - YAML structured state          |
|  roles/         - YAML declarations (minus       |
|                   provider-specific fields)       |
+--------------------------------------------------+
        |
        | thin adapter
        v
+--------------------------------------------------+
|  PROVIDER-SPECIFIC ADAPTER LAYER                 |
|                                                  |
|  Claude Code:                                    |
|    agents/*.md            (agent definitions)    |
|    .claude/               (settings, symlinks)   |
|    output-styles/*.md     (interaction modes)    |
|    .mcp.json              (MCP server config)    |
|    templates/mcp-profiles/  (profile templates)  |
|                                                  |
|  Gemini CLI:                                     |
|    .gemini/settings.json  (server config)        |
|    GEMINI.md              (context file)         |
|    templates/gemini-profiles/                    |
|                                                  |
|  Codex CLI:                                      |
|    ~/.codex/config.toml   (server config)        |
+--------------------------------------------------+
```

### 7.2 What Translates Across Providers

| Component | Claude Code | Gemini CLI | Codex CLI |
|---|---|---|---|
| Workflows | Read YAML, follow phases | Same | Same |
| Skills (SKILL.md) | `.claude/skills/` symlink | Included in GEMINI.md context | Included in prompt |
| Knowledge base | Read `.sciagent/` at session start | Same | Same |
| Roles (YAML) | Agents + skills + style | Skills only (no agent concept) | Skills only |
| Output styles | `.claude/output-styles/` | Prepended to GEMINI.md | Included in prompt |
| MCP servers | `.mcp.json` | `.gemini/settings.json` | `~/.codex/config.toml` |
| Agents | `.claude/agents/*.md` | No equivalent | No equivalent |

### 7.3 Provider-Specific Adapters

**activate-role.sh** already handles multi-provider output. In v3, it also
handles output styles and workflow-phase context:

```bash
activate-role.sh mentor --project-dir .
    |
    +--> Claude Code: symlink agents, skills, output-style to .claude/
    +--> Gemini CLI: update .gemini/settings.json, prepend style to GEMINI.md
    +--> Codex CLI: update ~/.codex/config.toml context
```

**activate-workflow.sh** is fully provider-agnostic. It writes to `.sciagent/`
which any provider can read. The role activation step calls the provider-
specific `activate-role.sh`.

### 7.4 Abstraction Boundary

The boundary sits between "what to do" and "how to tell the AI":

```
WHAT TO DO (agnostic)          HOW TO TELL THE AI (specific)
----------------------------   ----------------------------
workflow.yaml                  .claude/settings.local.json
phase definitions              .mcp.json
skill content (SKILL.md)       .claude/agents/*.md
knowledge state                output-styles/*.md
gate definitions               .gemini/settings.json
decision points                ~/.codex/config.toml
```

A workflow can be authored once and used with any provider. The adaptation
cost is limited to the thin adapter layer, which is already implemented in
`switch-mcp-profile.sh` for MCP configs and would extend to cover output
styles and workflow context injection.

---

## 8. Migration Path

### 8.1 Backward Compatibility Guarantees

| v2 Feature | v3 Status | Notes |
|---|---|---|
| `roles/*.yaml` | Works unchanged | New fields (output_style, verification_mode) are optional |
| `skills/*.md` (flat files) | Works unchanged | Directory format is additive, not replacing |
| `agents/*.md` | Works unchanged | `workflow_affinity` field is optional |
| `activate-role.sh` | Works unchanged | New flags are optional |
| `switch-mcp-profile.sh` | Works unchanged | No breaking changes |
| `templates/mcp-profiles/` | Works unchanged | New templates added, existing untouched |
| `templates/mcp-addons/` | Works unchanged | New addon types are additive |
| `active-addons.json` | Works unchanged | New `type` field defaults to `"mcp"` if absent |

### 8.2 Incremental Adoption Path

The migration is designed so each step is independently useful. Projects can
adopt any subset of v3 features without requiring the full stack.

**Step 0: Nothing changes.**
All v2 configurations continue to work. No action required.

**Step 1: Adopt directory-based skills.**
Convert high-value skills (anndata, scanpy, scvi-basic) to directory format.
Flat files continue to work alongside directories. This is the simplest win --
richer skill content with zero risk.

```
# Before (v2)
skills/anndata.md

# After (v3, both work simultaneously)
skills/anndata/SKILL.md
skills/anndata/references/memory-efficiency.md
skills/scvi-basic.md         # <-- still works as flat file
```

Update `activate-role.sh` to detect and symlink directories. This is a
backward-compatible change to the script.

**Step 2: Add output styles to roles.**
Create `templates/output-styles/` with 2-3 styles. Add `output_style` field
to role YAML files. Update `activate-role.sh` to handle styles.

Roles without `output_style` work as before (no style symlinked).

**Step 3: Add knowledge base to a project.**
Create `.sciagent/knowledge/` in one project. Start recording findings and
decisions manually. This requires no toolkit changes -- it is just a directory
convention.

Write `knowledge-update.sh` to make appending easier.

**Step 4: Add workflow definitions.**
Create the first workflow (e.g., `workflows/scrna-atlas/`). Write
`activate-workflow.sh`. This is the largest new component but it builds on
all previous steps.

Workflows are opt-in. Projects without workflows continue using roles directly.

**Step 5: Add verification gates.**
Add gate definitions to workflow phases. Write `verify-gate.sh`. This is
only useful after Step 4.

**Step 6: Add toolkit knowledge.**
Curate `knowledge/patterns/` and `knowledge/pitfalls/` from lessons learned
across projects. This is ongoing curation, not a one-time migration.

### 8.3 Implementation Priority

| Priority | Component | Effort | Impact |
|---|---|---|---|
| P0 | Directory-based skills | Low | High -- richer context for AI |
| P0 | Output styles in roles | Low | High -- consistent interaction modes |
| P1 | Knowledge base structure | Medium | High -- session continuity |
| P1 | knowledge-update.sh | Low | Medium -- easy knowledge recording |
| P2 | Workflow definitions | Medium | High -- phase-based progression |
| P2 | activate-workflow.sh | Medium | High -- workflow activation |
| P3 | Verification gates | Medium | Medium -- enforcement |
| P3 | verify-gate.sh | Low | Medium -- gate checking |
| P4 | Skill addons | Low | Low -- convenience |
| P4 | Toolkit knowledge curation | Ongoing | Medium -- cross-project learning |

### 8.4 File Changes Required

**Modified files:**

| File | Change |
|---|---|
| `scripts/activate-role.sh` | Detect skill dirs, handle output_style, read verification_mode |
| `roles/*.yaml` | Add optional output_style, verification_mode fields |
| `active-addons.json` schema | Add optional `type` field (defaults to "mcp") |

**New files:**

| File | Purpose |
|---|---|
| `scripts/activate-workflow.sh` | Workflow activation and phase advancement |
| `scripts/verify-gate.sh` | Verification gate checker |
| `scripts/knowledge-update.sh` | Append to knowledge base |
| `templates/output-styles/*.md` | Output style templates |
| `templates/skill-addons/*.addon.yaml` | Skill addon definitions |
| `templates/knowledge/*.template` | Knowledge base templates |
| `workflows/scrna-atlas/workflow.yaml` | First workflow definition |
| `workflows/scrna-atlas/phases/*.yaml` | Phase definitions |
| `knowledge/patterns/*.md` | Toolkit-level patterns |
| `knowledge/pitfalls/*.md` | Toolkit-level pitfall docs |

**No files deleted.** All v2 files remain functional.

---

## Appendix A: Complete Workflow Activation Sequence

```
User runs: activate-workflow.sh scrna-atlas --project-dir /workspaces/DC_hum_verse

[1] Parse workflows/scrna-atlas/workflow.yaml
[2] Determine current phase:
    - If .sciagent/active-workflow.yaml exists: read current phase
    - If not: create .sciagent/ structure, set phase to first phase (explore)
[3] Read phase definition: workflows/scrna-atlas/phases/01-explore.yaml
[4] Determine role: phase.default_role = "mentor"
[5] Call: activate-role.sh mentor --project-dir /workspaces/DC_hum_verse
    [5a] Symlink role's agents to .claude/agents/
    [5b] Symlink role's skills to .claude/skills/
    [5c] Symlink role's output style to .claude/output-styles/
    [5d] Switch MCP profile
[6] Merge phase_skills on top of role skills:
    - Phase declares: anndata, scanpy, single-cell-rna-qc
    - Role may already include some; duplicates are no-ops
    - Additional phase skills are symlinked
[7] Merge skill addons (from active-addons.json, type=skill)
[8] Write .sciagent/active-workflow.yaml:
    workflow: scrna-atlas
    phase: explore
    activated: 2026-03-07T14:00:00
    role: mentor
    profile: coding
[9] Print summary:
    "Workflow: scrna-atlas
     Phase: explore (1 of 4)
     Role: mentor (socratic style)
     Skills: 6 active (3 from role, 3 from phase, 2 overlap)
     Agents: 7 active
     Gates to pass: datasets-characterized, batch-structure-understood
     Knowledge: 0 findings, 0 decisions, 0 checkpoints"
```

## Appendix B: Comparison With Existing Patterns

| Pattern | Source | v3 Equivalent |
|---|---|---|
| Waypoint pattern | life-sciences (clinical-trial-protocol) | Knowledge base checkpoints + session logs |
| Validation gates | life-sciences (nextflow-development) | Verification gate system |
| Decision tree | life-sciences (scvi-tools) | Decision points in phase definitions |
| Progressive complexity | life-sciences (single-cell-rna-qc) | Small-sample directives + phase progression |
| Directory-based skills | claude-scientific-skills | SKILL.md + references/ + scripts/ + assets/ |
| Complementary skills | claude-scientific-skills | `complementary_skills` in SKILL.md frontmatter |
| Separation of concerns | RNAseq-toolkit | Phase = scope boundary, role = interaction mode boundary |
| Checkpoint caching | RNAseq-toolkit saveRDS() | checkpoints.yaml registry with dependency tracking |

## Appendix C: Glossary

| Term | Definition |
|---|---|
| **Workflow** | A named, phase-based progression for a complete analysis project |
| **Phase** | A stage within a workflow with defined scope, skills, and gates |
| **Gate** | A verification requirement that must pass before phase transition |
| **Decision point** | A moment where the AI must stop and present options to the human |
| **Knowledge base** | Structured YAML files recording findings, decisions, and checkpoints |
| **Finding** | A validated observation about the data or analysis |
| **Decision** | A parameter or method choice with recorded rationale |
| **Checkpoint** | A saved data artifact with provenance metadata |
| **Failed approach** | A documented attempt that did not work, with lessons |
| **Role** | A persona configuration: agents, skills, output style, verification mode |
| **Skill** | Domain knowledge in Markdown format, optionally with references and scripts |
| **Agent** | A sub-agent definition with specific capabilities and model assignment |
| **Profile** | An MCP server configuration template |
| **Addon** | A component (MCP server, skill, or style) that persists across role switches |
| **Output style** | An instruction document that shapes AI interaction behavior |
| **Provider** | An AI CLI tool (Claude Code, Gemini CLI, Codex CLI) |
| **Adapter** | Provider-specific configuration mapping |

---

*End of architecture design document.*
