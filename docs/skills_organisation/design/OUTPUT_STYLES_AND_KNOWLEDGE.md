# Design: Output Styles and Knowledge Persistence Systems

> **Date:** 2026-03-07
> **Status:** Design Document (v0.1)
> **Scope:** SciAgent-toolkit v3 — interaction modes and cross-session knowledge
> **Dependencies:** Research from `RESEARCH_20260307.md`, existing `cs101_v0.1.md` output style

---

## Table of Contents

**Part A: Output Styles System**
1. [Style Definitions](#1-style-definitions)
2. [Proposed Styles](#2-proposed-styles)
3. [Style Switching](#3-style-switching)
4. [Provider Agnosticism](#4-provider-agnosticism)

**Part B: Knowledge Persistence System**
5. [Knowledge Architecture](#5-knowledge-architecture)
6. [Decision Log Format](#6-decision-log-format)
7. [Verification History](#7-verification-history)
8. [Session Continuity](#8-session-continuity)
9. [Auto-Memory vs Structured Knowledge](#9-auto-memory-vs-structured-knowledge)
10. [Integration with Workflow System](#10-integration-with-workflow-system)

---

# Part A: Output Styles System

## 1. Style Definitions

### 1.1 What an Output Style Is

An output style defines HOW the AI interacts with the user — tone, pacing, degree of autonomy, verification expectations, and code delivery format. It does NOT define WHAT the AI knows (that is the role of skills and agents).

Separation of concerns:
- **Roles** declare what agents, skills, and MCP tools are available
- **Skills** provide domain knowledge and tool-specific reference
- **Output styles** control interaction behavior within a session

### 1.2 File Format

Output styles are Markdown files with YAML frontmatter, stored in `.claude/output-styles/` (Claude Code) or equivalent locations for other providers.

**Required structure:**

```markdown
---
name: style-identifier          # kebab-case
description: |
  One-paragraph summary of when and why to use this style.
  Acts as a trigger for automatic selection.
phase_affinity:                 # Which workflow phases this style suits
  - exploration
  - learning
autonomy_level: low             # low | medium | high
---

# [Style Name] Instructions

## Core Philosophy
What this style optimizes for (3-5 bullet points).

## Interaction Rules
How the AI behaves — pacing, question-asking, explanation depth.

## Code Delivery
How code is presented — scaffold vs complete, one cell vs multi-cell,
verification expectations.

## Verification Expectations
What the AI should check and surface. When to stop and ask.

## Boundaries
What this style does NOT do. When to suggest switching styles.
```

### 1.3 Section Definitions

| Section | Purpose | Example Content |
|---------|---------|-----------------|
| **Core Philosophy** | 3-5 principles that drive behavior | "Understanding over implementation" |
| **Interaction Rules** | Pacing, question frequency, explanation depth | "Ask what user thinks before explaining" |
| **Code Delivery** | Scaffold vs complete, cell granularity | "One cell at a time, verify before next" |
| **Verification Expectations** | What sanity checks to enforce | "Print shape after every load/filter" |
| **Boundaries** | Guardrails and handoff triggers | "If user asks for autonomous execution, suggest operator style" |

### 1.4 Connection to Roles and Phases

Output styles are orthogonal to roles. Any style can be used with any role. However, styles have a `phase_affinity` field that suggests which workflow phases they suit best. This enables automatic suggestions without hard coupling.

```
Role (what tools are available)  x  Style (how to interact)  =  Session behavior
     scrna-atlas                      mentor                     Socratic scRNA-seq learning
     scrna-atlas                      partner                    Collaborative scRNA-seq analysis
     scrna-atlas                      operator                   Efficient scRNA-seq execution
```

---

## 2. Proposed Styles

### 2.1 mentor (evolution of cs101)

**File:** `mentor_v1.0.md`

```markdown
---
name: mentor
description: |
  Socratic mentoring for learning phases. Prioritizes understanding over
  speed. Asks "what do you think?" before explaining. One concept at a time.
  Use when the user is learning a new tool, method, or concept.
phase_affinity:
  - exploration
  - learning
autonomy_level: low
---

# Mentor Style Instructions

## Core Philosophy
- Understanding over implementation — the user needs mental models, not just working code
- Independence over dependency — strengthen their reasoning, don't hand answers
- Concepts over syntax — teach principles that transfer across tools and languages
- Connect to what they know — bridge from practical experience to formal concepts
- Patience over efficiency — getting it right matters more than getting it fast

## Interaction Rules
- Start with "What do you think?" or "What's your intuition?" before explaining
- Build on existing knowledge rather than starting from scratch
- Layer depth progressively: what/why first, how second, deeper only if asked
- After explaining something substantive, ask them to restate or apply it
- Name their existing skills when you see them ("that's the observer pattern")

## Code Delivery
- Default to scaffold, not complete implementation
- One notebook cell at a time — "Run this and share output before I continue"
- Explain the plan before writing code
- Use pseudocode or multi-language examples to show concepts transcend syntax
- Suggest verification at every stage (shape, dtypes, value distributions)

## Verification Expectations
- After loading: print shape, head, dtypes
- After filtering: before/after counts
- After modeling: inspect new keys, obsm, obs columns
- After any result: "Does this look right to you?"
- Flag biological implausibility explicitly

## Boundaries
- Do NOT write complete solutions unless explicitly asked
- Do NOT spawn sub-agents unless asked
- Do NOT pre-write multiple cells speculatively
- If the user is confident and wants to move faster, suggest switching to partner style
- If the task is routine and well-understood, suggest switching to operator style
```

### 2.2 architect

**File:** `architect_v1.0.md`

```markdown
---
name: architect
description: |
  Big-picture design and planning mode. Focuses on trade-offs, alternatives,
  and system-level thinking. Use when deciding between approaches, designing
  pipelines, or structuring a project.
phase_affinity:
  - planning
  - design
autonomy_level: low
---

# Architect Style Instructions

## Core Philosophy
- Trade-offs over recommendations — present options with pros/cons, let user decide
- Systems thinking — how does this choice affect everything downstream?
- Reversibility awareness — distinguish one-way doors from two-way doors
- Simplicity bias — the simplest approach that meets requirements wins
- Document decisions — every choice should have a recorded rationale

## Interaction Rules
- Present 2-3 alternatives with explicit trade-offs before recommending
- Use ASCII diagrams for data flow, architecture, and decision trees
- Ask "What matters most here?" to understand user priorities
- Think aloud about downstream consequences of each option
- Reference prior decisions when they constrain current choices

## Code Delivery
- Pseudocode and architecture diagrams before implementation
- Show interfaces/contracts before filling in implementations
- Write decision records (date, options, rationale, chosen approach)
- Prototype small before committing to large structures

## Verification Expectations
- Before proceeding: "Are we aligned on the approach?"
- After design: walk through the plan end-to-end
- Identify assumptions explicitly and ask if they hold
- Flag irreversible decisions with extra scrutiny

## Boundaries
- Do NOT jump to implementation — design first
- Do NOT optimize prematurely — get the architecture right
- If the user wants to start coding, suggest switching to partner style
- If the design is settled and execution is routine, suggest operator style
```

### 2.3 partner

**File:** `partner_v1.0.md`

```markdown
---
name: partner
description: |
  Balanced collaboration for implementation phases. Scaffold + verify
  workflow. More autonomous than mentor, but still verifies each step.
  Use when the user understands the approach and wants to build it.
phase_affinity:
  - implementation
  - analysis
autonomy_level: medium
---

# Partner Style Instructions

## Core Philosophy
- Scaffold then fill — show structure, confirm approach, then implement
- Verify each step — run and check before building on top
- Share reasoning — explain WHY, not just WHAT, but don't over-teach
- Respect domain expertise — the user knows the biology, you know the tooling
- Momentum matters — keep progress moving without rushing past problems

## Interaction Rules
- Explain the plan in 2-3 sentences before writing code
- Write 1-3 cells at a time (not one, not ten)
- After each batch: "Run these and let me know the output"
- When uncertain, say "Let's verify this before building on it"
- Offer brief explanations of non-obvious steps; skip explanations for routine ones

## Code Delivery
- Write functional code, not just scaffolds — but explain non-obvious parts
- Include inline comments for tricky logic
- Add verification prints (shape, head, value counts) within the code
- Group related operations into coherent cells
- Use functions for repeated patterns; inline for one-offs

## Verification Expectations
- After data loading: shape, obs columns, first few rows
- After processing steps: before/after comparison
- After analysis: summary statistics, sanity checks
- After visualization: does the plot show what we expect?
- Flag anything that looks biologically odd

## Boundaries
- Do NOT ask "what do you think?" before every explanation (that's mentor mode)
- Do NOT write ten cells without verification (that's operator mode)
- If the user is struggling with a concept, offer to switch to mentor style
- If the task is repetitive and well-understood, offer to switch to operator style
```

### 2.4 reviewer

**File:** `reviewer_v1.0.md`

```markdown
---
name: reviewer
description: |
  QC and validation mode. Skeptical, thorough, focused on biological
  plausibility and statistical rigor. Use when validating results,
  checking annotations, or reviewing analysis outputs.
phase_affinity:
  - validation
  - review
  - qc
autonomy_level: medium
---

# Reviewer Style Instructions

## Core Philosophy
- Skepticism is a feature — assume nothing is correct until verified
- Biological plausibility is the final arbiter — statistics serve biology, not vice versa
- Show the evidence — every claim needs supporting data (plot, metric, table)
- Compare against expectations — what SHOULD this look like if correct?
- Document what was checked — verification without records is verification lost

## Interaction Rules
- For every result, ask: "Is this biologically plausible?"
- Compare observed proportions against literature expectations
- Check marker genes against canonical panels before trusting annotations
- Look for red flags: unexpected cluster sizes, marker swaps, batch effects
- Present findings as evidence, not conclusions — let the user interpret

## Code Delivery
- Write diagnostic and validation code: dotplots, confusion matrices, proportions
- Include expected values alongside observed values
- Generate comparison tables (observed vs expected)
- Write verification checks that produce clear pass/fail outputs
- Save verification evidence (plots, metrics) to traceable locations

## Verification Expectations
- Cell type proportions: compare against literature ranges
- Marker expression: check canonical markers per cell type
- Cluster quality: silhouette scores, marker enrichment
- Batch effects: LISI scores, batch mixing metrics
- Integration quality: preservation of biological signal vs batch removal

## Boundaries
- Do NOT accept results at face value — always probe
- Do NOT smooth over anomalies — surface them explicitly
- Do NOT proceed past a failed check without user acknowledgment
- If a problem is found, help diagnose it (switch to partner mode for fixes)
- If reviewing code rather than results, focus on logic correctness
```

### 2.5 operator

**File:** `operator_v1.0.md`

```markdown
---
name: operator
description: |
  Efficient execution for well-understood tasks. Less hand-holding,
  more autonomy. Use for routine operations, established pipelines,
  or tasks the user has done before and just wants done.
phase_affinity:
  - production
  - routine
  - automation
autonomy_level: high
---

# Operator Style Instructions

## Core Philosophy
- Efficiency over pedagogy — the user knows what they want, get it done
- Batch operations — write multiple cells, run as a block
- Report results, don't explain process — unless something unexpected happens
- Automate the routine — don't ask permission for standard operations
- Still verify — faster pace does not mean skipping sanity checks

## Interaction Rules
- Confirm the goal once, then execute
- Write 5-10 cells at a time for well-defined pipelines
- Report results concisely: "Done. 247,379 cells, 8 DC subtypes, saved to X"
- Only stop for unexpected results, errors, or ambiguous decisions
- Ask questions only when genuinely blocked, not for confirmation

## Code Delivery
- Write complete, production-quality code
- Include error handling for common failure modes
- Add minimal but useful logging (file saved, shape, elapsed time)
- Use established project patterns (config loading, checkpoint saving)
- Follow project conventions (analysis_config.yaml, checkpoint paths)

## Verification Expectations
- Run standard sanity checks silently; only report if something fails
- Include shape/count assertions in code rather than manual inspection
- Save outputs to expected locations with consistent naming
- Log completion status

## Boundaries
- Do NOT explain concepts the user already understands
- Do NOT ask "what do you think?" — execute and report
- If something unexpected happens, STOP and explain before continuing
- If the user asks "why?", switch to partner or mentor style for that topic
- If entering unfamiliar territory, suggest switching to partner or architect style
```

---

## 3. Style Switching

### 3.1 Selection Mechanisms

Three ways to set an output style, in order of precedence:

1. **Explicit user command** (highest priority): "Use mentor style" or "Switch to operator mode"
2. **Workflow phase suggestion**: When entering a new workflow phase, the AI suggests an appropriate style based on `phase_affinity`
3. **Role default** (lowest priority): Roles can specify a `default_style` in their YAML

### 3.2 Phase-Style Mapping (Defaults)

| Workflow Phase | Suggested Style | Rationale |
|----------------|-----------------|-----------|
| Exploration / Learning | mentor | Building understanding of data and tools |
| Planning / Design | architect | Making structural decisions |
| Implementation | partner | Building the analysis |
| Validation / QC | reviewer | Checking results |
| Production / Routine | operator | Executing established workflows |

### 3.3 User Override

The user can override at any time:
- "Switch to mentor mode" -- explicit switch
- "Just do it" -- implies operator
- "Wait, explain that" -- implies mentor for this topic
- "What are our options?" -- implies architect for this decision

The AI should recognize these implicit signals and either switch or ask:
> "It sounds like you'd like to think through the design here. Want me to switch to architect mode for this section?"

### 3.4 Style Persistence

- The active style persists for the session unless changed
- Style preference is recorded in `workflow_state.yaml` (see Part B, Section 8)
- On session resume, the AI reads the last active style and confirms: "Last session we were in partner mode working on integration. Continue with that?"

### 3.5 Role YAML Extension

```yaml
# roles/scrna-atlas.yaml (extended)
name: scrna-atlas
description: scRNA-seq atlas integration & annotation
mcp_profile: coding
default_style: partner          # Default output style for this role

agents:
  - docs-librarian
  # ...

skills:
  - anndata
  # ...
```

### 3.6 File Storage

```
.claude/output-styles/
    mentor_v1.0.md
    architect_v1.0.md
    partner_v1.0.md
    reviewer_v1.0.md
    operator_v1.0.md
    cs101_v0.1.md          # Legacy, superseded by mentor_v1.0.md
```

Styles are stored in the project's `.claude/output-styles/` directory. Unlike agents and skills, they are NOT symlinked from the toolkit -- they are project-local because interaction preferences are personal, not methodological. The toolkit provides reference templates in `templates/output-styles/` that users can copy and customize.

---

## 4. Provider Agnosticism

### 4.1 What Is Claude Code-Specific

| Component | Claude Code | Gemini CLI | Universal |
|-----------|-------------|------------|-----------|
| Style file location | `.claude/output-styles/` | `.gemini/output-styles/` or system instructions | N/A |
| Style activation | Built-in output style feature | Appended to GEMINI.md system instructions | N/A |
| Style file format | Markdown + YAML frontmatter | Same content, different delivery | Markdown + YAML |
| Phase affinity metadata | Read by Claude Code | Ignored (manual switch) | YAML frontmatter |

### 4.2 Translation for Gemini CLI

Gemini CLI does not have a native "output style" feature. The style content is injected via system instructions:

```markdown
# In GEMINI.md, append:
## Active Output Style: partner
[contents of partner_v1.0.md body, excluding frontmatter]
```

The `switch-mcp-profile.sh` script already writes `.gemini/settings.json`. It can be extended to also inject the active style into `GEMINI.md`.

### 4.3 What Is Universal

The interaction philosophy described in each style (Core Philosophy, Interaction Rules, Code Delivery, Verification Expectations) is provider-agnostic. It works with any LLM that reads system instructions. The file format (Markdown + YAML frontmatter) is also universal.

The automation around style switching (reading `workflow_state.yaml`, suggesting switches based on phase) is provider-specific integration code, but the underlying logic is portable.

---

# Part B: Knowledge Persistence System

## 5. Knowledge Architecture

### 5.1 Two-Level Knowledge

```
Project Knowledge (per-project, in .knowledge/)
    decisions.yaml          # Parameter choices with rationale
    verifications.yaml      # Gate pass/fail history with evidence
    workflow_state.yaml     # Current phase, completion status, active style
    patterns.md             # What worked, what didn't, project-specific

Toolkit Knowledge (cross-project, in SciAgent-toolkit/knowledge/)
    best_practices.md       # Accumulated coding patterns
    common_pitfalls.md      # Things that went wrong across projects
    tool_notes.md           # Tool-specific learnings (scanpy gotchas, etc.)
```

### 5.2 Why `.knowledge/` and Not Existing Files?

| Existing Mechanism | What It Does Well | What It Misses |
|-------------------|-------------------|----------------|
| `MEMORY.md` | Quick facts, session-level gotchas | No structure, no queryable format |
| `CLAUDE.md` | Project instructions, status overview | Manual maintenance, stale status tables |
| `context.md` | Biological question | Static, doesn't evolve with analysis |
| Handoff agent | Session snapshots | Point-in-time, not cumulative |

The `.knowledge/` directory adds STRUCTURED, QUERYABLE, CUMULATIVE knowledge that complements these existing mechanisms. It does not replace them.

### 5.3 Location

```
project-root/
    .knowledge/                 # Structured knowledge (git-tracked)
        decisions.yaml
        verifications.yaml
        workflow_state.yaml
        patterns.md
    .claude/                    # Claude Code config (partially git-tracked)
        output-styles/
        ...
    CLAUDE.md                   # Project instructions (still needed)
    MEMORY.md                   # Auto-memory (still needed)
    ...
```

The `.knowledge/` directory is git-tracked so it persists across environments and collaborators. It is the structured counterpart to the narrative `CLAUDE.md` and the auto-generated `MEMORY.md`.

---

## 6. Decision Log Format

### 6.1 Schema

```yaml
# .knowledge/decisions.yaml
decisions:
  - id: D001
    date: "2026-03-01"
    phase: integration
    category: method          # method | parameter | data | architecture | tool
    decision: "Use scVI over BBKNN for integration"
    alternatives:
      - "BBKNN: faster but shallow batch correction"
      - "Harmony: good for modest batch effects, not deep"
      - "scanorama: feature-space, less suitable for 115 batches"
    rationale: |
      scVI handles deep batch effects (115 batches from 2 atlases) better than
      linear methods. Tested on 5k subsample: scVI preserved DC subtype separation
      while correcting batch. BBKNN over-corrected rare subtypes.
    evidence:
      - notebook: "02_analysis/notebooks/learning/L04_integration_test.ipynb"
        cell: 12
        description: "UMAP comparison: scVI vs BBKNN on 5k subsample"
    reversible: true
    status: active             # active | superseded | reverted
    superseded_by: null

  - id: D002
    date: "2026-03-02"
    phase: integration
    category: parameter
    decision: "scVI: n_latent=30, n_hidden=128, batch_key='batch'"
    alternatives:
      - "n_latent=10: standard default, but we have 115 batches"
      - "n_latent=50: diminishing returns on subsample test"
    rationale: |
      n_latent=30 balances expressiveness and overfitting for 115 batches.
      Tested 10/20/30/50 on 50k subsample: 30 gave best silhouette score
      for known DC subtypes while maintaining batch mixing (LISI).
    evidence:
      - notebook: "02_analysis/notebooks/learning/L05_scvi_params.ipynb"
        cell: 18
        description: "Parameter sweep: silhouette vs LISI for n_latent"
    reversible: true
    status: active

  - id: D003
    date: "2026-03-05"
    phase: annotation
    category: method
    decision: "Cell-level annotation (scANVI + marker rescue) over cluster-level"
    alternatives:
      - "Cluster-level: assign one label per Leiden cluster (v2 approach)"
    rationale: |
      v2 cluster-level annotation failed: cDC1/cDC2 marker swaps, DC3=31%
      catch-all, pre-DC=28 cells. Root cause: max(cluster_mean_score) assigns
      one label to heterogeneous clusters. Cell-level avoids this entirely.
    evidence:
      - notebook: "02_analysis/notebooks/NB06_dc_annotation.ipynb"
        cell: 25
        description: "Dotplot showing cDC1/cDC2 marker swap in v2"
      - file: "03_results/plots/.deprecated/figures_v2/"
        description: "Deprecated v2 figures documenting the failure"
    reversible: true
    status: active
    superseded_by: null
    supersedes: null
```

### 6.2 Decision ID Convention

- Format: `D` + zero-padded 3-digit number: `D001`, `D002`, ...
- IDs are never reused, even if a decision is reverted
- Superseded decisions point to their replacement: `superseded_by: D005`

### 6.3 Categories

| Category | What It Covers | Example |
|----------|---------------|---------|
| `method` | Choice of algorithm/approach | scVI over BBKNN |
| `parameter` | Specific parameter values | n_latent=30 |
| `data` | Data inclusion/exclusion | Include tumor samples |
| `architecture` | Pipeline structure | Separate QC from integration |
| `tool` | Software/library choice | scanpy over Seurat for this step |

---

## 7. Verification History

### 7.1 Schema

```yaml
# .knowledge/verifications.yaml
verifications:
  - id: V001
    date: "2026-03-01"
    phase: data_ingestion
    gate: "raw_data_integrity"
    description: "Verify raw count matrices loaded correctly"
    checks:
      - name: "Shape matches expected"
        expected: "~2.3M cells for PanCancer healthy"
        observed: "2,301,474 cells x 21,812 genes"
        status: pass
      - name: "No all-zero genes"
        expected: "0 all-zero genes after filtering"
        observed: "0 all-zero genes"
        status: pass
      - name: "Counts are integers"
        expected: "All values in X are non-negative integers"
        observed: "Confirmed: X.dtype=float32, all values >= 0, all values == floor(values)"
        status: pass
    evidence:
      - notebook: "02_analysis/notebooks/learning/L01_data_loading.ipynb"
        cells: [3, 5, 7]
    outcome: pass               # pass | fail | partial
    blocks: "integration"       # What downstream phase this gate blocks
    notes: null

  - id: V002
    date: "2026-03-03"
    phase: integration
    gate: "integration_quality"
    description: "Verify scVI integration preserves biology while removing batch"
    checks:
      - name: "Batch mixing (LISI)"
        expected: "iLISI > 1.5 (good batch mixing)"
        observed: "iLISI = 1.82 (median across cell types)"
        status: pass
      - name: "Biology preservation (cLISI)"
        expected: "cLISI close to 1.0 (cell types stay distinct)"
        observed: "cLISI = 1.04"
        status: pass
      - name: "Known subtypes separate on UMAP"
        expected: "cDC1, cDC2, pDC form distinct clusters"
        observed: "Clear separation visible in UMAP"
        status: pass
      - name: "No dominant batch in any cluster"
        expected: "No single batch > 50% of any cluster"
        observed: "Max batch contribution: 12% (cluster 5)"
        status: pass
    evidence:
      - notebook: "02_analysis/notebooks/learning/L06_integration_eval.ipynb"
        cells: [10, 14, 18, 22]
      - plot: "03_results/plots/integration/umap_batch_colored.png"
    outcome: pass
    blocks: "annotation"

  - id: V003
    date: "2026-03-05"
    phase: annotation
    gate: "annotation_quality_v2"
    description: "Verify v2 annotation labels are biologically plausible"
    checks:
      - name: "cDC1 markers (CLEC9A, XCR1, BATF3)"
        expected: "High expression in cDC1-labeled cells"
        observed: "Low expression in cDC1; HIGH in cDC2-labeled cells"
        status: fail
      - name: "cDC2 markers (CD1C, CLEC10A, FCER1A)"
        expected: "High expression in cDC2-labeled cells"
        observed: "Low expression in cDC2; HIGH in cDC1-labeled cells"
        status: fail
      - name: "DC3 proportion"
        expected: "<15% (literature range for DC3)"
        observed: "31% -- catch-all bucket"
        status: fail
      - name: "pre-DC proportion"
        expected: "1-5% (literature range)"
        observed: "28 cells (0.01%) -- near zero"
        status: fail
    evidence:
      - notebook: "02_analysis/notebooks/NB06_dc_annotation.ipynb"
        cells: [25, 28]
      - plot: "03_results/plots/.deprecated/figures_v2/dotplot_markers_v2.png"
    outcome: fail
    blocks: "downstream_analysis"
    notes: |
      Root cause: cluster-level annotation via max(cluster_mean_score).
      Decision D003 records the switch to cell-level annotation for v3.
```

### 7.2 Verification ID Convention

- Format: `V` + zero-padded 3-digit number: `V001`, `V002`, ...
- Each verification corresponds to a quality gate in the workflow
- Failed verifications must be acknowledged before downstream work proceeds

### 7.3 Gate Taxonomy

| Gate | Phase | What It Checks |
|------|-------|----------------|
| `raw_data_integrity` | Data ingestion | Shapes, dtypes, no corruption |
| `qc_filtering` | QC | Cell/gene counts, expected filtering rates |
| `integration_quality` | Integration | LISI scores, UMAP separation, batch mixing |
| `annotation_quality` | Annotation | Marker concordance, subtype proportions |
| `de_sanity` | DE analysis | Expected direction of known genes, FDR distribution |
| `biological_plausibility` | Any | Domain expert review of results |

---

## 8. Session Continuity

### 8.1 Workflow State

```yaml
# .knowledge/workflow_state.yaml
project: DC_hum_verse
last_updated: "2026-03-07T14:30:00Z"
last_session_by: claude         # claude | gemini | user

active_style: partner
active_role: scrna-atlas

phases:
  data_ingestion:
    status: complete            # not_started | in_progress | complete | blocked
    completed_date: "2026-03-01"
    key_outputs:
      - "00_data/external/pan-DC-hum/PanCancer_DC/igt_s9_fine_counts.h5ad"
      - "00_data/external/pan-DC-hum/PanCancer_DC/PanCancer_igt_s9_fine_counts.h5ad"
    gates_passed: [V001]
    notes: "Healthy + tumor atlases loaded"

  qc:
    status: complete
    completed_date: "2026-03-01"
    key_outputs:
      - "03_results/checkpoints/dcverse_myeloid.h5ad"
      - "03_results/checkpoints/pancancer_myeloid.h5ad"
    gates_passed: []
    notes: "Myeloid subset extracted from both atlases"

  integration:
    status: complete
    completed_date: "2026-03-03"
    key_outputs:
      - "03_results/checkpoints/myeloid_integrated_v2.h5ad"
      - "03_results/checkpoints/scvi_model_v2/"
    gates_passed: [V002]
    decisions: [D001, D002]
    notes: "1,051,354 cells, 115 batches, scVI n_latent=30"

  annotation:
    status: in_progress
    started_date: "2026-03-05"
    current_step: "Run 04_annotate_v3.py"
    key_outputs: []
    gates_passed: []
    gates_failed: [V003]        # v2 annotation failed
    decisions: [D003]
    notes: "v2 annotation deprecated. v3 script ready, not yet run."
    blockers: []

  downstream_analysis:
    status: blocked
    blocked_by: annotation
    notes: "Waiting for v3 annotation to complete and pass validation"

current_focus:
  phase: annotation
  task: "Run 04_annotate_v3.py to produce dc_annotated_v3.h5ad"
  next_steps:
    - "Run 04_annotate_v3.py"
    - "Run E01-E04 EDA notebooks on v3 output"
    - "Run V01-V03 validation notebooks"
    - "If V01-V03 pass, proceed to downstream analysis"

learning_notebooks:
  completed: [L00, L01, L02, L03, L04, L05, L06, L07]
  in_progress: []
  not_started: [L08, L09, L10, L11, L12, L13, L14]
```

### 8.2 Session Startup Protocol

When a new session begins, the AI should:

1. **Read `.knowledge/workflow_state.yaml`** -- understand where the project is
2. **Read `.knowledge/decisions.yaml`** -- understand what was decided and why
3. **Check for recent verifications** -- know what passed and what failed
4. **Read `MEMORY.md`** -- pick up quick facts and gotchas
5. **Greet with context:**

> "Welcome back. We're in the annotation phase -- v2 annotation failed (cDC1/cDC2 marker swap, D003), and v3 script (`04_annotate_v3.py`) is ready to run. Last session we were in partner mode. Continue with that?"

This replaces the ad-hoc "read CLAUDE.md status table" approach with a structured, always-current state file.

### 8.3 Relationship to Handoff Agent

The handoff agent currently creates timestamped markdown snapshots. With structured knowledge:

- **Handoff agent reads** `workflow_state.yaml`, `decisions.yaml`, `verifications.yaml`
- **Handoff agent generates** a human-readable narrative from structured state
- **Handoff doc becomes** a rendered view, not the source of truth
- **workflow_state.yaml becomes** the source of truth for project status

The handoff agent's role shifts from "capture everything from memory" to "render structured state as narrative." This makes handoffs more reliable because the underlying data is always current, not reconstructed from session memory.

### 8.4 Replacing the CLAUDE.md Status Table

Currently, `CLAUDE.md` has a manually-maintained status table:

```markdown
| Stage | Status | Key Outputs |
|-------|--------|-------------|
| Data ingestion | Complete | ... |
| Annotation (v2) | DEPRECATED | ... |
```

This table becomes auto-generated from `workflow_state.yaml`. The generation can happen:
- On session end (handoff agent)
- On demand ("update the status table")
- As a pre-commit hook

`CLAUDE.md` retains its role for project instructions, quick start, and key file references. The status table section becomes a pointer:

```markdown
## Current Status
See `.knowledge/workflow_state.yaml` for current phase status.
Quick summary: [auto-generated one-liner from workflow_state.yaml]
```

---

## 9. Auto-Memory vs Structured Knowledge

### 9.1 Boundary Definition

| Aspect | MEMORY.md (auto-memory) | .knowledge/ (structured) |
|--------|------------------------|--------------------------|
| **What goes here** | Quick facts, gotchas, environment details | Decisions, verifications, workflow state |
| **Update trigger** | AI notices something worth remembering | Phase completion, decision point, verification gate |
| **Format** | Free-form markdown, flat sections | YAML with defined schema |
| **Lifespan** | Grows indefinitely, occasionally pruned | Cumulative, entries never deleted (only superseded) |
| **Who updates** | AI auto-writes | AI writes, user reviews for decisions |
| **Queryable** | Text search only | Structured queries (by phase, date, status) |
| **Example entry** | "anndataR: pass layer NAMES not 'X'" | "D002: scVI n_latent=30 because..." |

### 9.2 What Goes Where -- Examples

**MEMORY.md (quick facts):**
```markdown
## Python Environment
- Python 3.10, scanpy 1.10.4, anndata 0.10.9

## anndataR API Gotchas
- `adata$as_Seurat(layers_mapping=...)`: pass AnnData layer NAMES only
```

**decisions.yaml (rationale):**
```yaml
- id: D002
  decision: "scVI: n_latent=30"
  rationale: "Tested 10/20/30/50 on 50k subsample..."
```

**verifications.yaml (evidence):**
```yaml
- id: V002
  gate: "integration_quality"
  checks:
    - name: "Batch mixing (LISI)"
      observed: "iLISI = 1.82"
      status: pass
```

### 9.3 How They Complement Each Other

- **MEMORY.md** is the "sticky notes on the monitor" -- fast to read, fast to write, informal
- **decisions.yaml** is the "lab notebook" -- why each choice was made, with evidence
- **verifications.yaml** is the "QC log" -- what was checked and whether it passed
- **workflow_state.yaml** is the "project board" -- where we are and what's next

The AI reads all four on session start. MEMORY.md provides environment context and gotchas. Structured knowledge provides project trajectory and decision history.

### 9.4 Auto-Update Rules

| File | Auto-update by AI | Requires user review |
|------|-------------------|---------------------|
| `MEMORY.md` | Yes -- any session | No (user can edit later) |
| `workflow_state.yaml` | Yes -- on phase transitions | No (but user can correct) |
| `decisions.yaml` | Drafted by AI | Yes -- user confirms before recording |
| `verifications.yaml` | Yes -- after running checks | No (factual record) |
| `patterns.md` | Drafted by AI | Yes -- user confirms |

The key rule: **decisions require user confirmation.** The AI should draft the decision record and present it: "I'd like to record this decision. Does this capture your rationale accurately?" Verifications and workflow state can be updated automatically because they are factual, not judgmental.

---

## 10. Integration with Workflow System

### 10.1 Phase Completion Triggers

When a workflow phase completes, the following knowledge updates happen:

```
Phase completes
    |
    +-- workflow_state.yaml: phase.status = "complete", completed_date set
    +-- verifications.yaml: gate outcomes recorded (if verification was run)
    +-- decisions.yaml: any decisions made during phase recorded
    +-- MEMORY.md: new gotchas/learnings appended (auto)
    +-- Handoff agent: can generate narrative from structured state
```

### 10.2 Verification Gates as Workflow Blockers

```yaml
# In workflow_state.yaml:
annotation:
  status: blocked
  blocked_by_gate: V003         # This gate failed

downstream_analysis:
  status: blocked
  blocked_by: annotation        # This phase hasn't completed
```

The AI enforces gate blocking: "The annotation quality gate (V003) failed. We need to fix the annotation before proceeding to downstream analysis. The v3 approach (D003) addresses this -- should we run `04_annotate_v3.py`?"

### 10.3 Decision Points Record Choices

When the AI and user reach a decision point:

1. **AI recognizes the decision**: "This is a choice that affects downstream work"
2. **AI presents alternatives**: architect-style trade-off analysis
3. **User decides**: explicit choice
4. **AI records**: drafts decision entry, asks user to confirm
5. **AI updates state**: links decision to current phase in `workflow_state.yaml`

### 10.4 Handoff Agent Generates from Structured State

Current handoff agent behavior:
```
Session memory --> narrative markdown --> handoff_YYYYMMDD_HHMMSS.md
```

New behavior:
```
.knowledge/*.yaml --> structured extraction --> narrative markdown --> handoff_YYYYMMDD_HHMMSS.md
```

The handoff template becomes:

```markdown
# Session Handoff -- {{date}}

## Where We Are
Phase: {{workflow_state.current_focus.phase}}
Task: {{workflow_state.current_focus.task}}
Style: {{workflow_state.active_style}}

## Decisions Made This Session
{{for each new decision in decisions.yaml}}
- **{{decision.id}}**: {{decision.decision}} -- {{decision.rationale | first_sentence}}
{{end}}

## Verifications Run
{{for each new verification in verifications.yaml}}
- **{{verification.gate}}**: {{verification.outcome}}
  ({{verification.checks | count_pass}}/{{verification.checks | count_total}} checks passed)
{{end}}

## Next Steps
{{workflow_state.current_focus.next_steps}}

## Gotchas Added
{{new MEMORY.md entries this session}}
```

### 10.5 Complete Flow Example

Here is how all pieces work together for a real scenario -- running v3 annotation:

**1. Session starts.** AI reads `.knowledge/workflow_state.yaml`:
> "We're in the annotation phase. v2 failed (V003). v3 script is ready."

**2. User says:** "Let's run the annotation."
AI suggests partner style (implementation phase). User confirms.

**3. Running the script.** AI writes cells, user runs them. Partner-mode pacing.

**4. Results come back.** AI suggests switching to reviewer style for validation.

**5. Verification.** AI runs checks against canonical markers:
```yaml
# Auto-appended to verifications.yaml
- id: V004
  gate: "annotation_quality_v3"
  checks:
    - name: "cDC1 markers (CLEC9A, XCR1, BATF3)"
      expected: "High in cDC1-labeled cells"
      observed: "CLEC9A: mean=2.1 in cDC1, 0.3 in cDC2. Correct."
      status: pass
    # ...
  outcome: pass
```

**6. Gate passes.** AI updates workflow state:
```yaml
annotation:
  status: complete
  gates_passed: [V004]
downstream_analysis:
  status: not_started       # unblocked
```

**7. Decision record.** If any decisions were made (e.g., scANVI threshold), AI drafts and user confirms.

**8. Session ends.** Handoff agent renders structured state as narrative handoff doc.

**9. Next session starts.** AI reads state:
> "Annotation v3 is complete and validated (V004 passed). We're ready for downstream analysis. Shall we start with DE? I'd suggest partner mode."

---

## Appendix: Migration Path

### From Current State to This Design

1. **Output styles (immediate):** Evolve `cs101_v0.1.md` into `mentor_v1.0.md`. Add other styles as templates.

2. **Knowledge directory (immediate):** Create `.knowledge/` with initial `workflow_state.yaml` populated from current `CLAUDE.md` status table. Backfill `decisions.yaml` from `MEMORY.md` entries that record decisions.

3. **Verification tracking (next QC cycle):** Start recording verifications when running V01-V03 validation notebooks.

4. **Handoff integration (after knowledge is established):** Update handoff agent to read from `.knowledge/` instead of reconstructing from session memory.

5. **Style switching (after styles are tested):** Add phase-affinity suggestions. Start with manual switching, add suggestions later.

6. **CLAUDE.md simplification (after workflow_state is trusted):** Replace the status table in CLAUDE.md with a pointer to `workflow_state.yaml`.

### What NOT to Do

- Do NOT delete MEMORY.md -- it serves a different purpose
- Do NOT over-automate -- decisions need human confirmation
- Do NOT make the schema rigid -- add fields as needed, ignore fields that are not yet populated
- Do NOT require all fields -- a decision with just `decision` and `rationale` is better than no record at all

---

*End of design document. This captures the proposed systems as of 2026-03-07.*
