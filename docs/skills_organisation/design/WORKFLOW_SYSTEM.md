# Workflow System and Verification Gates

> **Date:** 2026-03-07
> **Status:** Design proposal
> **Addresses:** Gap Analysis items "No phase-based workflow skills" and "No verification/checkpoint skills" from RESEARCH_20260307.md
> **Informed by:** life-sciences waypoint pattern, nextflow validation gates, clinical-trial-protocol staged execution

---

## 1. Workflow Definition Format

A workflow is a YAML file in `workflows/` that declares the full phase progression for a project type. It is **not executable** -- it is a contract that the AI reads to know where the user is, what comes next, and what must be true before advancing.

```yaml
# workflows/scrna-atlas.yaml
name: scrna-atlas
description: scRNA-seq atlas integration and annotation
version: 1.0

# Where workflow state lives in the project
state_file: 02_analysis/config/workflow_state.yaml

# Phases in order
phases:
  - id: explore
    name: Independent Dataset Exploration
    description: Explore each dataset independently -- QC, structure, contents
    prerequisites: []
    skills: [anndata, scanpy, single-cell-rna-qc]
    output_style: cs-101
    scale: full  # no subsample needed, datasets are explored as-is
    verification_gates:
      - id: datasets_explored
        description: Each dataset explored with basic QC
        type: checklist
        items:
          - "Cell counts and feature counts documented per dataset"
          - "Annotation columns examined (what labels exist?)"
          - "HVG overlap checked between datasets"
          - "Basic UMAP/distribution plots generated"
      - id: data_understanding
        description: User can articulate what each dataset contributes
        type: human_judgment
        prompt: "Can you summarize what each dataset brings to the atlas?"
    decision_points:
      - id: datasets_to_include
        description: Which datasets to carry forward
        prompt: "Based on exploration, which datasets will you integrate? Any to drop?"
        must_document: true
    checkpoints:
      - description: Exploration summary per dataset
        suggested_format: markdown_in_notebook

  - id: integration_planning
    name: Integration Planning
    description: Decide what to integrate, what features to use, what cells to subset
    prerequisites: [explore.datasets_explored, explore.data_understanding]
    skills: [anndata, scanpy]
    output_style: cs-101
    scale: full
    verification_gates:
      - id: subset_validated
        description: Subsetting targets the right cells
        type: checklist
        items:
          - "Marker genes confirm target population identity"
          - "Cell counts after subsetting are reasonable"
          - "No garbage annotation leaking into subset"
    decision_points:
      - id: integration_scope
        description: What cell types to integrate
        prompt: "Integrate all myeloid? Just DCs? Include tumor cells for correlation?"
        must_document: true
      - id: batch_definition
        description: What constitutes a batch
        prompt: "What is your batch key? Dataset? Donor? Sample?"
        must_document: true
    checkpoints:
      - description: Subsetted object saved
        suggested_format: h5ad

  - id: integration_understanding
    name: Integration Understanding (Small Scale)
    description: Learn integration methods on subsamples before committing
    prerequisites: [integration_planning.subset_validated]
    skills: [scvi-basic, scvi-framework, scanpy]
    output_style: cs-101
    scale: small
    scale_spec:
      strategy: random_subsample
      target_cells: 5000-20000
      must_preserve: "all batches represented"
    verification_gates:
      - id: method_understood
        description: User understands what integration does and doesn't change
        type: human_judgment
        prompt: |
          Can you answer these:
          1. What does scVI change -- expression values? embeddings? both?
          2. Are raw counts modified after integration?
          3. If you integrate myeloid, can you still correlate with tumor cells?
          4. If you subset DCs from integrated myeloid, do you re-integrate?
      - id: small_run_succeeded
        description: Integration ran successfully on subsample
        type: checklist
        items:
          - "scVI/scANVI trained on subsample without errors"
          - "UMAP of latent space shows batch mixing"
          - "Known biology preserved (cell types still separate)"
    decision_points:
      - id: method_choice
        description: Which integration method to use
        prompt: "scVI, BBKNN, Harmony, or scANVI? Based on what you learned."
        must_document: true
      - id: parameters
        description: Key parameters
        prompt: "n_latent, n_hidden, batch_key -- what values and why?"
        must_document: true
    checkpoints:
      - description: Subsample integrated object + model
        suggested_format: h5ad + model_dir

  - id: full_integration
    name: Full Integration
    description: Run chosen method on full data with validated parameters
    prerequisites:
      - integration_understanding.method_understood
      - integration_understanding.small_run_succeeded
    skills: [scvi-basic, scvi-framework, scanpy]
    output_style: cs-101
    scale: full
    verification_gates:
      - id: integration_qc
        description: Full integration passes quality checks
        type: checklist
        items:
          - "UMAP shows batch mixing (no dataset islands)"
          - "Known marker genes still separate expected populations"
          - "Original annotations map to plausible UMAP regions"
          - "Cell counts match pre-integration object"
    checkpoints:
      - description: Integrated object + trained model
        suggested_format: h5ad + model_dir
        save_as: "myeloid_integrated.h5ad"

  - id: annotation_manual
    name: Manual Annotation
    description: Build annotation understanding through feature plots and marker genes
    prerequisites: [full_integration.integration_qc]
    skills: [scanpy, anndata]
    output_style: cs-101
    scale: full
    verification_gates:
      - id: markers_validated
        description: Marker genes match expected patterns
        type: checklist
        items:
          - "Feature plots show expected marker patterns per cluster"
          - "Dotplot confirms marker specificity (no swaps)"
          - "Known rare populations are identifiable"
          - "No catch-all cluster >25% of cells"
      - id: annotation_confidence
        description: User is confident in manual labels
        type: human_judgment
        prompt: >
          Are you confident these manual annotations are correct?
          Any ambiguous clusters?
    decision_points:
      - id: annotation_granularity
        description: How fine-grained to annotate
        prompt: "Broad types only (cDC1, cDC2, pDC)? Or subtypes (cDC2A, cDC2B)?"
        must_document: true
      - id: ambiguous_clusters
        description: How to handle unclear clusters
        prompt: "What to do with clusters that don't match known markers?"
        must_document: true
    checkpoints:
      - description: Manually annotated object
        suggested_format: h5ad
        save_as: "dc_annotated_manual.h5ad"

  - id: annotation_semisupervised
    name: Semi-Supervised Annotation
    description: Test automated methods against manual annotations
    prerequisites:
      - annotation_manual.markers_validated
      - annotation_manual.annotation_confidence
    skills:
      - scvi-scanvi
      - scvi-scarches-reference-mapping
      - cellxgene-census-annotation
    output_style: cs-101
    scale: small_then_full
    scale_spec:
      strategy: subsample_first
      small_target: 10000-30000
      advance_when: "small-scale results match manual annotation >90%"
    verification_gates:
      - id: matches_manual
        description: Automated labels agree with manual on known types
        type: checklist
        items:
          - "Confusion matrix between manual and automated labels"
          - "Agreement >90% on well-defined types"
          - "Marker dotplot for automated labels shows correct patterns"
          - "Rare types (pDC, AS-DC) correctly identified"
      - id: no_marker_swaps
        description: No systematic marker swaps (the v2 failure mode)
        type: checklist
        items:
          - "cDC1 markers (CLEC9A, XCR1) enriched in cDC1 label"
          - "cDC2 markers (CD1C, FCER1A) enriched in cDC2 label"
          - "No catch-all type >25%"
    decision_points:
      - id: final_annotation_source
        description: Which annotations to use going forward
        prompt: "Use manual, automated, or merged annotations?"
        must_document: true
    checkpoints:
      - description: Final annotated object
        suggested_format: h5ad
        save_as: "dc_annotated_final.h5ad"

  - id: deep_analysis
    name: Deep Analysis
    description: DE, GSEA, iron biology -- only after annotation confidence
    prerequisites:
      - annotation_semisupervised.matches_manual
      - annotation_semisupervised.no_marker_swaps
    skills: [scanpy, anndata, genenmf-metaprogram-discovery]
    output_style: null  # user's preference at this point
    scale: full
    verification_gates: []  # analysis-specific, defined per question
    decision_points:
      - id: analysis_questions
        description: What biological questions to answer
        prompt: "What are your specific hypotheses? What comparisons matter?"
        must_document: true
```

### Format Rationale

- **YAML, not code.** Workflows are contracts, not programs. The AI reads the YAML and knows what phase the user is in, what skills to load, and what gates to check.
- **Flat phase list with prerequisite references.** No nested trees. Prerequisites point to specific gate IDs using `phase_id.gate_id` syntax. This makes it trivial to check "can I advance?"
- **Human-readable.** A user can open this file, scan it, and understand the project progression without running anything.

---

## 2. Phase Structure

Each phase has these fields:

| Field | Required | Purpose |
|-------|----------|---------|
| `id` | Yes | Machine-readable identifier, used in prerequisite references |
| `name` | Yes | Human-readable name |
| `description` | Yes | What this phase accomplishes (1-2 sentences) |
| `prerequisites` | Yes | List of `phase_id.gate_id` that must be passed. Empty list `[]` for the first phase |
| `skills` | Yes | Which skills the AI should have loaded for this phase |
| `output_style` | No | Override interaction style (e.g., `cs-101` for pedagogical, `null` for default) |
| `scale` | Yes | One of: `full`, `small`, `small_then_full` |
| `scale_spec` | No | Details when `scale` is `small` or `small_then_full` (see Section 4) |
| `verification_gates` | Yes | List of gates (see Section 3). Can be empty `[]` for open-ended phases |
| `decision_points` | No | List of decisions requiring human input (see Section 5) |
| `checkpoints` | No | What to save at phase completion (see Section 6) |

### How the AI Uses Phase Structure

When a session starts, the AI:
1. Reads `workflow_state.yaml` to find the current phase
2. Reads the workflow YAML to get that phase's definition
3. Loads the listed skills
4. Adopts the output style
5. Respects the scale constraint
6. When the user seems ready to advance, presents the verification gates
7. Does NOT auto-advance -- waits for the user to confirm gate passage

---

## 3. Verification Gate Design

### Gate Types

```yaml
# Type 1: Checklist -- concrete items that can be checked
- id: datasets_explored
  type: checklist
  items:
    - "Cell counts and feature counts documented per dataset"
    - "Annotation columns examined"

# Type 2: Human judgment -- requires the user to answer a question
- id: data_understanding
  type: human_judgment
  prompt: "Can you summarize what each dataset brings?"

# Type 3: Computational -- a code snippet that returns pass/fail
- id: cell_counts_match
  type: computational
  check: |
    assert adata.n_obs == expected_n_obs, f"Expected {expected_n_obs}, got {adata.n_obs}"
```

### How Gates Are Presented

When the user is approaching the end of a phase (or explicitly asks "am I ready to move on?"), the AI presents gates as a checklist:

```
--- Phase Gate: Independent Dataset Exploration ---

Before moving to Integration Planning, let's verify:

Checklist: datasets_explored
  [ ] Cell counts and feature counts documented per dataset
  [ ] Annotation columns examined (what labels exist?)
  [ ] HVG overlap checked between datasets
  [ ] Basic UMAP/distribution plots generated

Human judgment: data_understanding
  Can you summarize what each dataset brings to the atlas?

Which items are done? Any you want to revisit?
```

### Anti-Auto-Pass Rules

The AI must follow these rules to prevent rubber-stamping gates:

1. **Never mark a checklist item as passed without evidence.** The AI must point to a specific notebook cell, plot, or output that demonstrates the item. "I believe we did this" is not sufficient.

2. **Never answer a human_judgment gate on behalf of the user.** The AI presents the prompt and waits. If the user says "yeah it's fine," the AI accepts that. But the AI never says "I think you understand this, so we can move on."

3. **Log all gate outcomes.** Every gate pass/fail is recorded in `workflow_state.yaml` with a timestamp and brief evidence note.

4. **Failed gates are information, not blockers.** If an item isn't done, the AI suggests what to do, but does not refuse to proceed if the user explicitly chooses to skip. The skip is logged.

### Gate Outcome Logging

In `workflow_state.yaml`:

```yaml
gates:
  explore.datasets_explored:
    status: passed
    timestamp: 2026-03-07T14:30:00
    evidence:
      - "Cell counts in L01 cell 12 output"
      - "Annotation columns in L01 cell 15"
      - "HVG overlap in L02 cell 8"
      - "UMAPs in L02 cells 10-12"
  explore.data_understanding:
    status: passed
    timestamp: 2026-03-07T14:35:00
    evidence:
      - "User summarized: DC-VERSE provides healthy DC atlas, PanCancer adds tumor context"
  integration_planning.subset_validated:
    status: skipped
    timestamp: 2026-03-08T09:00:00
    reason: "User chose to proceed without marker validation -- will revisit after integration"
```

---

## 4. Small-Sample Testing Pattern

### Why This Exists

The user's workflow explicitly calls for "test on SMALL subsamples first" before committing to expensive full-scale runs. This prevents wasting hours on a bad parameter choice and builds understanding incrementally.

### Scale Values

| Value | Meaning |
|-------|---------|
| `full` | Use the entire dataset. No subsample needed. |
| `small` | Work only on a subsample. The phase is about learning, not producing final results. |
| `small_then_full` | Start with subsample, graduate to full when small-scale results are satisfactory. |

### Scale Spec

```yaml
scale_spec:
  strategy: random_subsample          # or: stratified, per_batch
  target_cells: 5000-20000            # range, not exact
  must_preserve: "all batches represented"  # constraint on subsample
  advance_when: "small-scale results match manual annotation >90%"  # graduation criterion
```

### How the AI Uses Scale Constraints

When `scale: small`:
- The AI suggests subsample code at the start of the phase
- The AI reminds the user that results are preliminary
- The AI does NOT suggest running on full data until the phase's gates are passed

When `scale: small_then_full`:
- The AI starts with small-scale work
- After small-scale gates pass, the AI asks: "Ready to run on full data?"
- The full-scale run reuses parameters validated on the subsample
- The AI notes any differences between small and full results

### Carrying Learnings from Small to Full

The key outputs from a small-scale phase are:
1. **Validated parameters** (e.g., `n_latent=30, n_hidden=128`)
2. **Expected behavior** (e.g., "batch mixing should improve, cell types should stay separate")
3. **Known issues** (e.g., "dataset X has poor mixing even with integration")

These are captured in the checkpoint and referenced when the full-scale phase begins.

---

## 5. Decision Points

### Structure

```yaml
decision_points:
  - id: datasets_to_include
    description: Which datasets to carry forward
    prompt: "Based on exploration, which datasets will you integrate? Any to drop?"
    must_document: true
```

### How Decisions Are Handled

1. **The AI presents the decision point** with the prompt text and relevant context.
2. **The AI waits for the user's answer.** It may offer options or tradeoffs, but does not decide.
3. **The decision is logged** in `workflow_state.yaml`:

```yaml
decisions:
  explore.datasets_to_include:
    timestamp: 2026-03-07T15:00:00
    choice: "Include DC-VERSE (healthy) and PanCancer (tumor). Drop the mouse dataset."
    rationale: "Mouse data would require cross-species integration, not worth the complexity"
  integration_understanding.method_choice:
    timestamp: 2026-03-08T10:00:00
    choice: "scVI"
    rationale: "Deep learning handles large batch counts better than Harmony"
```

### Why `must_document: true`

Some decisions are critical for reproducibility. If `must_document: true`, the AI asks the user to explain their rationale, not just state the choice. This creates a decision log that makes the project reproducible by a third party.

---

## 6. Checkpoint/Waypoint Integration

### What Gets Saved

At the end of each phase (when gates are passed), the workflow state is updated:

```yaml
# 02_analysis/config/workflow_state.yaml
workflow: scrna-atlas
current_phase: full_integration
started: 2026-03-07T10:00:00
last_updated: 2026-03-08T11:00:00

phases:
  explore:
    status: completed
    started: 2026-03-07T10:00:00
    completed: 2026-03-07T16:00:00
    notebooks: [L01, L02, L03]
  integration_planning:
    status: completed
    started: 2026-03-07T16:30:00
    completed: 2026-03-08T09:00:00
    notebooks: [L04]
    artifacts:
      - path: 03_results/checkpoints/myeloid_subset.h5ad
        description: Subsetted myeloid cells for integration
  integration_understanding:
    status: completed
    started: 2026-03-08T09:30:00
    completed: 2026-03-08T10:30:00
    notebooks: [L05, L06]
    artifacts:
      - path: 03_results/checkpoints/subsample_integrated.h5ad
  full_integration:
    status: in_progress
    started: 2026-03-08T11:00:00
    notebooks: [L07]

gates:
  # ... (as shown in Section 3)

decisions:
  # ... (as shown in Section 5)
```

### Resume Protocol

When a new session starts:

1. AI reads `workflow_state.yaml`
2. Identifies `current_phase` and its status
3. If `in_progress`: resumes from where the user left off, loads relevant skills
4. If all phases `completed`: reports project status, asks what's next
5. AI reports: "You're in Phase 3 (Full Integration), started in notebook L07. Last session you trained the scVI model. Gates not yet checked."

### Going Back to a Previous Phase

Sometimes results in a later phase reveal problems in an earlier one (the v2 annotation failure is a perfect example).

```yaml
# User or AI can mark a phase as "needs_revisit"
phases:
  annotation_manual:
    status: needs_revisit
    reason: "Dotplot showed cDC1/cDC2 marker swap -- annotation is wrong"
    revisit_from: integration_understanding  # go back this far
```

The AI suggests which phase to revisit based on the nature of the problem:
- Bad annotations -> revisit `annotation_manual`
- Bad integration -> revisit `integration_understanding` (re-test parameters)
- Bad data -> revisit `explore` (check input data quality)

### Relationship to Handoff Agent

The handoff agent already captures session state in markdown files. The workflow state YAML is complementary:

| Concern | Handoff Agent | Workflow State |
|---------|---------------|----------------|
| What happened this session | Yes | No |
| Where the project is overall | Partially | Yes |
| What gates have passed | No | Yes |
| What decisions were made | Partially | Yes (structured) |
| Session-specific notes | Yes | No |

The handoff agent should read `workflow_state.yaml` when generating its handoff document, so the handoff includes phase context automatically.

---

## 7. Workflow Status Tracking

### Where State Lives

```
02_analysis/config/workflow_state.yaml    # machine-readable state
```

This is a single file. Not scattered across CLAUDE.md, not in multiple locations. One file, one truth.

### What CLAUDE.md Shows

CLAUDE.md has a "Current Status" section (it already does). The workflow system does not replace this -- it feeds it. When workflow state is updated, the CLAUDE.md status table can be regenerated:

```markdown
## Current Status

| Phase | Status | Key Outputs |
|-------|--------|-------------|
| 0. Explore | Completed | L01-L03 |
| 1. Integration Planning | Completed | L04, myeloid_subset.h5ad |
| 2. Integration Understanding | Completed | L05-L06, subsample_integrated.h5ad |
| 3. Full Integration | In Progress | L07 |
| 4. Manual Annotation | Not Started | -- |
| 5. Semi-Supervised Annotation | Not Started | -- |
| 6. Deep Analysis | Not Started | -- |
```

### Updating State

The AI updates `workflow_state.yaml` when:
- A gate is passed (with evidence)
- A decision is made (with rationale)
- A phase is completed
- A phase needs revisiting

The AI does NOT update state silently. It tells the user: "I'm marking the exploration phase as complete in workflow_state.yaml."

---

## 8. Example: Complete scRNA-Atlas Workflow

This is the full workflow from the user's description, expressed in the format above. The YAML in Section 1 IS this example -- it is not abstract, it is the actual `scrna-atlas.yaml` workflow definition.

### How a Session Looks in Practice

**Session 1: Starting fresh**

```
AI reads: workflow_state.yaml does not exist
AI reads: workflows/scrna-atlas.yaml
AI says: "This project uses the scrna-atlas workflow. You're starting at Phase 0:
         Independent Dataset Exploration. I'll load anndata and scanpy skills.
         Let's start by looking at your datasets."
```

**Session 3: Mid-exploration**

```
AI reads: workflow_state.yaml shows explore phase in_progress
AI says: "You're in Phase 0 (Exploration). Last session you explored the DC-VERSE
         dataset (L01). The PanCancer dataset hasn't been explored yet."
User: "Let's look at PanCancer"
AI: [works through exploration in L02]
User: "I think I've seen enough. What's next?"
AI: [presents verification gates for explore phase]
```

**Session 5: Gate check**

```
AI presents gates:
  Checklist: datasets_explored
    [x] Cell counts documented (L01 cell 12, L02 cell 8)
    [x] Annotations examined (L01 cell 15, L02 cell 11)
    [x] HVG overlap (L02 cell 20)
    [ ] Basic UMAP plots

  "Looks like we haven't generated UMAP plots for the raw datasets yet.
   Want to do that now, or skip it and move on?"

User: "Skip it, I'll look at UMAPs after integration"
AI: "OK, logging that as skipped. Moving to the human judgment gate..."
AI: "Can you summarize what each dataset brings to the atlas?"
User: "DC-VERSE has healthy DC reference with good cell type labels.
       PanCancer adds tumor-infiltrating DCs across 29 cancer types."
AI: "Great. I'm marking Phase 0 as complete and logging your summary.
     Phase 1 is Integration Planning. Ready?"
```

**Session 8: The v2 failure scenario**

```
AI reads: workflow_state.yaml shows annotation_manual in_progress
User: "The dotplot shows cDC1 and cDC2 markers are swapped"
AI: "That's a serious problem. Looking at the workflow, this could mean:
     - The integration mixed the populations (revisit Phase 3)
     - The clustering resolution is wrong (revisit within Phase 4)
     - The annotation logic has a bug (fix within Phase 4)

     What do you think is the root cause?"
User: "The integration looks fine. The annotation assigned cluster-level
       labels and some clusters had mixed populations."
AI: "So the fix is cell-level annotation, not re-integration.
     I'll mark annotation_manual as needs_revisit and we'll redo
     the annotation approach. Your integration checkpoint is still valid."
```

---

## Implementation Plan

This design needs three things to work:

1. **A workflow YAML file** (`workflows/scrna-atlas.yaml`) -- the contract
2. **A state file** (`02_analysis/config/workflow_state.yaml`) -- the tracker
3. **Instructions in CLAUDE.md** telling the AI to read both files at session start

That is it. No new scripts, no new infrastructure, no new MCP servers. The AI already knows how to read YAML and present checklists. The workflow file just gives it structure.

### What This Does NOT Include

- **Automatic skill loading based on phase.** The role system already handles skill activation. The workflow YAML lists skills as a suggestion, not as an automation trigger. The user can activate whatever role they want.
- **Automatic gate enforcement.** The AI presents gates and logs outcomes, but does not block the user. If the user wants to skip ahead, they can. The skip is logged.
- **Workflow editor UI.** The YAML is edited by hand or by the AI. No special tooling needed.
- **Multiple concurrent workflows.** One project, one workflow. If you need a different workflow, swap the YAML file.

### What This Replaces

- The ad-hoc "Current Status" table in CLAUDE.md that gets out of date
- The informal "where were we?" conversation at the start of each session
- The undocumented decisions that led to the v2 annotation failure

---

*End of design document.*
