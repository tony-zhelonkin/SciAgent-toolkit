# Skills Format v3: Design Specification

> **Date:** 2026-03-07 (adopted 2026-04-14)
> **Status:** Adopted — all skills migrated to directory format
> **Scope:** SciAgent-toolkit skill format evolution
> **Predecessor:** Flat `.md` skill files (v2), informed by claude-scientific-skills and life-sciences reference repos
> **Starter template:** `skills/_TEMPLATE/` — copy this when authoring new skills
> **Reference implementation:** `skills/skill-creator/` — canonical Rich-tier example
> **Activator support:** `scripts/activate-role.sh` resolves `skills/<name>/SKILL.md` (dir) with `skills/<name>.md` (flat) fallback

---

## Table of Contents

1. [Design Principles](#1-design-principles)
2. [Skill Directory Format](#2-skill-directory-format)
3. [SKILL.md Specification](#3-skillmd-specification)
4. [Skill Tiers](#4-skill-tiers)
5. [Skill Categories](#5-skill-categories)
6. [Practice Skills](#6-practice-skills-new-category)
7. [Workflow Skills](#7-workflow-skills-new-category)
8. [Provider Agnosticism](#8-provider-agnosticism)
9. [Migration Strategy](#9-migration-strategy)
10. [Worked Example: scvi-basic](#10-worked-example-scvi-basic)

---

## 1. Design Principles

These principles govern all format decisions. When in doubt, refer back here.

1. **Pedagogical first.** Skills exist to help the user *understand*, not just execute. Every skill must answer "why" before "how."
2. **Verification is mandatory, not optional.** Every skill that produces output must include verification steps. "Trust but verify" is insufficient; "verify then proceed" is the standard.
3. **Progressive disclosure.** Start with the simplest correct approach. Reveal complexity only when the user needs it. Never front-load advanced topics.
4. **Provider-agnostic by default.** SKILL.md content must work with any AI assistant (Claude Code, Gemini CLI, Codex, or a human reading the file). Provider-specific instructions go in clearly marked sections or separate files.
5. **Composable over monolithic.** Skills should declare what they pair with and what they replace. A user should never wonder "do I use skill A or skill B?"
6. **Contraindications are as important as indications.** "When NOT to use" prevents more errors than "when to use" enables.
7. **Executable over aspirational.** Code in `scripts/` must run. Checklists in `checks/` must be actionable. Templates in `assets/` must be fillable.

---

## 2. Skill Directory Format

### 2.1 Directory Structure

Each skill becomes a directory. The directory name matches the skill identifier (kebab-case).

```
skill-name/
    SKILL.md              # REQUIRED: Entry point, metadata, navigation
    references/           # OPTIONAL: Deep-dive topic guides (*.md)
    scripts/              # OPTIONAL: Executable Python/Bash/R utilities
    assets/               # OPTIONAL: Templates, configs, example data schemas
    checks/               # OPTIONAL: Verification scripts and checklists
```

### 2.2 File Naming Conventions

| Directory | Naming Pattern | Example |
|-----------|---------------|---------|
| `references/` | `topic_name.md` (snake_case) | `differential_expression.md`, `batch_correction.md` |
| `scripts/` | `verb_noun.py` or `verb_noun.sh` (snake_case) | `validate_adata.py`, `train_model.py` |
| `assets/` | Descriptive, format-specific | `config_template.yaml`, `marker_panel.csv` |
| `checks/` | `check_noun.md` or `check_noun.py` | `check_integration_quality.md`, `check_counts_layer.py` |

### 2.3 SKILL.md as Entry Point

SKILL.md is the *only* file an AI assistant reads by default. It must be self-contained enough to handle common use cases. References, scripts, and checks are loaded on demand when the assistant needs deeper guidance.

The SKILL.md file should:
- Be readable and useful on its own (Simple tier skills have nothing else)
- Reference deeper files by relative path when needed: "See `references/batch_correction.md` for advanced batch handling"
- Never exceed ~400 lines (if it does, factor content into references/)

---

## 3. SKILL.md Specification

### 3.1 Frontmatter (YAML)

```yaml
---
name: skill-identifier                  # REQUIRED: kebab-case, matches directory name
description: |                          # REQUIRED: 1-3 sentence trigger description
  When to activate this skill and what it covers.
  Mention key trigger words for AI assistant matching.

# Taxonomy
category: foundation | integration | annotation | analysis | workflow | practice
tier: simple | standard | rich
tags: [single-cell, python, scverse]    # Free-form tags for discovery

# Relationships
complementary-skills:                   # Skills that pair well with this one
  - anndata                             # "Use anndata for format questions"
  - scanpy                              # "Use scanpy for downstream clustering"
contraindications:                      # When NOT to use this skill
  - "Do not use for multi-modal data (CITE-seq, multiome) -- use scvi-multivi"
  - "Do not use if you have fewer than 1000 cells -- classical methods suffice"

# Optional metadata
version: 1.0.0                          # Skill version (semver)
upstream-docs: https://docs.example.org # Official tool documentation
license: MIT
metadata:
  skill-author: SciAgent-toolkit
  last-reviewed: 2026-03-07
---
```

**Required fields:** `name`, `description`, `category`, `tier`
**Recommended fields:** `tags`, `complementary-skills`, `contraindications`
**Optional fields:** `version`, `upstream-docs`, `license`, `metadata`

### 3.2 Body Sections

Every SKILL.md follows this section order. Sections marked REQUIRED must be present; others are included based on tier.

#### Section 1: Title and Overview (REQUIRED)

```markdown
# Skill Title

## Overview

What problem does this skill solve? Who needs it and why?
One paragraph, max 4 sentences.
```

#### Section 2: Decision Tree (REQUIRED for Standard and Rich tiers)

When to use *this* skill versus alternatives. Uses a text-based decision tree.

```markdown
## Decision Tree

```
Have scRNA-seq data needing batch correction?
|
+-- Have cell type labels for some cells?
|   +-- Yes --> scANVI (scvi-scanvi skill)
|   +-- No  --> scVI (this skill)
|
+-- Single dataset, no batches?
    +-- Use scanpy directly (scanpy skill)
```
```

The decision tree is the single most important routing mechanism. It prevents the user from using the wrong skill entirely. Rules for decision trees:

- Start from the user's *situation*, not the tool's features
- Every leaf node names a specific skill
- Include a "you might not need this at all" branch where appropriate

#### Section 3: Quick Start (REQUIRED)

Minimal working example. Must be copy-pasteable and correct. Should take less than 10 lines of code.

```markdown
## Quick Start

```python
# Minimal working example - 5 lines to a result
import scvi
adata.layers["counts"] = adata.X.copy()
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")
model = scvi.model.SCVI(adata)
model.train()
adata.obsm["X_scVI"] = model.get_latent_representation()
```

**Verify it worked:**
```python
assert "X_scVI" in adata.obsm, "Latent representation not found"
assert adata.obsm["X_scVI"].shape == (adata.n_obs, 10), f"Unexpected shape: {adata.obsm['X_scVI'].shape}"
print(f"scVI latent space: {adata.obsm['X_scVI'].shape}")
```
```

Quick Start MUST include a verification step. This is non-negotiable. The verification block is labelled "Verify it worked:" and contains assertions that would catch the most common failure modes.

#### Section 4: Progressive Depth (REQUIRED for Standard and Rich tiers)

Three levels of depth, clearly delineated. The user reads only as far as they need.

```markdown
## Basic Usage

Core functionality that covers 80% of use cases.
Parameters with sensible defaults. Simple code examples.

## Intermediate Usage

Parameter tuning, multiple approaches, conditional logic.
When defaults are not enough. Cross-references to references/.

## Advanced Usage

Edge cases, custom architectures, performance optimization.
Links to references/ for full treatment.
For expert guidance, see `references/advanced_training.md`.
```

The boundary between levels:
- **Basic:** "I want the default, correct approach." No parameter decisions required.
- **Intermediate:** "The default didn't work well enough, or I have a specific requirement." Requires understanding parameter trade-offs.
- **Advanced:** "I need to customize the internals or handle an unusual situation." Requires understanding the tool's architecture.

#### Section 5: Verification Checklist (REQUIRED)

What to check after using this skill. Each item is a concrete, executable check -- not a vague instruction.

```markdown
## Verification Checklist

After completing this workflow, verify:

- [ ] **Data integrity**: `adata.layers["counts"]` contains integer counts
- [ ] **No NaN in latent space**: `np.isnan(adata.obsm["X_scVI"]).sum() == 0`
- [ ] **Batch mixing**: UMAP colored by batch shows interleaving, not isolated islands
- [ ] **Biological signal preserved**: Known marker genes still separate expected populations on UMAP
- [ ] **Training converged**: Loss curve plateaued (check `model.history["elbo_train"]`)

For automated verification, run: `python checks/check_integration_quality.py`
```

Each checklist item must be:
- **Observable**: You can look at something to confirm it
- **Falsifiable**: There is a concrete failure condition
- **Executable**: There is code or a visual check, not just a hope

#### Section 6: Common Pitfalls (REQUIRED)

Mistakes that are frequently made and costly. Each pitfall uses a consistent three-part format: Symptom, Cause, Fix.

```markdown
## Common Pitfalls

### Pitfall: Passing normalized data to scVI
**Symptom:** Model trains but integration is poor; loss is unusually low.
**Cause:** scVI expects raw integer counts. Normalized/log-transformed data violates the model's assumptions.
**Fix:** Always use `layer="counts"` pointing to raw counts. Verify with:
```python
assert adata.layers["counts"].data.dtype in [np.int32, np.int64, np.float32]
# If float, check they're actually integers:
assert np.allclose(adata.layers["counts"].data, adata.layers["counts"].data.astype(int))
```
```

The Symptom/Cause/Fix structure matters because users usually encounter the symptom first and need to work backwards. Never describe just the cause without explaining what it looks like in practice.

#### Section 7: Complementary Skills (RECOMMENDED)

Explicit cross-references showing how this skill connects to others. Uses a table with a "Relationship" column that clarifies directionality.

```markdown
## Complementary Skills

| When you need... | Use skill | Relationship |
|-----------------|-----------|-------------|
| Data format operations | `anndata` | Prerequisite |
| Downstream clustering/viz | `scanpy` | Next step |
| Semi-supervised annotation | `scvi-scanvi` | Extension |
| Query-to-reference mapping | `scvi-scarches` | Extension |
| Verify QC before integration | `single-cell-rna-qc` | Prerequisite |
```

Valid relationship values: `Prerequisite`, `Next step`, `Extension`, `Alternative`, `Foundation`.

#### Section 8: Resources (OPTIONAL)

Links to official docs, papers, tutorials. Keep this short (3-5 links max).

---

## 4. Skill Tiers

### 4.1 Tier Definitions

| Tier | Contents | When to Use | Example |
|------|----------|-------------|---------|
| **Simple** | `SKILL.md` only | Straightforward tools with one primary use case, few parameters, minimal decision-making | `louper-seurat-conversion`, `starsolo-spliced-unspliced` |
| **Standard** | `SKILL.md` + `references/` | Tools with multiple use cases requiring deeper guides, but no executable automation | `anndata`, `scanpy`, `scvi-basic` |
| **Rich** | `SKILL.md` + `references/` + `scripts/` + optional `assets/`, `checks/` | Complex workflows with executable pipelines, templates, and automated verification | `single-cell-rna-qc`, workflow skills |

### 4.2 Tier Selection Criteria

Use **Simple** when:
- The tool does one thing
- No decision tree is needed (or a trivial one)
- Quick Start covers 90%+ of use cases
- The entire skill fits comfortably in <200 lines

Use **Standard** when:
- Multiple distinct use cases exist
- Decision trees help users pick the right approach
- Deep-dive guides would exceed 400 lines in SKILL.md
- No executable scripts are needed (all guidance is textual)

Use **Rich** when:
- Executable scripts save significant user effort
- Automated verification catches common errors
- Templates or configuration files are needed
- The skill is a multi-step workflow with validation gates

### 4.3 Tier Progression

Skills can be promoted from Simple to Standard to Rich over time without breaking anything. The directory structure is additive:
1. Start Simple: just `SKILL.md`
2. Add `references/` to become Standard
3. Add `scripts/` and/or `checks/` to become Rich

Update the `tier` field in SKILL.md frontmatter to match actual contents. The tier is descriptive (what exists), not prescriptive (what must exist).

### 4.4 What Each Tier Requires in SKILL.md

| Section | Simple | Standard | Rich |
|---------|--------|----------|------|
| Frontmatter | REQUIRED | REQUIRED | REQUIRED |
| Overview | REQUIRED | REQUIRED | REQUIRED |
| Decision Tree | Optional | REQUIRED | REQUIRED |
| Quick Start | REQUIRED | REQUIRED | REQUIRED |
| Progressive Depth | Optional (single level OK) | REQUIRED (3 levels) | REQUIRED (3 levels) |
| Verification Checklist | REQUIRED | REQUIRED | REQUIRED |
| Common Pitfalls | REQUIRED | REQUIRED | REQUIRED |
| Complementary Skills | Optional | RECOMMENDED | REQUIRED |

---

## 5. Skill Categories

### 5.1 Category Definitions

| Category | Purpose | What It Contains | Typical Tier | Examples |
|----------|---------|-----------------|-------------|---------|
| **Foundation** | Data structures, QC, format conversion | How to read, write, subset, validate data | Simple-Standard | `anndata`, `single-cell-rna-qc`, `anndatar-seurat-scanpy-conversion` |
| **Integration** | Batch correction, cross-dataset harmonization | How to combine datasets while preserving biology | Standard-Rich | `scvi-basic`, `scvi-scanvi`, `scglue-unpaired-multiomics-integration` |
| **Annotation** | Cell typing, label transfer, reference mapping | How to assign biological identity to cells | Standard-Rich | `cellxgene-census-annotation`, `scvi-scarches-reference-mapping` |
| **Analysis** | DE, GSEA, trajectory, GRN inference, scoring | How to extract biological insight from annotated data | Standard-Rich | `scanpy`, `scenic-grn-inference`, `rna-velocity-trajectory` |
| **Workflow** | Multi-step procedures composing other skills | How to connect steps into a coherent pipeline with gates | Rich (always) | `scrna-eda-workflow`, `integration-workflow`, `annotation-workflow` |
| **Practice** | Coding habits, verification patterns, mental models | How to *think about* analysis, not how to *run* a tool | Simple-Standard | `checkpoint-verify`, `small-sample-test`, `batch-aware-analysis` |

### 5.2 Category in Frontmatter

```yaml
category: foundation   # One of: foundation, integration, annotation, analysis, workflow, practice
```

### 5.3 Category-Based Discovery

Roles can optionally filter skills by category. This enables role definitions like:

```yaml
# In roles/learning.yaml (future enhancement)
skills:
  categories:
    - foundation
    - practice
  additional:
    - scvi-basic     # One integration skill for learning
```

This is a future enhancement to `activate-role.sh`. Current activation continues to work with explicit skill name lists.

---

## 6. Practice Skills (NEW Category)

Practice skills encode *how to work*, not *how to use a tool*. They are the codified habits of a careful computational biologist.

### 6.1 Purpose

Practice skills solve a different problem than tool-reference skills:

| Tool skill answers | Practice skill answers |
|---|---|
| "How do I use scVI?" | "How do I avoid wasting a day on a doomed run?" |
| "What parameters does scanpy.tl.leiden take?" | "How do I verify that my clusters mean something?" |
| "How do I run DE?" | "How do I know my DE results reflect biology, not batch?" |

Practice skills are shorter, more opinionated, and less code-heavy than tool skills. They encode lessons learned from real analysis failures.

### 6.2 Proposed Practice Skills

#### `checkpoint-verify`

**Problem:** Expensive computations (scVI training, large DE) fail silently or produce subtly wrong results. Without checkpoints, the user loses hours of work; without verification, they build on a broken foundation.

**Key content:**
- When and how to save intermediate results (`adata.write_h5ad()`, `saveRDS()`, `model.save()`)
- Verification patterns: shape checks, value range checks, distribution sanity
- The "checkpoint contract": what must be true before saving, and what to check after loading
- Resume patterns: loading a checkpoint and confirming its integrity before proceeding

**Tier:** Standard (SKILL.md + `references/checkpoint_patterns.md` with language-specific examples)

#### `small-sample-test`

**Problem:** Running a 1M-cell pipeline end-to-end before knowing if the code works wastes GPU-hours and calendar time.

**Key content:**
- How to create representative subsamples (stratified by batch, cell type)
- What "representative" means for different analyses (QC needs all batches, DE needs both conditions)
- Size recommendations: 5K cells for code validation, 50K for parameter tuning, full for production
- How to validate that subsample results predict full-scale behavior
- The one-line subsample: `adata_sub = sc.pp.subsample(adata, n_obs=5000, copy=True)`

**Tier:** Simple (SKILL.md only -- the practice is simple even if the motivation is nuanced)

#### `declarative-exploration`

**Problem:** Jumping to complex analysis before understanding the data leads to misinterpretation. EDA is not a luxury; it is a prerequisite.

**Key content:**
- The "look before you compute" principle
- Standard EDA checklist: shape, sparsity, batch sizes, label distributions, marker expression
- How to explore each dataset independently before integration
- Visual verification: UMAP colored by batch, dotplot of known markers, QC metric distributions
- Connection to the `scrna-eda-workflow` workflow skill

**Tier:** Standard (SKILL.md + `references/eda_checklists.md`)

#### `separation-of-concerns`

**Problem:** Notebooks that mix data loading, computation, and visualization are hard to debug, slow to iterate on, and impossible to checkpoint.

**Key content:**
- The three-phase pattern: compute -> cache -> visualize
- Why functions should never both compute AND plot
- How to structure notebook cells: load cell, compute cell, checkpoint cell, viz cell
- Connection to RNAseq-toolkit's design principles (separation of processing and visualization)

**Tier:** Simple (SKILL.md only)

#### `batch-aware-analysis`

**Problem:** Comparing gene expression across batches without accounting for batch effects finds batch artifacts, not biology. This is the single most common analytical error in multi-sample scRNA-seq.

**Key content:**
- Why naive cross-batch count comparisons detect batch effects, not biology
- The three correct approaches: pseudobulk DE with batch covariate, scVI native DE (`model.differential_expression()`), MAST with batch random effect
- When batch effects are and are NOT a concern (within-batch comparisons are safe)
- How to verify that DE results are not driven by a single batch/sample
- The pseudobulk sanity check: aggregate to sample level, check if effect survives

**Tier:** Standard (SKILL.md + `references/batch_de_methods.md` with code for each approach)

### 6.3 Practice Skill SKILL.md Template

Practice skills follow a modified section structure. They replace the tool-oriented sections with problem-oriented ones:

```yaml
---
name: practice-name
description: ...
category: practice
tier: simple
tags: [methodology, verification]
complementary-skills: []
contraindications: []
---
```

```markdown
# Practice Name

## The Problem

What goes wrong without this practice? A concrete example from real analysis.
Not abstract -- describe an actual failure mode the user might experience.

## The Principle

One sentence. The core idea. Example: "Always verify intermediate results
before using them as input to the next step."

## The Pattern

Step-by-step procedure. Code examples where applicable.
This is the "how to do it" section.

## When to Apply

Situations where this practice is essential.
Tied to specific analysis stages or data types.

## When to Skip

Situations where this practice adds overhead without value.
Be honest -- not every practice applies everywhere. Example:
"Skip subsampling for datasets under 10K cells."

## Signs You Skipped It

Observable symptoms that indicate this practice was not followed.
These help the user diagnose problems retroactively.
Example: "If your DE results contain 5000+ significant genes,
you may have batch confounding."
```

---

## 7. Workflow Skills (NEW Category)

Workflow skills are multi-step procedures that compose other skills into a coherent pipeline. They are the "recipes" that connect individual tool skills into analysis workflows.

### 7.1 Purpose

A workflow skill answers: "I have data type X and want result Y. What sequence of steps, using which tools, with what verification gates between them?"

Workflow skills do NOT duplicate the content of tool skills. They reference them. A workflow skill says "run QC using the `single-cell-rna-qc` skill" -- it does not re-explain how QC works.

### 7.2 Key Differentiator: Validation Gates

The defining feature of workflow skills is explicit gates between steps. A gate is a verification check that MUST pass before proceeding to the next step. This prevents cascading errors -- the most expensive category of analysis failure.

```
Step 1: Load and QC data
    | GATE: cell count > 500, median genes > 200, no empty batches
Step 2: Integrate with scVI
    | GATE: UMAP shows batch mixing, training loss converged
Step 3: Cluster and annotate
    | GATE: dotplot confirms marker-label correspondence
Step 4: Differential expression
    | GATE: results are batch-aware, not driven by single sample
```

Each gate specifies:
1. What to check (concrete assertion or visual inspection)
2. What passing looks like
3. What failing looks like
4. What to do if it fails (which step to revisit, what parameter to change)

### 7.3 Proposed Workflow Skills

#### `scrna-eda-workflow`

**Scope:** From raw .h5ad to understanding what you have, before any integration.

**Steps and gates:**

| Step | Skill Used | Gate |
|------|-----------|------|
| 1. Load and inspect shape/metadata | `anndata` | Data loads, shape is reasonable, expected metadata columns exist |
| 2. QC metrics and filtering | `single-cell-rna-qc` | Cell/gene counts reasonable, QC metrics plotted, filtering thresholds justified |
| 3. Normalize, HVG, PCA, UMAP | `scanpy` | UMAP renders, clusters visible, no single-cell outlier islands |
| 4. Batch structure assessment | `declarative-exploration` | Batch sizes documented, batch effect severity assessed visually |
| 5. Marker gene survey | `scanpy` | Key markers detected in dataset, expression patterns plausible |

**Output:** Understanding of what integration strategy is needed (or if none is needed).

#### `integration-workflow`

**Scope:** From QC'd datasets to an integrated atlas with batch-corrected latent space.

**Steps and gates:**

| Step | Skill Used | Gate |
|------|-----------|------|
| 1. Assess integration need | `batch-aware-analysis` | Decision documented: integrate vs. analyse separately |
| 2. Prepare AnnData for scVI | `anndata` | Counts layer verified, batch key set, HVGs selected |
| 3. Subsample test run | `small-sample-test` | 5K-cell test completes, UMAP shows mixing |
| 4. Full scVI training | `scvi-basic` | Training converged (loss plateau), no NaN in latent |
| 5. Evaluate integration quality | `checkpoint-verify` | Batch mixing confirmed, biology preserved |
| 6. Downstream: neighbors, UMAP, Leiden | `scanpy` | Clusters are biologically interpretable |

**Output:** Integrated atlas saved as checkpoint, ready for annotation.

#### `annotation-workflow`

**Scope:** From integrated atlas to validated cell type annotations.

**Steps and gates:**

| Step | Skill Used | Gate |
|------|-----------|------|
| 1. Survey marker expression | `scanpy` | Key markers detected, expression localised to UMAP regions |
| 2. Choose annotation strategy | Decision tree in SKILL.md | Strategy justified based on available labels and reference data |
| 3. Execute annotation | `scvi-scanvi` or `cellxgene-census-annotation` | Labels assigned to all cells |
| 4. Validate with dotplot | Verification checklist | Marker-label correspondence confirmed, no swaps |
| 5. Handle ambiguous cells | Manual review | No massive catch-all category (>30% of cells in one label) |

**Output:** Annotated atlas with validated labels, saved as checkpoint.

### 7.4 Workflow Skill SKILL.md Template

```yaml
---
name: workflow-name
description: ...
category: workflow
tier: rich
tags: [workflow, multi-step]
complementary-skills: []    # Usually empty -- the workflow itself references skills
contraindications: []
---
```

```markdown
# Workflow Name

## Overview

What this workflow achieves. Input -> Output. When to use it.

## Prerequisites

- [ ] Skills needed: `skill-a`, `skill-b`, `skill-c`
- [ ] Data requirements: format, minimum size, required metadata
- [ ] Environment: packages, GPU availability, memory estimate

## Steps

### Step 1: [Step Name]

**Skill:** `skill-a`
**Input:** [what you start with]
**Output:** [what this step produces]

[Concise instructions, referencing skill-a's Quick Start or specific section.
Do NOT re-explain what skill-a already covers.]

#### Gate 1

```python
# These assertions must pass before proceeding to Step 2
assert condition_1, "Failure message explaining what went wrong"
assert condition_2, "Failure message explaining what went wrong"
print("Gate 1 passed -- proceed to Step 2")
```

**If gate fails:**
- condition_1 failure: [what to do -- e.g., revisit Step 1 with different parameters]
- condition_2 failure: [what to do]

### Step 2: [Step Name]

[... pattern repeats ...]

## Checkpoints

| After Step | Checkpoint File | Contents | How to Resume |
|-----------|----------------|----------|---------------|
| Step 2 | `checkpoints/qc_filtered.h5ad` | QC'd AnnData | `adata = sc.read_h5ad(...)` |
| Step 4 | `checkpoints/integrated.h5ad` + `scvi_model/` | Integrated atlas + model | `adata = sc.read_h5ad(...); model = scvi.model.SCVI.load(...)` |
| Step 6 | `checkpoints/annotated.h5ad` | Final annotated object | `adata = sc.read_h5ad(...)` |

## Resume Protocol

To resume from a checkpoint:
```python
import scanpy as sc
adata = sc.read_h5ad("checkpoints/integrated.h5ad")

# Verify checkpoint integrity before continuing
assert "X_scVI" in adata.obsm, "Missing latent representation"
assert adata.obs["batch"].nunique() > 1, "Batch info missing"
print(f"Resuming from integrated checkpoint: {adata.shape}")
# Continue from Step 5
```
```

---

## 8. Provider Agnosticism

### 8.1 What is Universal (in SKILL.md)

Everything in SKILL.md is provider-agnostic. It uses standard Markdown, standard code blocks, and plain-text decision trees. Any AI assistant -- or a human reader -- can use it.

Specifically:
- Decision trees use text-based ASCII art, not provider-specific syntax
- Code examples use standard Python/R/Bash with no provider tool invocations
- Verification checklists use standard markdown checkboxes
- References use relative file paths within the skill directory
- No provider-specific YAML keys in frontmatter (the `allowed-tools` key from claude-scientific-skills is NOT adopted -- it is provider-specific and belongs elsewhere)

### 8.2 What is Provider-Specific

Provider-specific content is rare and isolated. When needed, it goes in one of two places:

**Option A: A clearly labelled section at the end of SKILL.md** (for minor notes):

```markdown
## Provider Notes

### Claude Code
- Skills are activated via `activate-role.sh` which creates symlinks in `.claude/skills/`
- References are loaded on-demand when Claude reads relative paths from SKILL.md

### Gemini CLI
- Skills are mapped via `.gemini/settings.json` by `switch-mcp-profile.sh`
- References may need explicit inclusion in context window
```

**Option B: A separate file** (for substantial provider-specific content):

```
skill-name/
    SKILL.md
    provider_claude.md      # Claude Code-specific activation/usage
    provider_gemini.md      # Gemini-specific configuration
```

This separation is only needed for skills with significant provider differences (rare).

### 8.3 Scripts and Environment

Scripts in `scripts/` are standard, portable executables. They:
- Use standard library imports (no provider SDKs)
- Accept command-line arguments (`argparse` for Python, `getopts` for Bash)
- Print to stdout/stderr with structured output
- Return standard exit codes (0 for success, 1 for validation failure, 2 for usage error)
- Include a shebang line (`#!/usr/bin/env python3`)
- Document their own usage in a docstring or `--help` output

Scripts do NOT:
- Import provider-specific libraries (`anthropic`, `google.generativeai`)
- Assume a specific working directory (use paths relative to the script's own location, or accept paths as arguments)
- Require interactive input (must be fully scriptable)
- Write to hardcoded absolute paths

### 8.4 Activation Mechanism Adaptation

The `activate-role.sh` script needs a minor update for directory-based skills. The key insight is that the result in `.claude/skills/` remains the same -- a `.md` file named after the skill:

**v2 (current):** Symlinks `toolkit/skills/skill-name.md` to `.claude/skills/skill-name.md`
**v3 (new):** Symlinks `toolkit/skills/skill-name/SKILL.md` to `.claude/skills/skill-name.md`

The symlink target changes, but the consumer (Claude Code, Gemini) sees the same thing: a file called `skill-name.md` in the skills directory. When the AI assistant follows relative paths mentioned in SKILL.md (e.g., `references/topic.md`), the symlink resolves to the real directory in the toolkit, making all subdirectories accessible.

For Gemini CLI, the `switch-mcp-profile.sh` script applies the same logic when writing `.gemini/settings.json`.

---

## 9. Migration Strategy

### 9.1 Principles

- **Non-breaking:** Current flat `.md` skills continue to work throughout migration
- **Incremental:** Skills are migrated one at a time, on any schedule
- **Backward-compatible:** `activate-role.sh` handles both formats simultaneously
- **No content loss:** Existing skill content is preserved and reorganized, not rewritten

### 9.2 Phase 1: Infrastructure (do first, once)

Update `activate-role.sh` to handle both formats transparently:

```bash
# In the skill symlinking loop:
for skill in "${skills[@]}"; do
    if [ -d "$TOOLKIT_DIR/skills/$skill" ] && [ -f "$TOOLKIT_DIR/skills/$skill/SKILL.md" ]; then
        # v3 directory format: symlink SKILL.md as skill-name.md
        ln -sf "$TOOLKIT_DIR/skills/$skill/SKILL.md" "$PROJECT_DIR/.claude/skills/$skill.md"
    elif [ -f "$TOOLKIT_DIR/skills/$skill.md" ]; then
        # v2 flat format: symlink the .md file directly
        ln -sf "$TOOLKIT_DIR/skills/$skill.md" "$PROJECT_DIR/.claude/skills/$skill.md"
    else
        log_warn "Skill not found: $skill"
    fi
done
```

This is the only blocking prerequisite. Once deployed, skills can be migrated at any pace without coordination.

### 9.3 Phase 2: Migrate Existing Skills (incremental)

For each skill to migrate:

1. Create directory: `skills/skill-name/`
2. Move content: `mv skills/skill-name.md skills/skill-name/SKILL.md`
3. Add missing frontmatter fields: `category`, `tier`, `complementary-skills`, `contraindications`
4. Add `## Verification Checklist` section (if missing)
5. Add `## Decision Tree` section (if Standard or Rich tier and missing)
6. Ensure Quick Start includes a "Verify it worked" block
7. Restructure Common Pitfalls to use Symptom/Cause/Fix format (if not already)
8. If the file exceeds ~400 lines, factor deep-dive sections into `references/`
9. Update `tier` in frontmatter to match actual directory contents
10. Test: `activate-role.sh <role> --project-dir <project>` and verify symlink resolves

**Migration priority** (highest value, migrate first):

| Priority | Skill | Reason |
|----------|-------|--------|
| 1 | `scvi-basic` | Most complex, most used, benefits most from references/ and checks/ |
| 2 | `scanpy` | Long file (326 lines), benefits from factoring into references/ |
| 3 | `anndata` | Foundational, referenced by everything |
| 4 | `single-cell-rna-qc` | Already has a Rich equivalent in life-sciences repo to draw from |
| 5 | `scvi-scanvi` | Closely related to scvi-basic, migrate together |
| 6-N | Remaining skills | In order of usage frequency within active roles |

### 9.4 Phase 3: Add New Skills (after migration infrastructure exists)

Create new practice and workflow skills directly in v3 format:

| Priority | Skill | Category | Reason |
|----------|-------|----------|--------|
| 1 | `checkpoint-verify` | practice | Highest immediate value -- prevents lost work |
| 2 | `small-sample-test` | practice | Prevents wasted compute time |
| 3 | `batch-aware-analysis` | practice | Prevents the most common analytical error |
| 4 | `scrna-eda-workflow` | workflow | Structures the exploration phase |
| 5 | `separation-of-concerns` | practice | Quick win, simple tier |
| 6 | `declarative-exploration` | practice | Supports EDA workflow |
| 7 | `integration-workflow` | workflow | Structures the integration phase |
| 8 | `annotation-workflow` | workflow | Structures the annotation phase |

### 9.5 Phase 4: Cleanup (not urgent)

Once all skills are migrated to v3 format:
- Remove the v2 flat-format fallback from `activate-role.sh` (optional -- the dual-format support has zero maintenance cost)
- Update role YAML files if any skill names changed (they should not -- only the internal structure changes)
- Update `skills/README.md` to document the v3 format

---

## 10. Worked Example: scvi-basic

This shows the complete transformation of the current `scvi-basic.md` (156 lines, flat file) into v3 format.

### 10.1 Directory Structure

```
scvi-basic/
    SKILL.md                              # Reorganized from scvi-basic.md (~300 lines)
    references/
        training_parameters.md            # Deep dive: n_latent, n_layers, gene_likelihood, dispersion
        differential_expression.md        # Bayesian DE with scVI (moved from SKILL.md)
        library_size_handling.md          # Observed vs modeled library size (moved from SKILL.md)
    scripts/
        validate_before_scvi.py           # Pre-training data validation
    checks/
        check_integration_quality.md      # Visual + quantitative post-training checklist
```

### 10.2 SKILL.md

```yaml
---
name: scvi-basic
description: |
  scVI unsupervised integration and batch correction for scRNA-seq.
  Use when integrating multiple batches/samples without cell type labels.
  For semi-supervised annotation transfer, use scvi-scanvi instead.
category: integration
tier: rich
tags: [single-cell, scvi-tools, batch-correction, integration, deep-learning]
complementary-skills:
  - anndata
  - scanpy
  - scvi-scanvi
  - scvi-framework
  - single-cell-rna-qc
contraindications:
  - "Do not use if you have cell type labels for most cells -- use scvi-scanvi for better integration"
  - "Do not use for single-sample data with no batch structure -- use scanpy PCA directly"
  - "Do not use for multi-modal data (CITE-seq) -- use scvi-totalvi"
  - "Requires GPU for practical runtimes on >100K cells"
version: 1.0.0
upstream-docs: https://docs.scvi-tools.org/en/stable/user_guide/models/scvi.html
---
```

```markdown
# scVI: Unsupervised Integration and Batch Correction

## Overview

scVI learns a shared latent representation across batches, removing technical variation
while preserving biological signal. It is the standard first step for multi-sample
scRNA-seq integration when you do NOT have cell type labels. The latent space replaces
PCA for all downstream analysis (neighbors, UMAP, clustering).

**Foundation:** See `scvi-framework` skill for installation, data prep, and core patterns
shared across all scvi-tools models.

## Decision Tree

```
Need to integrate multi-batch scRNA-seq?
|
+-- Have cell type labels for >50% of cells?
|   +-- Yes --> scvi-scanvi (semi-supervised, better integration)
|   +-- No  --> scVI (this skill)
|
+-- Single batch, no integration needed?
|   +-- Yes --> scanpy PCA (no deep learning needed)
|
+-- Multi-modal data?
|   +-- CITE-seq (RNA + protein) --> scvi-totalvi
|   +-- Multiome (RNA + ATAC)   --> scvi-multivi
|   +-- scATAC-seq only         --> scvi-peakvi
|
+-- Have a pre-trained reference model?
    +-- Yes --> scvi-scarches (query-to-reference mapping)
```

## Quick Start

```python
import scvi
import scanpy as sc

scvi.settings.seed = 42

# PREREQUISITE: adata must have raw integer counts
adata.layers["counts"] = adata.X.copy()
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3",
                            layer="counts", batch_key="batch", subset=True)

# Train
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")
model = scvi.model.SCVI(adata, n_latent=30, n_layers=2, gene_likelihood="nb")
model.train()

# Extract and use latent space
adata.obsm["X_scVI"] = model.get_latent_representation()
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata)
```

**Verify it worked:**
```python
import numpy as np
assert "X_scVI" in adata.obsm, "Latent representation missing"
assert not np.isnan(adata.obsm["X_scVI"]).any(), "NaN values in latent space"
assert adata.obsm["X_scVI"].shape == (adata.n_obs, 30), \
    f"Unexpected shape: {adata.obsm['X_scVI'].shape}"
# Visual: batches should interleave on UMAP
sc.pl.umap(adata, color="batch", title="Batch mixing check")
```

## Basic Usage

### Key Model Parameters

```python
model = scvi.model.SCVI(
    adata,
    n_hidden=128,          # Hidden layer size (default 128)
    n_latent=30,           # Latent dimensions (10-50 typical; start with 30)
    n_layers=2,            # Encoder/decoder depth (1-3)
    dropout_rate=0.1,      # Regularization
    gene_likelihood="nb",  # "nb" often better than default "zinb" for integration
    dispersion="gene",     # "gene", "gene-batch", or "gene-label"
)
```

**Gene likelihood tip:** Try `"nb"` first. Switch to `"zinb"` only if you see zero-inflation artifacts.

### Outputs

```python
# Latent representation (posterior mean) -- use this for neighbors/UMAP
latent = model.get_latent_representation()

# Denoised, batch-corrected expression (for visualization, not DE)
norm = model.get_normalized_expression(library_size=1e4)

# Counterfactual: expression as if all cells were in one batch
norm_batch1 = model.get_normalized_expression(transform_batch="batch1")
```

### Transitioning to scANVI

Always initialize scANVI from a trained scVI model:

```python
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    model,
    labels_key="cell_type",
    unlabeled_category="Unknown"
)
scanvi_model.train(max_epochs=20)
```

See the `scvi-scanvi` skill for full scANVI workflows.

## Intermediate Usage

### Library Size Handling

For most datasets, the default (observed library size) works well:

```python
model = scvi.model.SCVI(adata, use_observed_lib_size=True)  # default
```

For details on when to model library size as a latent variable or use manual size
factors, see `references/library_size_handling.md`.

### Parameter Tuning

For guidance on choosing `n_latent`, `n_layers`, and `gene_likelihood` based on
dataset size and structure, see `references/training_parameters.md`.

## Advanced Usage

### Differential Expression

scVI provides Bayesian DE with uncertainty quantification. This is the correct
way to do cross-condition DE with integrated data (not `sc.tl.rank_genes_groups`
on raw counts, which confounds batch with biology).

```python
de = model.differential_expression(
    groupby="cell_type",
    group1="T_cells",
    group2="B_cells",
    mode="change",          # Tests if LFC > delta
    delta=0.25,
    batch_correction=True   # Account for batch effects
)
markers = de[(de["lfc_mean"] > 0) & (de["bayes_factor"] > 3)]
```

For full DE workflows including 1-vs-all marker discovery and result
interpretation, see `references/differential_expression.md`.

## Verification Checklist

After scVI integration, verify ALL of the following:

- [ ] **Raw counts input**: `adata.layers["counts"]` contains integers, not floats
- [ ] **No NaN in latent space**: `np.isnan(adata.obsm["X_scVI"]).sum() == 0`
- [ ] **Training converged**: Plot `model.history["elbo_train"]` -- loss should plateau in the last 10-20% of epochs
- [ ] **Batch mixing on UMAP**: Color UMAP by batch -- batches should interleave, not form isolated islands
- [ ] **Biology preserved**: Color UMAP by known marker genes -- expected populations should still separate
- [ ] **Latent dimensionality**: If clusters are too coarse, try increasing `n_latent`; if too fine, decrease it
- [ ] **Checkpoint saved**: `model.save("scvi_model/")` and `adata.write_h5ad("integrated.h5ad")`

For a guided visual checklist, see `checks/check_integration_quality.md`.
For pre-training validation, run: `python scripts/validate_before_scvi.py data.h5ad`

## Common Pitfalls

### Pitfall: Passing normalized data to scVI
**Symptom:** Training loss is unusually low from the start; integration looks poor or over-corrected.
**Cause:** scVI assumes raw integer counts with a negative binomial likelihood. Normalized or log-transformed data violates this assumption fundamentally.
**Fix:** Always use `layer="counts"` pointing to raw counts. Verify:
```python
data = adata.layers["counts"]
vals = data.data if hasattr(data, 'data') else data.ravel()
assert np.allclose(vals, vals.astype(int)), "Counts layer contains non-integer values"
```

### Pitfall: Using scVI latent space for naive cross-batch DE
**Symptom:** DE results dominated by thousands of genes, many of which are housekeeping.
**Cause:** scVI corrects the latent space, but running `sc.tl.rank_genes_groups` on raw counts still sees batch effects. The integration helps visualization but does not fix the counts matrix.
**Fix:** Use `model.differential_expression()` (batch-aware Bayesian DE) or pseudobulk DE with batch as a covariate. See the `batch-aware-analysis` practice skill.

### Pitfall: Too few highly variable genes
**Symptom:** Rare cell types or subtypes with low-abundance markers are merged on UMAP.
**Cause:** HVG selection with `n_top_genes` < 1500 may exclude markers critical for distinguishing subtypes.
**Fix:** Use 2000-4000 HVGs. Check that key markers survive selection:
```python
critical_markers = ["CLEC9A", "CD1C", "CLEC4C", "AXL"]
missing = [g for g in critical_markers if g not in adata.var_names]
if missing:
    print(f"WARNING: Markers excluded by HVG filter: {missing}")
```

### Pitfall: Not saving the model
**Symptom:** Need to recompute latent space or run DE, but model was garbage-collected.
**Cause:** scVI models exist only in memory unless explicitly saved.
**Fix:** Always save after training:
```python
model.save("scvi_model/")
# To reload later:
model = scvi.model.SCVI.load("scvi_model/", adata=adata)
```

## Complementary Skills

| When you need... | Use skill | Relationship |
|-----------------|-----------|-------------|
| Understand .h5ad structure | `anndata` | Prerequisite |
| QC before integration | `single-cell-rna-qc` | Prerequisite |
| Installation and shared setup | `scvi-framework` | Foundation |
| Clustering, UMAP, visualization | `scanpy` | Next step |
| Add labels for better integration | `scvi-scanvi` | Extension |
| Map new data to trained model | `scvi-scarches-reference-mapping` | Extension |
| Pre-trained models from HuggingFace | `scvi-hub-models` | Alternative starting point |
| Verify DE is not batch-driven | `batch-aware-analysis` | Complementary practice |

## Resources

- [scVI User Guide](https://docs.scvi-tools.org/en/stable/user_guide/models/scvi.html)
- [scVI DE Tutorial](https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scrna/scVI_DE_worm.html)
- [Lopez et al. 2018 (original paper)](https://www.nature.com/articles/s41592-018-0229-2)
```

### 10.3 Example Reference: `references/training_parameters.md`

This file would contain ~200 lines expanding on:
- How `n_latent` affects cluster granularity (with rules of thumb by dataset size)
- `n_layers` trade-offs (deeper = more expressive but slower, risk of overfitting)
- `gene_likelihood`: when `nb` vs `zinb` matters (datasets with true zero inflation)
- `dispersion` modes: per-gene vs per-gene-batch (when batches have very different protocols)
- Training hyperparameters: `max_epochs`, `batch_size`, `early_stopping`
- GPU memory estimation: `~4 bytes * n_cells * n_genes * n_layers`
- How to compare two parameter configurations systematically

### 10.4 Example Script: `scripts/validate_before_scvi.py`

```python
#!/usr/bin/env python3
"""Validate that an AnnData object is ready for scVI training.

Checks: counts layer exists and contains integers, batch key exists,
HVGs are selected, sufficient cells and batches.

Usage:
    python validate_before_scvi.py data.h5ad
    python validate_before_scvi.py data.h5ad --layer counts --batch-key batch

Exit codes:
    0: Validation passed
    1: Validation failed (issues found)
    2: Usage error
"""
import argparse
import sys
import anndata as ad
import numpy as np


def validate(path: str, layer: str, batch_key: str) -> int:
    """Validate AnnData for scVI. Returns 0 if valid, 1 if not."""
    adata = ad.read_h5ad(path)
    errors = []
    warnings = []

    print(f"Shape: {adata.n_obs} cells x {adata.n_vars} genes")

    # Check counts layer
    if layer not in adata.layers:
        errors.append(f"Layer '{layer}' not found. Available: {list(adata.layers.keys())}")
    else:
        data = adata.layers[layer]
        vals = data.data if hasattr(data, "data") else data.ravel()
        if (vals < 0).any():
            errors.append(f"Layer '{layer}' contains negative values")
        if not np.allclose(vals, vals.astype(int)):
            errors.append(f"Layer '{layer}' contains non-integer values")
        else:
            print(f"Counts layer '{layer}': OK (integers, non-negative)")

    # Check batch key
    if batch_key not in adata.obs.columns:
        errors.append(
            f"Batch key '{batch_key}' not in obs. "
            f"Available: {adata.obs.columns.tolist()}"
        )
    else:
        n_batches = adata.obs[batch_key].nunique()
        batch_sizes = adata.obs[batch_key].value_counts()
        print(f"Batches: {n_batches} ('{batch_key}')")
        print(f"  Smallest batch: {batch_sizes.min()} cells")
        print(f"  Largest batch:  {batch_sizes.max()} cells")
        if n_batches < 2:
            errors.append("Only 1 batch found -- scVI integration is unnecessary")
        if batch_sizes.min() < 50:
            warnings.append(
                f"Smallest batch has only {batch_sizes.min()} cells -- "
                "may cause training instability"
            )

    # Check HVGs
    if "highly_variable" in adata.var.columns:
        n_hvg = adata.var["highly_variable"].sum()
        print(f"HVGs: {n_hvg}")
        if n_hvg < 1000:
            warnings.append(f"Only {n_hvg} HVGs -- recommend at least 1500")
    else:
        warnings.append("No HVG selection detected in adata.var")

    # Check cell count
    if adata.n_obs < 500:
        warnings.append(
            f"Only {adata.n_obs} cells -- scVI may not train well. "
            "Consider classical methods."
        )

    # Report
    if warnings:
        print("\nWARNINGS:")
        for w in warnings:
            print(f"  - {w}")

    if errors:
        print("\nVALIDATION FAILED:")
        for e in errors:
            print(f"  - {e}")
        return 1
    else:
        print("\nVALIDATION PASSED: Data is ready for scVI")
        return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Validate AnnData for scVI training"
    )
    parser.add_argument("path", help="Path to .h5ad file")
    parser.add_argument(
        "--layer", default="counts", help="Layer containing raw counts (default: counts)"
    )
    parser.add_argument(
        "--batch-key", default="batch", help="Batch column in obs (default: batch)"
    )
    args = parser.parse_args()
    sys.exit(validate(args.path, args.layer, args.batch_key))
```

### 10.5 Example Check: `checks/check_integration_quality.md`

```markdown
# Integration Quality Checklist

Run these checks after scVI training. Each section is a self-contained verification.

## 1. Training Convergence

```python
import matplotlib.pyplot as plt

train_loss = model.history["elbo_train"]["elbo_train"]
plt.plot(train_loss)
plt.xlabel("Epoch")
plt.ylabel("ELBO")
plt.title("Training Loss")
plt.show()

# Quantitative check: last 10% of epochs should be stable
last_10pct = train_loss[int(len(train_loss) * 0.9):]
variation = (last_10pct.max() - last_10pct.min()) / abs(last_10pct.mean())
print(f"Loss variation in last 10% of training: {variation:.4f}")
assert variation < 0.01, "Training may not have converged -- consider more epochs"
print("PASS: Training converged")
```

## 2. Batch Mixing

```python
import scanpy as sc

sc.pl.umap(adata, color="batch", title="Colored by batch -- should be interleaved")
```

**Pass:** Batches interleave across clusters. No batch forms its own isolated island.
**Fail:** Batches form separate islands. Try: increase `n_latent`, switch to `gene_likelihood="nb"`,
check that `batch_key` was correctly specified during `setup_anndata`.

## 3. Biological Signal Preservation

```python
# Adjust markers for your system
marker_genes = ["CD3D", "CD14", "MS4A1", "FCGR3A"]
available = [g for g in marker_genes if g in adata.var_names]
sc.pl.umap(adata, color=available, use_raw=True, ncols=2,
           title=[f"{g} expression" for g in available])
```

**Pass:** Markers localise to distinct UMAP regions matching expected biology.
**Fail:** Markers are spread uniformly or co-localise unexpectedly. Integration
may have over-corrected. Try reducing `n_latent` or using fewer batches.

## 4. Latent Space Sanity

```python
import numpy as np

latent = adata.obsm["X_scVI"]
print(f"Shape: {latent.shape}")
print(f"NaN count: {np.isnan(latent).sum()}")
print(f"Range: [{latent.min():.2f}, {latent.max():.2f}]")
print(f"Mean: {latent.mean():.4f}")
print(f"Std:  {latent.std():.4f}")

assert not np.isnan(latent).any(), "FAIL: NaN values in latent space"
assert not np.isinf(latent).any(), "FAIL: Inf values in latent space"
print("PASS: Latent space is finite and complete")
```

## 5. Checkpoint

```python
# Save everything needed to resume
model.save("scvi_model/")
adata.write_h5ad("integrated.h5ad")
print(f"Model saved to scvi_model/")
print(f"AnnData saved to integrated.h5ad ({adata.shape})")
```
```

### 10.6 What Changed from v2 to v3

| Aspect | v2 (`scvi-basic.md`) | v3 (`scvi-basic/`) |
|--------|---------------------|-------------------|
| Format | Single 156-line .md file | Directory with SKILL.md (~300 lines) + 3 reference files + 1 script + 1 check |
| Frontmatter | `name`, `description` only | Adds `category`, `tier`, `tags`, `complementary-skills`, `contraindications` |
| Decision tree | Absent | Present: routes user to correct skill before they start |
| Quick Start verification | Absent | Present: assertions + visual check |
| Progressive depth | Flat (all sections at same level) | Three tiers: Basic, Intermediate, Advanced |
| Verification checklist | Absent | 7-item checklist with executable checks |
| Pitfalls format | Table (Issue / Solution) | Symptom / Cause / Fix with code |
| Complementary skills | One-line "See X.md" references | Structured table with relationship types |
| Executable validation | None | `scripts/validate_before_scvi.py` |
| Deep-dive content | All inline | Factored to `references/` (training params, DE, library size) |
| Post-training guide | None | `checks/check_integration_quality.md` |

---

## Appendix A: Full Frontmatter Schema

```yaml
---
# REQUIRED
name: string                     # kebab-case identifier, matches directory name
description: string              # 1-3 sentences, include trigger words
category: enum                   # foundation | integration | annotation | analysis | workflow | practice
tier: enum                       # simple | standard | rich

# RECOMMENDED
tags: [string]                   # Free-form discovery tags
complementary-skills: [string]   # Skill names that pair well
contraindications: [string]      # When NOT to use (1 sentence each)

# OPTIONAL
version: string                  # semver (e.g., "1.0.0")
upstream-docs: string            # URL to official tool documentation
license: string                  # License identifier (e.g., "MIT")
metadata:
  skill-author: string           # Who wrote this skill
  last-reviewed: date            # YYYY-MM-DD
  website: string                # URL
---
```

## Appendix B: Category-Tier Matrix

Expected distribution of tiers across categories. Not a hard rule, but a guide for consistency.

| Category | Simple | Standard | Rich |
|----------|--------|----------|------|
| Foundation | Common | Common | Rare |
| Integration | Rare | Common | Common |
| Annotation | Rare | Common | Common |
| Analysis | Common | Common | Rare |
| Workflow | Never | Never | Always |
| Practice | Common | Common | Never |

## Appendix C: Migration Checklist for a Single Skill

Use this checklist when migrating a flat `.md` skill to v3 directory format.

- [ ] Create `skills/skill-name/` directory
- [ ] Move `skills/skill-name.md` to `skills/skill-name/SKILL.md`
- [ ] Add missing frontmatter: `category`, `tier`, `complementary-skills`, `contraindications`
- [ ] Add `## Decision Tree` section (required for Standard and Rich tiers)
- [ ] Add `## Verification Checklist` section with executable checks
- [ ] Ensure Quick Start includes a "Verify it worked:" block
- [ ] Add or restructure `## Common Pitfalls` with Symptom/Cause/Fix format
- [ ] Add `## Complementary Skills` table with Relationship column
- [ ] If file exceeds ~400 lines, factor deep-dive content into `references/`
- [ ] Update `tier` in frontmatter to match actual directory contents
- [ ] Test: `activate-role.sh <role> --project-dir <project>` creates correct symlink
- [ ] Test: Symlink in `.claude/skills/` resolves and is readable by AI assistant
- [ ] Test: Relative paths in SKILL.md (e.g., `references/topic.md`) resolve through symlink

## Appendix D: Section Template Quick Reference

For quick copy-paste when creating new skills.

### Tool Skill (Foundation/Integration/Annotation/Analysis)

```markdown
---
name:
description:
category:
tier:
tags: []
complementary-skills: []
contraindications: []
---

# Title

## Overview

## Decision Tree

## Quick Start

**Verify it worked:**

## Basic Usage

## Intermediate Usage

## Advanced Usage

## Verification Checklist

## Common Pitfalls

## Complementary Skills

## Resources
```

### Practice Skill

```markdown
---
name:
description:
category: practice
tier: simple
tags: [methodology]
complementary-skills: []
contraindications: []
---

# Title

## The Problem

## The Principle

## The Pattern

## When to Apply

## When to Skip

## Signs You Skipped It
```

### Workflow Skill

```markdown
---
name:
description:
category: workflow
tier: rich
tags: [workflow]
complementary-skills: []
contraindications: []
---

# Title

## Overview

## Prerequisites

## Steps

### Step 1: [Name]
**Skill:** `skill-name`

#### Gate 1

### Step 2: [Name]
...

## Checkpoints

## Resume Protocol
```

---

*End of Skills Format v3 specification.*
