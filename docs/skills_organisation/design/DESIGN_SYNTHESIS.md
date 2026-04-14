# SciAgent-toolkit v3: Design Synthesis

> **Date:** 2026-03-07
> **Status:** Synthesis of design documents
> **Source documents:**
> - `ARCHITECTURE_v3.md` — Overall tier system, component architecture, migration
> - `SKILLS_FORMAT_v3.md` — Evolved skill format, categories, practice skills
> - `WORKFLOW_SYSTEM.md` — Phase-based workflows, verification gates, checkpoints
> - `OUTPUT_STYLES_AND_KNOWLEDGE.md` — Interaction modes, knowledge persistence

---

## What Problem Are We Solving?

AI-assisted bioinformatics sessions drift. They become unpredictable, over-engineered, and disconnected from prior work. The user wants:

1. **To stay in the loop** — AI proposes, human verifies. No castles of cards.
2. **To build confidence progressively** — Small samples first, understand before scaling.
3. **To accumulate knowledge** — Sessions build on each other, not start from zero.
4. **To have reproducible workflows** — Same analysis structure, regardless of AI provider.
5. **To avoid over-engineering** — Simple, declarative code with rigorous verification.

---

## Architecture Overview

### The 6-Tier System (evolved from 4-tier)

```
Tier 0: KNOWLEDGE BASE     persists across all sessions and workflows
Tier 1: WORKFLOWS           persists across sessions within a project
Tier 2: ROLES               persists within a workflow phase
Tier 3: SKILLS              loaded by role + phase
Tier 4: AGENTS              loaded by role
Tier 5: PROFILES + ADDONS   what MCP servers are connected
```

Key change: **Workflow** becomes the top-level organizing concept (above Role). A workflow defines a multi-phase progression with verification gates. Roles become the interaction mode *within* a phase.

### How the Cascade Works

```
Workflow (scrna-atlas)
  └── Phase (explore → integrate → annotate → analyse)
       └── Role (mentor / partner / operator)
            └── Skills (anndata, scanpy, scvi-basic...)
                 └── Agents (bio-interpreter, insight-explorer...)
                      └── Profile + Addons (coding + notebook-tools)

All tiers READ/WRITE → Knowledge Base
```

---

## The Four New Concepts

### 1. Workflows (WORKFLOW_SYSTEM.md)

A workflow is a YAML contract declaring phases, verification gates, and decision points. Not executable — the AI reads it to know where the user is and what must be true before advancing.

**Core elements per phase:**
- Prerequisites (which prior gates must pass)
- Skills activated
- Output style (how the AI interacts)
- Scale constraint (`small` / `full` / `small_then_full`)
- Verification gates (checklists, human judgments, computational checks)
- Decision points (where the human MUST choose)
- Checkpoints (what to save for resume)

**Anti-drift mechanism:** The AI reads `workflow_state.yaml` at session start to know: what phase we're in, what gates have passed, what decisions were made. It does not auto-advance.

**Small-sample testing:** Phases can declare `scale: small` with a `scale_spec` that defines subsample strategy and graduation criteria. The user tests on 5K cells before committing to 1M.

**Implementation:** Three files total — one workflow YAML, one state YAML, one CLAUDE.md instruction. No new scripts or MCP servers required.

### 2. Evolved Skills (SKILLS_FORMAT_v3.md)

Skills evolve from flat `.md` files to directories:

```
skill-name/
    SKILL.md          # Entry point (always read)
    references/       # Deep-dive guides (read on demand)
    scripts/          # Executable utilities
    assets/           # Templates, configs
    checks/           # Verification scripts/checklists
```

**Backward compatible:** Flat `.md` files still work. `activate-role.sh` detects format.

**Every skill must have:**
- Decision tree (when to use this vs alternatives)
- Quick start with verification step
- Common pitfalls (Symptom/Cause/Fix format)
- Verification checklist

**Two new skill categories:**

| Category | Purpose | Examples |
|----------|---------|---------|
| **Practice** | Coding habits and mental models | `checkpoint-verify`, `small-sample-test`, `batch-aware-analysis`, `separation-of-concerns`, `declarative-exploration` |
| **Workflow** | Multi-step procedures composing other skills | `scrna-eda-workflow`, `integration-workflow`, `annotation-workflow` |

Practice skills are the codified lessons from real failures (like the v2 annotation disaster). They answer "how to think" not "how to use a tool."

### 3. Output Styles (OUTPUT_STYLES_AND_KNOWLEDGE.md, Part A)

Five interaction modes, role-aware and phase-aware:

| Style | Autonomy | When | Key Behavior |
|-------|----------|------|-------------|
| **mentor** | Low | Learning, exploration | Socratic, one-cell-at-a-time, "what do you think?" |
| **architect** | Low | Planning, design | Trade-offs, ASCII diagrams, decision records |
| **partner** | Medium | Implementation, analysis | Scaffold+verify, 1-3 cells, explain non-obvious |
| **reviewer** | Medium | Validation, QC | Skeptical, biological plausibility, evidence-based |
| **operator** | High | Production, routine | Batch execution, report results, minimal explanation |

**Switching:** Explicit command, phase-based suggestion, or implicit signals ("just do it" → operator). Style persists for the session and is recorded in workflow state.

**Provider-agnostic:** Styles are Markdown content. Claude Code uses output-styles directory; Gemini CLI gets the body injected into system instructions.

### 4. Knowledge Base (OUTPUT_STYLES_AND_KNOWLEDGE.md, Part B + ARCHITECTURE_v3.md)

Structured persistence replacing ad-hoc MEMORY.md and handoff docs:

```
.knowledge/ (or .sciagent/knowledge/)
    decisions.yaml        # D001: "Use scVI over BBKNN because..."
    verifications.yaml    # V001: gate passed/failed with evidence
    workflow_state.yaml   # Current phase, gates, decisions, active style
    patterns.md           # What worked, what didn't
```

**Session startup protocol:** AI reads workflow_state → decisions → verifications → MEMORY.md, then greets with context: "Last session we were in Phase 2 (integration understanding), partner style. We'd decided on scVI with n_latent=30. Continue?"

**Boundary with MEMORY.md:**
- MEMORY.md = sticky notes (quick facts, API gotchas, path shortcuts)
- Knowledge base = lab notebook (decisions with rationale, verification evidence, phase tracking)
- They coexist. MEMORY.md for speed; knowledge base for rigor.

---

## Migration Path

Incremental, each step independently useful:

| Step | What Changes | Effort |
|------|-------------|--------|
| 0 | Nothing — v2 continues to work | None |
| 1 | Convert high-value skills to directory format | Medium |
| 2 | Add output styles to roles | Low |
| 3 | Create `.knowledge/` directory convention | Low |
| 4 | Define workflow YAML for current project | Medium |
| 5 | Add verification gates to workflow | Low |
| 6 | Curate toolkit-level knowledge | Ongoing |

No v2 files deleted. All new fields optional. New formats coexist with old.

---

## What This Means for the DC_hum_verse Project

With this design, the user's described workflow maps directly to:

1. **Activate workflow:** `scrna-atlas` with 7 phases
2. **Phase 0 (explore):** mentor style, `anndata` + `scanpy` + `declarative-exploration` skills, gate = "each dataset explored independently"
3. **Phase 1 (integration planning):** mentor style, decision point = "what to integrate?", gate = "subsetting targets verified"
4. **Phase 2 (integration understanding):** mentor style, `scale: small`, scVI on 5K cells, gate = "user can explain what scVI changes"
5. **Phase 3 (full integration):** partner style, validated parameters from Phase 2, gate = "batch mixing confirmed"
6. **Phase 4 (manual annotation):** mentor → partner style, gate = "marker genes match, no catch-all clusters"
7. **Phase 5 (semi-supervised annotation):** partner style, `scale: small_then_full`, gate = "automated matches manual >90%"
8. **Phase 6 (deep analysis):** partner → operator style, only after full annotation confidence

Every phase transition requires gate passage. Every decision is logged. Every session picks up where the last left off.

---

## Design Principles Summary

1. **Human-in-the-loop by design** — Verification gates, decision points, scale constraints
2. **Progressive complexity** — Small samples → understanding → full scale
3. **Reproducibility through structure** — Workflows + knowledge base + decision logs
4. **Knowledge accumulates** — `.knowledge/` persists decisions, findings, verification evidence
5. **Provider-agnostic** — Workflows, skills, knowledge in plain YAML/Markdown
6. **Small surface area** — 3 new files per project (workflow.yaml, state.yaml, .knowledge/)

---

## File Index

| Document | Lines | Covers |
|----------|-------|--------|
| `ARCHITECTURE_v3.md` | ~700 | 6-tier system, component specs, provider agnosticism, migration |
| `SKILLS_FORMAT_v3.md` | ~800 | Directory format, SKILL.md spec, tiers, categories, practice skills, worked example |
| `WORKFLOW_SYSTEM.md` | ~400 | Phase definitions, verification gates, small-sample testing, checkpoints, example sessions |
| `OUTPUT_STYLES_AND_KNOWLEDGE.md` | ~600 | 5 output styles, switching mechanics, knowledge architecture, decision logs, session continuity |
| `DESIGN_SYNTHESIS.md` | This file | Executive summary tying it all together |

---

*End of synthesis. Next step: review these designs and decide which components to implement first.*
