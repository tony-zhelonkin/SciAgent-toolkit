# Web Research: Provider-Agnostic Agentic Science — External Best Practice
**Date:** 2026-06-23  
**Scope:** Per-repo agent context management, Claude Code patterns, spec-driven dev, reproducibility, figure legibility  
**Purpose:** Inform SciAgent-toolkit design

---

## Area 1 — The AGENTS.md Standard and Siblings

### Background and Governance

**AGENTS.md was released by OpenAI in August 2025** as an open, Markdown-based standard for providing project-specific instructions to AI coding agents. On **December 9, 2025**, OpenAI donated it to the **Agentic AI Foundation (AAIF)**, a directed fund under the Linux Foundation, alongside Anthropic donating the Model Context Protocol (MCP) and Block donating the Goose agent framework. Platinum members include AWS, Anthropic, Block, Bloomberg, Cloudflare, Google, Microsoft, and OpenAI.

- **Source:** [Linux Foundation AAIF Announcement](https://www.linuxfoundation.org/press/linux-foundation-announces-the-formation-of-the-agentic-ai-foundation) — 2025-12-09  
  *Why it matters: AGENTS.md is now a neutral, foundation-governed open standard — not a single-vendor lock-in.*

**Adoption as of 2026:** 60,000+ open-source projects; read natively by Claude Code, OpenAI Codex CLI, Cursor, GitHub Copilot, Windsurf, Gemini CLI, Amazon Q, Devin, Aider, Sourcegraph Amp, Zed AI, Continue, Roo Code, Factory Droids, and Google Jules.

- **Source:** [BuildBetter AGENTS.md Complete Guide (2026)](https://blog.buildbetter.ai/agents-md-complete-guide-for-engineering-teams-in-2026/)  
  *Why it matters: AGENTS.md is the closest thing to a cross-tool standard — putting it in a repo means 28+ tools read it.*

### Tool-by-Tool Comparison (as of 2026)

| File | Owner/Origin | Status (2026) | Cross-Tool? |
|------|-------------|---------------|-------------|
| `AGENTS.md` | OpenAI → Linux Foundation (AAIF) | De facto standard, 60K+ repos | Yes (28+ tools) |
| `CLAUDE.md` | Anthropic | Active; richer 3-layer memory model | Claude Code only |
| `GEMINI.md` | Google | Active in Gemini CLI | Gemini CLI only |
| `.cursorrules` / `.cursor/rules/*.mdc` | Cursor | `.cursorrules` deprecated; `.mdc` is current | Cursor only |
| `.github/copilot-instructions.md` | GitHub/Microsoft | Supported; also `instructions/*.instructions.md` since July 2025 | Copilot only |
| `.windsurfrules` / `.windsurf/rules/` | Windsurf/Codeium | Directory system, 12K char limit per file | Windsurf only |
| `SKILL.md` | Open (AAIF) | Growing support across all major tools | Most major tools |

- **Source:** [DeployHQ Config Files Guide](https://www.deployhq.com/blog/ai-coding-config-files-guide) — 2025  
  *Why it matters: Each tool requires different files; knowing the matrix prevents maintaining 6 copies of the same instructions.*

### CLAUDE.md Specifics (Anthropic Official)

Claude Code uses a **3-layer memory model**:
1. `~/.claude/CLAUDE.md` — global user settings
2. `CLAUDE.md` at project root — project-level (committed to repo)
3. `CLAUDE.local.md` — personal local overrides (gitignored)

Later files override earlier ones. Claude Code also supports:
- `@import` directives to reference external files
- Path-scoped rules (frontmatter `paths:`) that load only when matching files are active
- Subdirectory `CLAUDE.md` files loaded on-demand
- Organizational deployment via MDM using centrally-managed CLAUDE.md files

**Anthropic's official recommendation:** Keep root CLAUDE.md under 200 lines; give it an owner; move 30+ line procedures into skills instead.

- **Source:** [Anthropic — Steering Claude Code blog](https://claude.com/blog/steering-claude-code-skills-hooks-rules-subagents-and-more)  
  *Why it matters: SciAgent-toolkit CLAUDE.md should be a thin pointer + project conventions, not a monolith.*

### One Source of Truth Strategy

The **canonical multi-tool pattern** adopted by practitioners in 2026:
1. Write `AGENTS.md` as the single source of truth (under 150 lines; imperative, specific language)
2. Make `CLAUDE.md` a one-liner: `Strictly follow the rules in ./AGENTS.md`  
   (or use a symlink: `ln -s AGENTS.md CLAUDE.md`)
3. Tool-specific files (`.cursor/rules/`, `.github/copilot-instructions.md`) contain only features unique to that tool (e.g., Cursor's glob-activated scoped rules)
4. For monorepos: thin root `AGENTS.md` for org-wide standards + package-specific `AGENTS.md` files

```
project-root/
├── AGENTS.md                          # Universal source of truth
├── AGENTS.override.md                 # Local-only overrides (gitignored)
├── CLAUDE.md                          # → symlink to AGENTS.md or one-liner reference
├── .github/
│   └── copilot-instructions.md        # Copilot-specific only
└── .cursor/
    └── rules/                         # Cursor glob-scoped rules only
```

- **Source:** [DeployHQ Blog](https://www.deployhq.com/blog/ai-coding-config-files-guide) + [Codex Knowledge Base (2026-05-27)](https://codex.danielvaughan.com/2026/05/27/agent-instruction-files-agents-md-claude-md-cross-tool-portability-codex-cli/)  
  *Why it matters: Prevents maintenance drift across tools; one update propagates everywhere.*

### Key Pitfalls

- Instructions at the bottom of long files are **ignored 43% more often** than top instructions ("lost in the middle" effect)
- Files exceeding ~1,200 instruction tokens show "sharp enforcement degradation"
- Keep AGENTS.md under 150–500 lines; prioritize the most critical conventions first
- Per an ETH study, LLM-generated instruction files (not human-reviewed) reduced task success rates by 0.5–2% while increasing inference costs by 20%+

- **Source:** [AgentRuleGen — Cursor vs CLAUDE.md vs Copilot](https://www.agentrulegen.com/guides/cursorrules-vs-claude-md) — 2026  
  *Why it matters: SciAgent-toolkit should pre-populate a short, dense AGENTS.md rather than auto-generating a verbose one.*

### SKILL.md (Open Standard, AAIF)

As of 2026, `SKILL.md` files in `.claude/skills/<name>/SKILL.md` (or equivalent locations) are supported by Claude Code, Codex CLI, Gemini CLI, GitHub Copilot, Cursor, Cline, Windsurf, and OpenCode. Skills use YAML frontmatter (name, description) with Markdown instruction body. The metadata loads at session start; the full body loads only when invoked. This enables **progressive disclosure** of reusable workflows.

- **Source:** [AGENSI.io — SKILL.md Open Standard](https://www.agensi.io/learn/agent-skills-open-standard) — 2026  
  *Why it matters: Reusable science workflows (e.g., "run-analysis", "generate-figure") can be packaged as skills that work across all supported tools.*

---

## Area 2 — Claude Code: Skills, Subagents, Hooks, and Orchestration Patterns

### The Core Extensibility Stack

Claude Code provides three extensibility layers:
1. **Skills** — reusable procedural instruction packages
2. **Subagents** — isolated parallel worker contexts
3. **Hooks** — deterministic lifecycle automation

- **Source:** [Anthropic — Steering Claude Code](https://claude.com/blog/steering-claude-code-skills-hooks-rules-subagents-and-more)  
  *Why it matters: These three layers together replace ad-hoc prompting with structured, reusable, auditable workflows.*

### Skills

Skills live in `.claude/skills/<name>/SKILL.md` with YAML frontmatter (`name`, `description`) and a procedural Markdown body. Slash commands (`.claude/commands/`) are unified with skills as of 2026 — each skill gets a `/slash-command` interface.

**What belongs in a skill (not CLAUDE.md):** deploy workflows, release checklists, review processes, analysis pipelines — anything procedural and multi-step.

**Triggering:** via slash command (e.g., `/run-analysis`), auto-matched by task description, or explicit in a system prompt. Only the name + description loads at session start; full body loads on invocation (context-efficient).

**On compaction:** invoked skills re-inject up to a shared budget, dropping oldest first.

- **Source:** [AlexOp.dev — Claude Code Customization Guide](https://alexop.dev/posts/claude-code-customization-guide-claudemd-skills-subagents/) — 2025  
  *Why it matters: SciAgent-toolkit analysis workflows should be packaged as skills, not buried in a monolithic CLAUDE.md.*

### Subagents (Orchestrator-Worker Pattern)

Subagents are autonomous workers spawned in **fresh, isolated context windows**. Only their final message returns to the main session. They nest up to **5 levels deep** and can scale to **tens to hundreds** via dynamic `/workflows`.

Defined in `.claude/agents/<name>.md` with YAML frontmatter specifying: `model`, `tools`, `disallowedTools`, `permissionMode`, `mcpServers`, `hooks`, `maxTurns`, `skills`, and more.

**Model tiering (cost control):**
- Use `haiku` (cheap, fast) for exploration subagents: file search, quick questions, log scanning
- Escalate to `sonnet` or `opus` for planning, code generation, and review
- Opus coordinates multi-agent teams via agent templates (released May 5, 2026)

**"Master-Clone" architecture (practitioner pattern):** Put all key context in CLAUDE.md, let the main Opus agent delegate to self-copies. The orchestrator manages its own subagent topology dynamically.

**Massive Parallel Scripting:** For large-scale refactors, call `claude -p "..."` in parallel via bash scripts — more scalable than managing dozens of subagent tasks from within an agent.

- **Source:** [Anthropic — Steering Claude Code](https://claude.com/blog/steering-claude-code-skills-hooks-rules-subagents-and-more) + [ClaudeCode Full Stack — AlexOp.dev](https://alexop.dev/posts/understanding-claude-code-full-stack/) — 2025/2026  
  *Why it matters: Long-running science pipelines (QC → analysis → figure → caption) are natural candidates for subagent decomposition.*

### Hooks

Hooks are deterministic event-driven scripts that fire at **25 lifecycle points** including:
- `PreToolUse` — fires before any tool call; can block with exit code 2 (primary safety checkpoint)
- `UserPromptSubmit` — fires when user submits a prompt; can block or modify
- `PermissionRequest` — fires on permission requests; can auto-approve or deny
- `Stop` / `SubagentStop` — fires when main session or subagent finishes
- `PreCompact` — fires before context compaction

Hook types: `command`, `HTTP`, `mcp_tool`, `prompt`, `agent`. Hooks run deterministically outside the main context window — **they cannot hallucinate**.

**Use for SciAgent-toolkit:** Auto-run linting after edits, enforce "no ephemeral scripts" policy by blocking writes to `/tmp`, auto-commit analysis artifacts, trigger provenance logging on `Stop`.

- **Source:** [OFox.ai — Claude Code Hooks Subagents Skills Guide](https://ofox.ai/blog/claude-code-hooks-subagents-skills-complete-guide-2026/) — 2026  
  *Why it matters: Hooks are the mechanism to enforce reproducibility rules automatically without relying on the LLM's instruction-following.*

### Decision Framework

| Mechanism | Use When |
|-----------|----------|
| **CLAUDE.md / AGENTS.md** | Always-on context: repo conventions, tech stack, prohibited patterns |
| **Skill** | Reusable multi-step procedures that stay in main thread (visible, steerable) |
| **Subagent** | Side tasks that would pollute main context; parallel isolated work |
| **Hook** | Deterministic enforcement; lifecycle automation; cannot be overridden by LLM |
| **Slash command** | User-triggered, explicit workflow invocation |

---

## Area 3 — Spec-Driven / Plan-Driven Development with Agents

### The Core Methodology

**Spec-Driven Development (SDD)** emerged in 2025 as a formal response to "vibe coding" failure modes (hallucination, scope drift, decay at scale). Core principle: **separate planning from implementation**. AI agents analyze requirements and generate specification artifacts (Markdown files); a human reviews and validates; then a coding agent generates product code constrained by the spec.

- **Source:** [Thoughtworks — Spec-Driven Development (2025)](https://www.thoughtworks.com/en-us/insights/blog/agile-engineering-practices/spec-driven-development-unpacking-2025-new-engineering-practices)  
  *Why it matters: SDD is Thoughtworks' named methodology for AI-assisted engineering — gives SciAgent-toolkit a credible external framework to reference.*

### The Plan → Implement → Review Loop

Standard 3-phase structure:
1. **Plan:** Requirements → spec (Markdown) → implementation plan → step decomposition
2. **Implement:** Step-by-step execution, each step tested before the next
3. **Review:** Compare output against spec; deterministic CI/CD checks; human validation

**Sizing steps:** Each step should be small enough to implement safely with strong testing but big enough to move the project forward — no orphaned code; each prompt builds on the previous.

**TDD is the best anti-hallucination mechanism:** "The robots LOVE TDD — build test and mock first, then make the mock real." (Harper Reed, Feb 2025)

- **Source:** [Harper Reed — My LLM Codegen Workflow (2025-02-16)](https://harper.blog/2025/02/16/my-llm-codegen-workflow-atm/)  
  *Why it matters: The TDD + spec loop is the most-cited practical anti-hallucination technique for agentic coding.*

### Harper Reed's Artifact Pattern (Feb 2025)

Canonical artifacts committed to the repo root:
- `spec.md` — comprehensive developer-ready specification (output of conversational LLM + iterative Q&A)
- `prompt_plan.md` — detailed step-by-step blueprint with executable prompts for each step
- `todo.md` — checklist generated from `prompt_plan.md`

Planning timeline: ~15 minutes for greenfield spec; ~8–12 steps typical plan; 30–45 min end-to-end for plan + initial implementation regardless of project complexity.

- **Source:** [Harper Reed Blog (2025-02-16)](https://harper.blog/2025/02/16/my-llm-codegen-workflow-atm/)  
  *Why it matters: Concrete file-naming convention that can be adopted directly by SciAgent-toolkit as the planning artifact structure.*

### Named Methodologies and Tools (2025–2026)

| Name | Origin | Key Distinguishing Feature |
|------|---------|---------------------------|
| **Harper Reed workflow** | Individual practitioner (Feb 2025) | spec.md + prompt_plan.md; TDD; 3-phase personal workflow |
| **GitHub Spec Kit** | GitHub (Sept 2025) | Multi-step: Constitution → Specify → Plan → Tasks → Implement → PR; OSS; tool-agnostic |
| **BMAD Method** | OSS community (48K GitHub stars) | Full SDLC with named AI personas (analyst, architect, dev, QA); enterprise-grade orchestration |
| **Amazon Kiro** | AWS (Nov 2025) | Requirements → Design → Tasks → Implement → PR; integrated into AWS ecosystem |
| **SDD (Thoughtworks)** | Thoughtworks (2025) | Formalized methodology; Given/When/Then spec structure; BDD-as-few-shot-prompting |

- **Source:** [David Di Spenza — Spec-Driven Development Overview (2025)](https://www.davidedispenza.com/en-us/blog/agent-spec-driven-development)  
  *Why it matters: SciAgent-toolkit can implement a concrete subset of BMAD/spec-kit without reinventing the wheel.*

### SDD Maturity Spectrum

1. **Spec-first development** — specs precede code; code is primary deliverable (most teams)
2. **Spec-anchored development** — governance layers, constitutional constraints, supervision checkpoints (regulated domains)
3. **Spec-as-source development** — specs ARE source code (frontier; few production examples)

For scientific software: **Spec-anchored** is the target — specs serve as audit trail and reproducibility artifact simultaneously.

- **Source:** [BCMS — Spec-Driven Development Definitive Guide (2026)](https://thebcms.com/blog/spec-driven-development)  
  *Why it matters: Defines the appropriate maturity tier for science workflows where auditability matters.*

### Context Window Efficiency

Effective specs achieve **completeness yet conciseness** without enumerating all cases. BDD-style Given/When/Then spec-by-example is "essentially the few-shot prompt technique" — well-crafted examples are more token-efficient than exhaustive prose descriptions.

- **Source:** [Thoughtworks SDD article (2025)](https://www.thoughtworks.com/en-us/insights/blog/agile-engineering-practices/spec-driven-development-unpacking-2025-new-engineering-practices)  
  *Why it matters: Plans stored as AGENTS.md context must be concise; BDD format is a natural fit for science workflows with defined inputs/outputs.*

---

## Area 4 — Reproducibility in Computational Science with AI Agents

### Core Principle: Reproducibility as Infrastructure

**AI-driven scientific software demands constant monitoring, retraining, and governance** unlike classical codes that remain valid for decades once validated. Opaque AI pipelines risk magnifying the reproducibility crisis unless transparency is embedded from the outset.

Key insight: **Academic-grade MLOps differs fundamentally from industrial practice.** Research labs benefit from prioritizing "recomputability, transparent provenance, low-cost automation, and knowledge transfer across student turnover" rather than replicating enterprise-scale continuous deployment.

- **Source:** [Frontiers in Physics — Scientific Software in the AI Era (Nov 2025)](https://www.frontiersin.org/journals/physics/articles/10.3389/fphy.2025.1711356/full)  
  *Why it matters: Validates the "minimum viable MLOps" framing appropriate for academic science labs.*

### Research Compendium Structure

The **research compendium** (Gentleman & Lang, 2007; Marwick et al., 2018) is the canonical unit for organizing reproducible computational research:

```
project-root/
├── data/
│   ├── raw/          # Immutable; protected by cryptographic checksums
│   └── processed/    # All transformations documented; checksums tracked
├── code/ (or scripts/ or analysis/)
│   ├── pipeline/     # Declarative workflow (Snakemake/Nextflow DAG)
│   └── figures/      # Scripts that generate each figure
├── results/
│   ├── figures/      # Output figures; adjacent to generating scripts
│   └── tables/       # Source tables; adjacent to figure that uses them
├── docs/
│   └── notebooks/    # Executable narrative (Jupyter); parameterized; clean kernels enforced
├── AGENTS.md         # Agent conventions for this compendium
└── spec.md           # Analysis plan / hypothesis log
```

**Adjacency rule:** Figures and their generating source scripts and source tables should be co-located or cross-referenced by name. Captions should live in the same directory as the figure or in a companion `.caption.md` / `.caption.json` file.

- **Source:** [How to Read a Research Compendium — arXiv:1806.09525](https://arxiv.org/pdf/1806.09525)  
  *Why it matters: This is the established academic pattern for results organization; SciAgent-toolkit should enforce this structure.*

### FAIR Principles for Computational Research

FAIR (Findable, Accessible, Interoperable, Reusable) extended to research software and workflows:
- **Findable:** Persistent identifiers (DOIs) for datasets, workflows, and model checkpoints
- **Accessible:** Clear access protocols; no authentication barriers for reviewers
- **Interoperable:** Machine-readable metadata using W3C PROV standards
- **Reusable:** Explicit licensing; datasheets; model cards; provenance logs

For workflows specifically (2025 update): "FAIR for research software" principles apply to both pipeline definitions (e.g., Snakemake/Nextflow files) and workflow management systems.

- **Source:** [An Ecosystem of Services for FAIR Computational Workflows — arXiv:2505.15988 (2025)](https://arxiv.org/pdf/2505.15988)  
  *Why it matters: SciAgent-toolkit should emit FAIR-compliant metadata as part of every analysis run.*

### Provenance Tracking and Lab Notebook Patterns

**Recommended provenance stack:**
- **ReproZip** — captures all dependencies, libraries, config parameters automatically
- **W3C PROV** standards for provenance metadata
- **ProvBook** — provenance-based semantic enrichment of Jupyter notebooks
- **Content-addressed storage** — identical data hashes yield identical results

**Lab notebook patterns for AI-generated analyses:**
- Explicitly version, timestamp, and tag LLM contributions with: model identifier, prompt context, temperature setting, tool calls invoked, and hardware details (STARD-AI checklist)
- Parameterize notebooks for reusability across datasets; enforce clean kernel restart + sequential cell execution before publishing
- Archive intermediate artifacts: processed datasets, feature transformations, model checkpoints

**Anti-patterns to avoid:**
- Ephemeral throwaway scripts (write to `scripts/` not `/tmp`)
- Unreproducible random seeds (always set and log seeds)
- Undocumented model versions (pin model ID + date in every analysis)
- "Hidden technical debt" from undocumented configuration choices

- **Source:** [Frontiers in Physics (Nov 2025)](https://www.frontiersin.org/journals/physics/articles/10.3389/fphy.2025.1711356/full)  
  *Why it matters: The STARD-AI checklist provides a concrete audit template that SciAgent-toolkit should emit at the end of each analysis run.*

### Declarative Workflow Architecture

Express computational pipelines as **directed acyclic graphs (DAGs)** using Snakemake, Nextflow, or equivalent:
- Explicit input/output contracts at each pipeline stage
- Version-controlled workflow definitions
- Each stage coupled with its environment specification
- Periodic (not continuous) retraining cycles for ML components
- Feature stores with explicit, stable schema contracts

- **Source:** [Frontiers in Physics (Nov 2025)](https://www.frontiersin.org/journals/physics/articles/10.3389/fphy.2025.1711356/full) + [PLOS Computational Biology — Ten Simple Rules for RDM (2025)](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1013779)  
  *Why it matters: SciAgent-toolkit should generate Snakemake/Nextflow DAGs as output, not ad-hoc bash scripts.*

### LLM-Generated Code Governance

When AI agents generate scientific code:
1. Tag every LLM-generated file with: model ID, prompt hash, generation timestamp
2. Verify physical correctness (conservation laws, boundary conditions) through tests
3. Track training data sources and inherent biases in generated code
4. Maintain auditable human-machine authorship records
5. Use Model Cards documenting intended uses, performance metrics, and limitations

- **Source:** [Frontiers in Physics (Nov 2025)](https://www.frontiersin.org/journals/physics/articles/10.3389/fphy.2025.1711356/full)  
  *Why it matters: Science requires attributing AI contributions — SciAgent-toolkit should auto-generate these tags.*

### Agentic AI Reproducibility Challenges (2025 Survey)

Key obstacles identified:
- Flawed or incomplete data can propagate errors across agentic pipelines
- Lack of human oversight in autonomous agents compounds errors (critical in bio/chem)
- Literature review phases show significant performance drops in tested systems
- "Questions about system reliability, reproducibility, and ethical governance continue to pose significant hurdles"

Benchmark resources: **CORE-Bench** (arXiv:2409.11363) — computational reproducibility benchmark; **REPRO-BENCH** (ACL 2025) — automated reproducibility assessment for social science papers.

- **Source:** [Agentic AI for Scientific Discovery Survey — arXiv:2503.08979 (March 2025)](https://arxiv.org/html/2503.08979v1)  
  *Why it matters: Validates that reproducibility in agentic science is an unsolved problem — SciAgent-toolkit targets a real gap.*

---

## Area 5 — Scientific Figure Best Practice for Dual Legibility

### The Fundamental Problem

A figure designed for Nature print (5–7 pt font, fine line weights at 90 mm width) is **completely illegible when projected on a screen**. Dual legibility requires designing for both simultaneously — or producing two versions.

### Print Journal Specifications (Nature Portfolio — Official)

From the official [Nature Research Figure Guide](https://research-figure-guide.nature.com/figures/preparing-figures-our-specifications/):

| Parameter | Nature Specification |
|-----------|---------------------|
| Panel labels (a, b, c) | 8 pt, bold, lowercase, upright (not italic) |
| Max text size | 7 pt |
| Min text size | 5 pt |
| Font family | Sans-serif only (Arial or Helvetica preferred) |
| Single-column width | 89 mm (3.5 in) |
| 1.5-column width | 120–136 mm |
| Double-column width | 183 mm (7.2 in) |
| Max height | 247 mm (9.7 in) |
| Line art DPI | 1000+ DPI or vector (PDF/EPS preferred) |
| Photos DPI | 300–600 DPI |
| Combination panels | 600 DPI |
| Min line weight | 0.25–0.28 pt |
| Color mode | RGB (journals convert to CMYK for print) |
| Color accessibility | Wong's colorblindness palette; colorblind filter must retain contrast |

- **Source:** [Nature Research Figure Guide (official)](https://research-figure-guide.nature.com/figures/preparing-figures-our-specifications/) + [PlotiVy Nature Guidelines 2025–2026](https://plotivy.app/blog/nature-journal-figure-guidelines-2025)  
  *Why it matters: These are the hard constraints for publication-quality figures; SciAgent-toolkit's figure-generation templates must respect them.*

### Projection / Screen Specifications

Based on Penn State Accessibility Guidelines, Guy Kawasaki's 10/20/30 Rule, and Sidney Smith (1979) visual angle research:

| Room context | Title | Body text | Captions/labels |
|-------------|-------|-----------|-----------------|
| Classroom / small meeting room | 36–40 pt | 24 pt minimum | 18–20 pt |
| Large conference hall | 44+ pt | 28–32 pt | 20–24 pt |
| Online / Zoom presentation | 36 pt | 24–30 pt | 18 pt |

**Scientific basis:** Smith (1979) derived a minimum subtended visual angle of 0.007 radians; at 4 feet distance this requires a letter height of 0.34 in = **25 pt minimum** for posters. Penn State Accessibility: **24 pt absolute floor** for projected slides.

**Line weights for projection:** Significantly heavier than print — minimum 1.5–2 pt to be visible when projected. Fine lines (< 0.5 pt) disappear at projection distances.

**Scientific poster sizes (A0):**
- Title: 85–100 pt bold
- Section headings: 44–56 pt bold
- Body text: 28–36 pt
- Captions: 22–28 pt
- Rule: body text under 20 pt on A0 is unreadable from 1 meter

- **Source:** [ConceptViz — Scientific Poster Font Sizes](https://conceptviz.app/blog/best-fonts-for-scientific-posters-and-figures) + [Andrew Wheeler — Poster Font Size minimum 25pt (2015)](https://andrewpwheeler.com/2015/10/14/poster-presentations-should-have-a-minimum-font-size-of-25-points/) + [Penn State Accessibility — Font Size](https://accessibility.psu.edu/legibility/fontsize/)  
  *Why it matters: SciAgent-toolkit figures need a "dual-mode" design flag: generate print version AND screen version with scaled typography.*

### AIAA Journal Standards (Representative Technical Journal)

| Parameter | AIAA Specification |
|-----------|-------------------|
| Min font size (at column width) | 8 pt |
| Min line weight | 0.28 pt |
| Line art DPI | 600 DPI minimum |
| Photo DPI | 300 DPI minimum |
| Column width | 3.25 in (82.5 mm) |
| Max figure (single column) | 3.25 in × 3.5 in |
| Background screens (with overlaid text) | Max 25% darkness |
| Caption length | 20–25 words maximum |

- **Source:** [AIAA Journal Figures and Tables Guidelines](https://www.aiaa.org/publications/journals/journal-author/guidelines-for-journal-figures-and-tables/)  
  *Why it matters: AIAA provides the most explicit numeric standards for engineering/science journals; a good cross-check against Nature.*

### Decluttering and Accessibility Rules

- **Color:** Never use red-green for data encoding (top rejection reason at Nature and Cell). Use blue-orange, viridis, or cividis. Verify figures pass colorblind simulation filter.
- **Label truncation:** Axis labels and legend text must be complete at final print size; truncation at 5–7 pt is common when figures are scaled down after design.
- **Design at final size:** Set figure dimensions in code before creating the figure (90 mm single / 180 mm double for Nature). Resizing afterward changes text sizes and line thicknesses unpredictably.
- **Test:** Print a test section at full size; stand 1.5 m away; if straining to read, increase size.
- **Font:** Sans-serif universally preferred (Arial, Helvetica for print; Roboto, Inter for screen); no decorative or serif fonts in figure panels.
- **Redundancy:** Encode data with color + shape + label (not color alone).

- **Source:** [SciDraw — Best Fonts for Scientific Figures](https://sci-draw.com/blog/best-fonts-for-scientific-figures-and-posters) + [Nature Research Figure Guide](https://research-figure-guide.nature.com/figures/preparing-figures-our-specifications/)  
  *Why it matters: SciAgent-toolkit figure generator should enforce these rules automatically (color check, font check, size-at-creation).*

---

## Recommendations for SciAgent-toolkit

Mapping external best practice to the five owner pain points:

### Pain Point 1 — Figure Legibility (print journal column vs. projected screen)

**External finding:** Print and screen require fundamentally different specs. Nature mandates 5–7 pt fonts at 89 mm width; projected slides require 24 pt minimum (28–32 pt for large halls). Line weights below 1 pt are invisible when projected.

**Recommendation:**
- SciAgent-toolkit should generate figures in **two modes**: `print` (at exact journal column width with 5–7 pt labels, 0.25+ pt lines, 300–1000 DPI) and `screen` (scaled to 1920×1080 equivalent with 24–28 pt labels, 2 pt+ lines)
- Enforce a **checklist hook** (`PreToolUse` or `Stop`) that verifies: (a) font size ≥ 5 pt at final print size, (b) color palette passes colorblind filter, (c) no label truncation, (d) line weight ≥ 0.25 pt
- Use **vector output by default** (PDF/SVG) — fonts stay sharp at any zoom; rasterize only for final submission
- Font: Arial or Helvetica for print panels; Roboto/Inter as acceptable screen alternatives

### Pain Point 2 — Results Placement and Adjacency (figures + source tables + captions co-located)

**External finding:** Research compendium best practice (Marwick et al., 2018) mandates adjacency between figures, their generating scripts, and source tables. PLOS Computational Biology guidelines require figures to be citable by number, with captions placed adjacent.

**Recommendation:**
- Enforce a **canonical results directory structure** in AGENTS.md:
  ```
  results/
  ├── fig01/
  │   ├── fig01.pdf               # Vector figure
  │   ├── fig01_source_data.csv   # Source table
  │   ├── fig01_caption.md        # Caption text (citable)
  │   └── fig01_generate.py       # Generating script
  └── fig02/
      └── ...
  ```
- A `Stop` hook should **refuse to save a figure** unless `fig_source_data` and `fig_caption.md` are present in the same directory
- Include a `FAIR-manifest.json` at the results root listing each figure with its generating script hash, data hash, and model ID used

### Pain Point 3 — README / AGENTS.md Upkeep (conventions drift without maintenance)

**External finding:** AGENTS.md files over 500 lines show degraded enforcement. The most effective AGENTS.md files are under 150 lines, imperative, specific, and updated in the same PR that introduces new conventions. Assign an owner; conduct quarterly reviews.

**Recommendation:**
- SciAgent-toolkit's `AGENTS.md` should be under 150 lines covering: setup commands, project structure, prohibited patterns (no ephemeral scripts, no hardcoded paths), figure conventions, and provenance requirements
- Use the **referential architecture**: CLAUDE.md one-liner pointing to AGENTS.md; Cursor `.mdc` rules for scoped activation only
- Add a **quarterly-review hook** via a `CronCreate`-triggered reminder in `.claude/hooks/quarterly_review.sh`
- Encode science-specific skills (e.g., `/new-analysis`, `/finalize-figure`, `/write-caption`) as SKILL.md files — these stay current because they're invoked daily

### Pain Point 4 — Planning Decomposition (breaking analysis into context-window-sized phases)

**External finding:** Harper Reed's workflow (2025) targets 8–12 steps per plan, each implementable in one agentic context with tests, stored in `prompt_plan.md`. BMAD structures this at enterprise scale with named personas. Thoughtworks formalizes it as SDD with BDD-style specs.

**Recommendation:**
- Adopt **Harper Reed's artifact pattern** as the SciAgent-toolkit planning standard:
  - `spec.md` — analysis hypothesis, data inputs, expected outputs, success criteria (Given/When/Then format)
  - `prompt_plan.md` — 8–15 numbered steps, each a self-contained executable prompt with test criterion
  - `todo.md` — checkboxes generated from prompt_plan.md; committed to repo for traceability
- Each step should produce a **testable intermediate artifact** (e.g., a validated intermediate figure, a passing unit test on the data transformation)
- Use a `plan` skill (`/plan`) that: (a) asks clarifying questions one at a time, (b) outputs spec.md + prompt_plan.md, (c) requires human approval before execution begins
- Encode the **spec-anchored maturity tier** (not just spec-first): specs serve simultaneously as planning artifact AND reproducibility audit trail

### Pain Point 5 — Reproducibility / No Ephemeral Scripts

**External finding:** Frontiers in Physics (2025) and the research compendium pattern both require: version-controlled workflow definitions, content-addressed storage, pinned dependencies, logged random seeds, LLM contribution tagging (model ID + prompt hash + timestamp), and W3C PROV provenance records.

**Recommendation:**
- **Block ephemeral script patterns** via a `PreToolUse` hook: reject writes to `/tmp`, `~/scratch`, or any path outside the project directory structure
- **Auto-tag LLM-generated files** with a header comment: `# Generated by: claude-sonnet-4-6 | 2026-06-23 | prompt_hash: abc123`
- Emit a **STARD-AI provenance record** at the end of every analysis run (via `Stop` hook): dataset version, model ID, random seeds, tool calls invoked, hardware spec
- Use **Snakemake or Nextflow** DAG definitions (generated as a skill output, never as ad-hoc bash) for all pipeline steps
- Enforce **clean environment policy** via AGENTS.md: all dependencies pinned in `environment.yml` / `pyproject.toml`; Docker/Singularity container required for publication-grade runs
- The `spec.md` doubles as the lab notebook entry: it is committed before any code is written, updated after the run with actual outcomes, and never deleted

---

## Top Source URLs

1. [Linux Foundation — AAIF Formation (2025-12-09)](https://www.linuxfoundation.org/press/linux-foundation-announces-the-formation-of-the-agentic-ai-foundation)
2. [Anthropic — Steering Claude Code: skills, hooks, subagents (official blog)](https://claude.com/blog/steering-claude-code-skills-hooks-rules-subagents-and-more)
3. [Anthropic — Create Custom Subagents (Claude Code Docs)](https://code.claude.com/docs/en/sub-agents)
4. [Harper Reed — My LLM Codegen Workflow (2025-02-16)](https://harper.blog/2025/02/16/my-llm-codegen-workflow-atm/)
5. [Thoughtworks — Spec-Driven Development (2025)](https://www.thoughtworks.com/en-us/insights/blog/agile-engineering-practices/spec-driven-development-unpacking-2025-new-engineering-practices)
6. [Frontiers in Physics — Scientific Software in the AI Era (Nov 2025)](https://www.frontiersin.org/journals/physics/articles/10.3389/fphy.2025.1711356/full)
7. [Nature Research Figure Guide (official)](https://research-figure-guide.nature.com/figures/preparing-figures-our-specifications/)
8. [DeployHQ — AI Coding Config Files Guide](https://www.deployhq.com/blog/ai-coding-config-files-guide)
9. [Codex KB — AGENTS.md vs CLAUDE.md Cross-Tool Portability (2026-05-27)](https://codex.danielvaughan.com/2026/05/27/agent-instruction-files-agents-md-claude-md-cross-tool-portability-codex-cli/)
10. [BuildBetter — AGENTS.md Complete Guide (2026)](https://blog.buildbetter.ai/agents-md-complete-guide-for-engineering-teams-in-2026/)
11. [Penn State Accessibility — Font Size Guidelines](https://accessibility.psu.edu/legibility/fontsize/)
12. [AIAA Journal Figure Guidelines](https://www.aiaa.org/publications/journals/journal-author/guidelines-for-journal-figures-and-tables/)
13. [arXiv:1806.09525 — How to Read a Research Compendium (Marwick et al.)](https://arxiv.org/pdf/1806.09525)
14. [arXiv:2503.08979 — Agentic AI for Scientific Discovery Survey (March 2025)](https://arxiv.org/html/2503.08979v1)
15. [AGENSI.io — SKILL.md Open Standard](https://www.agensi.io/learn/agent-skills-open-standard)
