# Prior-Art Scout: mattpc/skills — Transferable Patterns for SciAgent-toolkit

**Date:** 2026-06-23  
**Source:** `/data1/users/antonz/pipeline/.ref/mattpc/skills/skills`  
**Purpose:** Identify transferable patterns for reproducible, agentic scientific computing — focusing on five owner concerns: figure quality, results layout, documentation upkeep, planning/orchestration, and reproducibility/persistence.

---

## 1. Inventory Table

| Skill | Path (relative to `skills/`) | Rel. | Concern # | Transferable Idea |
|---|---|---|---|---|
| `writing-great-skills` | `productivity/writing-great-skills/SKILL.md` | 3 | 4, 5 | Skill-authoring meta-discipline: predictability, completion criteria, progressive disclosure, pruning — directly governs how SciAgent skills should be written |
| `decision-mapping` | `in-progress/decision-mapping/SKILL.md` | 3 | 4 | Stateful, git-tracked decision map with typed tickets (Research/Prototype/Discuss); fog-of-war metaphor for multi-session planning |
| `prototype` | `engineering/prototype/SKILL.md` + `LOGIC.md` | 3 | 4, 5 | Throwaway code discipline: capture the *answer* in a durable note (NOTES.md / ADR), delete the shell; logic module stays portable |
| `diagnosing-bugs` | `engineering/diagnosing-bugs/SKILL.md` | 3 | 5 | Phase-gated loop with explicit "tight feedback loop" before hypothesis; completion checklists (tagged debug logs, throwaway cleanup); HITL fallback script |
| `triage` + `AGENT-BRIEF.md` | `engineering/triage/` | 3 | 4 | Behavioral (not procedural) spec writing; durable agent briefs that avoid file paths/line numbers; acceptance-criteria checklists; out-of-scope KB |
| `ask-matt` | `engineering/ask-matt/SKILL.md` | 3 | 4 | Router/meta-skill over all user-invoked skills; explicit "main flow: idea→ship" with phase discipline and context-hygiene rules |
| `to-prd` | `engineering/to-prd/SKILL.md` | 3 | 4 | PRD template with seam-first thinking; avoids file paths; recommends prototype snippets for decisions that can't be prose-captured |
| `to-issues` | `engineering/to-issues/SKILL.md` | 3 | 4 | Tracer-bullet vertical slices; dependency-ordered publishing; acceptance-criteria checklists per issue |
| `review` | `in-progress/review/SKILL.md` | 3 | 4 | Parallel sub-agents for two-axis review (Standards vs Spec); aggregation without reranking; fixed-point pinning before spawning |
| `implement` | `engineering/implement/SKILL.md` | 2 | 4 | Minimal skill: implement → tdd → review → commit in one compact chain; end-to-end phase discipline |
| `domain-modeling` | `engineering/domain-modeling/SKILL.md` | 2 | 3, 4 | Inline glossary updates as decisions crystallize; ADR offered only for hard-to-reverse, surprising, real trade-offs; CONTEXT.md stays devoid of implementation |
| `improve-codebase-architecture` + `HTML-REPORT.md` | `engineering/improve-codebase-architecture/` | 2 | 1, 2 | Visual HTML reports saved to temp dir (not repo); before/after diagrams with explicit style guide (editorial, not corporate); recommendation-strength badges |
| `handoff` | `productivity/handoff/SKILL.md` | 2 | 5 | Compact conversation → handoff document (saved to OS temp dir, not workspace); avoids duplicating content in artifacts; includes "suggested skills" for next session |
| `codebase-design` | `engineering/codebase-design/SKILL.md` | 2 | 4 | Shared vocabulary discipline; "design it twice" parallel sub-agents; deletion test; single-source-of-truth glossary |
| `tdd` | `engineering/tdd/SKILL.md` | 2 | 4, 5 | Tracer-bullet (vertical) not horizontal slices; behavior-only tests through public interfaces; per-cycle checklist |
| `teach` | `productivity/teach/SKILL.md` | 2 | 2, 3, 5 | Structured output workspace (MISSION.md, learning-records/, lessons/, assets/, reference/); learning records analogized to ADRs; lessons as beautiful, citable HTML |
| `setup-matt-pocock-skills` | `engineering/setup-matt-pocock-skills/SKILL.md` | 2 | 3, 4 | Per-repo configuration skill that writes `docs/agents/` files consumed by other skills; prompt-driven, not scripted |
| `grill-with-docs` / `grilling` | `engineering/grill-with-docs/`, `productivity/grilling/` | 1 | 4 | Relentless one-question-at-a-time interview; recommended answer per question; codebase exploration as alternative to asking |
| `prototype/LOGIC.md` | `engineering/prototype/LOGIC.md` | 1 | 5 | Logic-vs-UI prototype distinction; portable pure reducer / state machine stays; TUI shell is throwaway; NOTES.md captures answer |
| `git-guardrails-claude-code` | `misc/git-guardrails-claude-code/SKILL.md` | 1 | 5 | PreToolUse hook pattern for blocking destructive operations; bundled script + settings.json registration |
| `obsidian-vault` | `personal/obsidian-vault/SKILL.md` | 1 | 5 | Flat vault structure with index notes; search/link patterns; unit-of-learning philosophy |
| `handoff` | `productivity/handoff/SKILL.md` | 2 | 5 | Cross-session continuity via temp-dir document; reference artifacts by path, don't duplicate |
| `deprecated/design-an-interface` | `deprecated/design-an-interface/SKILL.md` | 1 | 4 | Superseded by `codebase-design/DESIGN-IT-TWICE.md`; parallel sub-agent interface exploration |

**Relevance scale:** 0 = irrelevant, 1 = low, 2 = medium, 3 = high

---

## 2. Skill-Format Craft

### 2.1 Frontmatter Fields

```yaml
---
name: skill-name                  # canonical invocation name
description: "..."                # MODEL-FACING trigger; omit → user-invoked only
disable-model-invocation: true    # strips description from agent context
argument-hint: "What to learn?"   # hint shown when user types the skill name
---
```

The critical distinction is `disable-model-invocation: true` — it creates a **user-invoked** skill that pays zero context load. Skills like `handoff`, `implement`, `to-prd`, and `to-issues` are all user-invoked: they require deliberate human initiation. Model-invoked skills (`diagnosing-bugs`, `review`, `grilling`, `domain-modeling`) include rich descriptions with trigger phrasing so the agent fires them autonomously.

**Key insight:** The description field does double duty — it IS the invocation trigger. Every word is context load. Write descriptions for the agent, not the human.

### 2.2 Trigger Phrasing in Descriptions

Good examples from mattpc:

- `diagnosing-bugs`: "Use when the user says 'diagnose'/'debug this', or reports something broken/throwing/failing/slow."
- `review`: "Use when the user wants to review a branch, a PR, work-in-progress changes, or asks to 'review since X'."
- `grilling`: "Use when the user wants to stress-test a plan before building, or uses any 'grill' trigger phrases."

Pattern: **enumerate surface forms** ("diagnose"/"debug this"/reports broken/failing/slow) not just paraphrases. This broadens invocation reliability.

### 2.3 Progressive Disclosure via Companion Files

The `SKILL.md` stays lean; heavy reference is pushed into sibling files reached by a pointer in the text:

```
engineering/diagnosing-bugs/
  SKILL.md                         ← steps only, pointer to script
  scripts/hitl-loop.template.sh    ← external asset
```

```
productivity/writing-great-skills/
  SKILL.md                         ← principles, points to GLOSSARY
  GLOSSARY.md                      ← full definitions
```

```
engineering/codebase-design/
  SKILL.md                         ← core vocab
  DEEPENING.md                     ← dependency categories
  DESIGN-IT-TWICE.md               ← sub-agent pattern
```

The wording of the pointer matters: "see [DEEPENING.md](DEEPENING.md): dependency categories, seam discipline, and replace-don't-layer testing" — specific enough that the agent knows *when* to follow it.

### 2.4 Completion Criteria as First-Class Citizens

Every step ends on an explicit, checkable completion criterion. From `diagnosing-bugs` Phase 1:

> Phase 1 is done when the loop is **tight** and **red-capable**: you can name **one command** — a script path, a test invocation, a curl — that you have **already run at least once** (paste the invocation and its output)

The criterion is:
- [ ] Red-capable (asserts user's exact symptom)
- [ ] Deterministic
- [ ] Fast (seconds)
- [ ] Agent-runnable

This is the strongest pattern in the library: **checkboxes embedded in the completion criterion**, not just at the end of a skill.

### 2.5 Leading Words (Leitwort)

From `writing-great-skills`: a "leading word" is a compact pretrained concept recruited to anchor behavior — _tight_ (loop), _red_ (bug), _tracer bullet_, _fog of war_. It accumulates definition by repetition and recruits priors the model already holds. Each skill has 1–3 leading words that appear multiple times.

In practice: `diagnosing-bugs` uses "tight" throughout; `decision-mapping` uses "fog of war"; `to-issues` uses "tracer bullet". These words do real behavioral work with minimal tokens.

### 2.6 Router Skills (Meta-Skills)

`ask-matt` is a pure routing skill — it describes the entire workflow graph (main flow + on-ramps + standalone skills) so users have one thing to remember. It uses flowchart-style prose to describe the "idea → ship" path. This is the pattern for a "meta-skill" that guides users to the right tool.

### 2.7 Scope Discipline

`AGENT-BRIEF.md` defines a template with an explicit **"Out of scope"** section. `to-issues` and `to-prd` both carry Out-of-Scope sections. This is a consistent pattern: every deliverable-generating skill explicitly names what will NOT be done.

### 2.8 Behavioral (not Procedural) Specs

From `AGENT-BRIEF.md`:
- **Good:** "The `SkillConfig` type should accept an optional `schedule` field of type `CronExpression`"
- **Bad:** "Open src/types/skill.ts and add a schedule field on line 42"
- Never reference file paths or line numbers in specs — they go stale.

### 2.9 Artifact Placement Discipline

`improve-codebase-architecture` and `handoff` both explicitly place outputs **outside the repo**:
- HTML reports → `$TMPDIR/architecture-review-<timestamp>.html`
- Handoff documents → "temporary directory of the user's OS — not the current workspace"

This is a direct analog to results-layout discipline: intentional choices about where artifacts land.

### 2.10 Single-Source-of-Truth Vocabulary

`domain-modeling` maintains `CONTEXT.md` as a pure glossary (no implementation details), updated inline as decisions crystallize. ADRs in `docs/adr/` record only hard-to-reverse, surprising, real trade-off decisions. Multiple skills (`tdd`, `diagnosing-bugs`, `improve-codebase-architecture`) all open by reading `CONTEXT.md` — shared vocabulary as a first-class architectural artifact.

---

## 3. Top Transferable Patterns

### Pattern A — Decision Mapping for Multi-Session Science Plans (Concern #4)

**Source:** `in-progress/decision-mapping/SKILL.md`

A stateful, git-tracked markdown file with typed tickets (Research / Prototype / Discuss), each sized to one agent session. Uses "fog of war" as leading word — the map is *deliberately* incomplete beyond the frontier.

> "The decision map is a single compact Markdown file, one per planning effort, git-tracked alongside the project. It is the canonical artifact — the whole map is loaded as context into every session, so it must stay compact."

**Application to SciAgent:** Directly adoptable for multi-session analysis plans. Research tickets map to literature/data exploration phases; Prototype tickets map to exploratory analysis scripts; Discuss tickets map to interpretation sessions with the PI. The "fog of war" metaphor matches real scientific uncertainty. The per-ticket size constraint (100K token session) is directly portable.

---

### Pattern B — Durable Agent Briefs / Behavioral Specs (Concern #4, #5)

**Source:** `engineering/triage/AGENT-BRIEF.md`

> "Durability over precision. The issue may sit in `ready-for-agent` for days or weeks. The codebase will change in the meantime. Write the brief so it stays useful even as files are renamed, moved, or refactored."

Rules:
- Describe interfaces and behavioral contracts, not file paths
- Acceptance criteria as independently-verifiable checkboxes
- Explicit out-of-scope section

**Application to SciAgent:** Science tasks delegated to subagents need the same durability discipline. A brief for "run differential expression analysis on cohort X" should describe data contracts and expected output shapes, not specific file paths that will change as the pipeline evolves.

---

### Pattern C — Prototype Answer Capture (Concern #5 — Reasoning Trace Persistence)

**Source:** `engineering/prototype/SKILL.md` + `LOGIC.md`

> "The answer is the only thing worth keeping from a prototype. Capture it somewhere durable (commit message, ADR, issue, or a NOTES.md next to the prototype) along with the question it was answering."

The logic module (pure reducer / state machine) is kept portable and potentially absorbed into production; the TUI shell is deleted. The throwaway is clearly marked from day one.

**Application to SciAgent:** Exploratory analysis scripts are the scientific equivalent of prototypes. The pattern: (1) state the question explicitly before writing code, (2) isolate the core logic in a portable module, (3) capture the answer in a NOTES.md or lab-notebook entry, (4) delete or absorb the scaffolding. This is the anti-ephemeral-scripts discipline the owner wants.

---

### Pattern D — Phase-Gated Feedback Loop Before Hypothesis (Concern #5)

**Source:** `engineering/diagnosing-bugs/SKILL.md`

> "Phase 1 is done when the loop is tight and red-capable: you can name one command — a script path, a test invocation, a curl — that you have already run at least once."

> "If you catch yourself reading code to build a theory before this command exists, stop — jumping straight to a hypothesis is the exact failure this skill prevents."

The completion criterion is a checklist with four binary properties. Debug logs are tagged with unique prefixes for deterministic cleanup. Throwaway harnesses must be deleted in Phase 6 cleanup.

**Application to SciAgent:** Replace "bug" with "unexpected result" and the structure transfers directly. Before hypothesizing about a failed pipeline run, require a reproducible minimal command that demonstrates the failure. The tagged-log / cleanup discipline maps to the no-ephemeral-scripts concern.

---

### Pattern E — Parallel Sub-Agent Review with Fixed Axes (Concern #4)

**Source:** `in-progress/review/SKILL.md`

> "Two-axis review... Both axes run as parallel sub-agents so they don't pollute each other's context, then this skill aggregates their findings."

Axes: Standards (does it follow repo conventions?) vs Spec (does it implement what was asked?). The fixed-point pinning step happens before spawning. Aggregation preserves both reports separately — no reranking across axes.

**Application to SciAgent:** Science outputs need a two-axis review: (1) Methods correctness vs documented analysis plan, (2) Presentation quality vs figure/reporting conventions. Spawning parallel subagents for each axis prevents one from contaminating the other's judgment.

---

### Pattern F — Shared Vocabulary as First-Class Artifact (Concern #3)

**Source:** `engineering/domain-modeling/SKILL.md` + multiple consumers

Every skill that touches the codebase opens by reading `CONTEXT.md`. The domain-modeling skill updates it *inline as decisions crystallize*, not in a batch afterward. ADRs are offered sparingly — only when all three hold: hard to reverse, surprising without context, result of a real trade-off.

> "`CONTEXT.md` should be totally devoid of implementation details. Do not treat `CONTEXT.md` as a spec, a scratch pad, or a repository for implementation decisions. It is a glossary and nothing else."

**Application to SciAgent:** A `CONTEXT.md` (or equivalent `GLOSSARY.md`) per analysis project, maintained by a domain-modeling skill, consumed by every other skill. The ADR pattern maps well to pipeline design decisions (e.g., "why we use Wilcoxon not t-test for this cohort").

---

### Pattern G — Teach Workspace as Lab Notebook (Concern #2, #5)

**Source:** `productivity/teach/SKILL.md`

The teach workspace has a prescribed artifact layout:
```
workspace/
  MISSION.md                 ← grounding document
  RESOURCES.md               ← external sources
  NOTES.md                   ← preferences / working notes
  reference/*.html           ← compressed cheat sheets
  learning-records/*.md      ← ADR-style key insights, titled 0001-*.md
  lessons/*.html             ← primary output, beautiful and citable
  assets/                    ← reusable components shared across lessons
```

> "Learning records capture non-obvious lessons and key insights that may need to be revised later... They are loosely equivalent to architectural decision records in software development."

**Application to SciAgent:** The `learning-records/` → `analysis-records/` pattern maps directly to lab notebook entries. The numbered naming scheme (`0001-*.md`) and ADR-style structure (question / decision / reasoning / alternatives) is directly adoptable for documenting analysis choices. The `assets/` component reuse pattern maps to shared figure stylesheets.

---

### Pattern H — Artifact Placement Outside Repo (Concern #2)

**Source:** `engineering/improve-codebase-architecture/SKILL.md`

> "Write a self-contained HTML file to the OS temp directory so nothing lands in the repo. Resolve the temp dir from `$TMPDIR`, falling back to `/tmp`."

`handoff` uses the same rule: "Save to the temporary directory of the user's OS — not the current workspace."

**Application to SciAgent:** Exploratory / intermediate artifacts (HTML reports, diagnostic plots, scratch summaries) belong in temp directories. Only reviewed, named outputs get committed to the project. This is the output-placement discipline: intentional, not accidental.

---

## 4. Skill-Authoring Craft Worth Adopting

1. **Model-invoked vs user-invoked split.** Most SciAgent skills are probably user-invoked (deliberate, typed by researcher). Only skills the agent should fire autonomously need descriptions. Using `disable-model-invocation: true` for the majority reduces context load dramatically.

2. **Completion criteria as checklists.** Every step should end with a checkable criterion — not "analysis is done" but a named command that was run and its output captured.

3. **Progressive disclosure.** Keep `SKILL.md` to steps and essential reference; push heavy methodology, templates, and examples into sibling files (`METHODS.md`, `REPORT-FORMAT.md`, `FIGURE-CONVENTIONS.md`) reachable by named context pointers.

4. **Leading words.** Pick 1–3 pretrained concepts per skill and repeat them as tokens: _tight_ (feedback loop), _durable_ (spec), _tracer_ (slice), _fog_ (unknown territory). These anchor behavior with minimal tokens.

5. **Single source of truth for vocabulary.** One `CONTEXT.md` per analysis domain, updated inline as terms crystallize, consumed by every other skill. The glossary is the architecture.

6. **Router skill.** A single user-invoked `ask-sciagent` (or equivalent) that describes the entire workflow — which skill for which situation — means users only need to remember one skill.

7. **Pruning discipline.** Hunt sediment (stale layers), duplication (same meaning twice), and no-ops (instructions the model already follows). Skills should shrink over time, not grow.

8. **Scope discipline.** Every deliverable-generating skill should have an explicit "Out of scope" section. This prevents gold-plating and scope creep in agentic runs.
