> **Updated** 2026-05-24 with grounded findings from local repo at `docs/.ref/science-skills/`.

# 04 — Anthropic Skills, `skill-creator`, and the "Antigravity Skills" attribution

## Negative finding first: there is no Anthropic "Antigravity Skills" paper

The spec's executive summary and ADRs 005 and 006 reference "the Antigravity Skills paper" as if it were an Anthropic publication describing a three-tier eval pattern, a `workflow-skill-creator` skill, and biology-domain numbers (49% → 93% internal reliability, BioReason VEP 41% → 61%). This conflates two different things:

1. **Google Antigravity** — a Google product (announced Nov 2025; built on Gemini) that supports an "Agent Skills" extension surface. See [How to Build Custom Skills in Google Antigravity](https://medium.com/google-cloud/tutorial-getting-started-with-antigravity-skills-864041811e0d) and [Authoring Google Antigravity Skills codelab](https://codelabs.developers.google.com/getting-started-with-antigravity-skills). This is the source of the term "Antigravity Skills."
2. **Anthropic Agent Skills** — Anthropic's SKILL.md-based skill system, documented at [Agent Skills — Claude API docs](https://platform.claude.com/docs/en/agents-and-tools/agent-skills/overview), with reference skills at [github.com/anthropics/skills](https://github.com/anthropics/skills).

No paper titled "Antigravity Skills" by Anthropic surfaces in web search or on arXiv. The biology numbers cited (49% → 93%, 41% → 61% on BioReason VEP) are not traceable to any Anthropic document. **Treat the spec's claim as unverified for Anthropic.** Possible upstream: the Google DeepMind science-skills technical report (see below) — but it ships from Google, not Anthropic, and the spec doesn't cite it either.

The closest published Anthropic artifact is [The Complete Guide to Building Skills for Claude](https://resources.anthropic.com/hubfs/The-Complete-Guide-to-Building-Skill-for-Claude.pdf), which is a how-to guide, not a benchmark paper.

## What `skill-creator` actually does (Anthropic's, not "workflow-skill-creator")

There is no skill called `workflow-skill-creator` in any Anthropic-owned repository. The skill that exists in `anthropics/skills` is [`skill-creator`](https://github.com/anthropics/skills/tree/main/skill-creator), described in [SKILL.md](https://github.com/anthropics/skills/blob/main/skills/skill-creator/SKILL.md).

Anthropic's `skill-creator` is an **interactive, interview-driven** skill — not a trace distiller. Its workflow:

1. **Capture intent** — ask clarifying questions about what the skill should do, when it should trigger, output format, test-case needs.
2. **Interview & research** — probe edge cases, IO formats, dependencies.
3. **Write SKILL.md** — generate frontmatter (name, description) + Markdown instructions.
4. **Test & evaluate** — spawn parallel subagents for with-skill vs baseline runs against test prompts; produce `benchmark.json` and `grading.json`.
5. **Improve** — iterate.
6. **Optimize description** — fine-tune the description string for trigger accuracy.
7. **Package** — emit a `.skill` file.

Inputs: user intent (conversational), existing skill code (if improving), 2–3 realistic test prompts, user feedback via HTML review UI.

Outputs: `SKILL.md`, `evals/evals.json`, per-iteration workspace directory with side-by-side benchmark artefacts.

**Anthropic's `skill-creator` does not ingest a trace or session — it interviews the user.** The trace-distiller framing in the sciagent spec has no analog in `anthropics/skills`.

## A `workflow-skill-creator` *does* exist — in Google DeepMind's science-skills repo

The local reference at `docs/.ref/science-skills/skills/workflow_skill_creator/SKILL.md` ships a skill named `workflow-skill-creator` with the description:

> "Distills a completed user workflow or interaction into a reusable agent skill. Use when the user asks to turn their workflow, interaction, or multi-step process into a skill, or when they say 'make this a skill', 'create a skill from what we just did', 'package this workflow' or similar. Do not use for creating skills from scratch without an existing workflow (use a generic skill-creator for that)."
> — `docs/.ref/science-skills/skills/workflow_skill_creator/SKILL.md:3-9`

This is the source of the spec's `workflow-skill-creator` reference — not an Anthropic skill. Ownership is Google DeepMind; bundle ships as the "Science" plugin in Google Antigravity (`docs/.ref/science-skills/README.md:29-42`).

What it actually does — read in full from the local copy — is **not** trace distillation in the sense the sciagent spec implies. It is a four-phase **interview-driven authoring workflow** that uses the user's *recently completed conversation* as the seed context, but does not consume any structured session log:

- **Phase 1 (Mandatory): Brainstorming** — iterative 2-3-questions-per-round Q&A across five rounds (workflow shape, flexibility, dependencies, scope/code-or-not, sample query) — `:23-112`.
- **Phase 2: Skill Design** — produce a design document, wait for explicit user approval — `:114-127`.
- **Phase 3: Implementation** — six numbered rules: reuse existing skills, rate-limit new APIs, default to CLI-script pattern (`references/cli_script_template.py`), default to file output, instruction-only fallback, fixed SKILL.md structure — `:129-259`.
- **Phase 4: Validation** — manual test invocation; if sample query provided, run it — `:261-270`.

The only "reference" file shipped alongside is `cli_script_template.py` — no JSONL parser, no trace consumer, no `--from-trace` flag, no eval harness. The skill works by relying on the agent's own recent conversation memory ("Here's my understanding of the workflow: [summary]. Is this accurate?" — `:32-34`), not on parsing a recorded session file.

**Conclusion: a `workflow-skill-creator` exists, but it does the same thing Anthropic's `skill-creator` does — interview the user — with the framing that an interaction has just happened.** Neither shipped skill is what ADR-005 proposes. ADR-005's trace-ingesting distiller is still net-new. The negative finding sharpens: *no published skill in either Anthropic or Google's reference bundles ingests a session JSONL and produces a SKILL.md draft.* The sciagent design has no shipped precedent.

What ADR-005 *could* borrow from the Google `workflow-skill-creator`:
- The **four-phase gate model** (brainstorm → design → implement → validate), with explicit "wait for user approval" between phases. Most distillers fail by going implementation-first.
- The **brainstorming completion checklist** (`:99-112`) — nine items the distiller must answer before writing any code. Useful as a sanity contract for the autonomous distiller too.
- The **fixed SKILL.md structure** template (`:226-259`) — a generated SKILL.md must have Overview / Dependencies / Quick Start / Utility Scripts / Workflow / Rate Limiting / Common Mistakes sections in order.
- The **CLI-pattern-by-default** rule (`:178-205`) — multi-subcommand argparse with `--output` and required arguments. Avoids the "summarized as prose, not as runnable code" failure mode.
- The **rate-limit lookup discipline** (`:155-176`) — search public docs first, default to 1 rps, use `time.monotonic()`, file-lock for cross-process safety.

None of these require code-level porting; they are authoring conventions worth encoding as part of ADR-005's distiller prompt.

## The "three-tier eval pattern" claim — confirmed absent in both repos

The local `docs/.ref/science-skills/` ships **no eval harness, no benchmark numbers, no three-tier pattern documentation, no `evals/` directory, no scoring infrastructure**. `CONTRIBUTING.md` (`docs/.ref/science-skills/CONTRIBUTING.md:7`) says external contributions are not accepted — issue-only reporting. The repo is a delivery vehicle for skills, not a methodology document.

There is a [linked Google DeepMind technical report](https://storage.googleapis.com/deepmind-media/papers/google_deepmind_science_skills_for_antigravity_towards_efficient_and_reliable_scientific_workflows.pdf) referenced from the README (`docs/.ref/science-skills/README.md:62-64`) that may describe evaluations — not pulled into the local ref. If the spec's "49% → 93%" and "41% → 61%" numbers exist anywhere, they would be in that PDF. The repo itself is silent. **The three-tier eval pattern is still not citable from anywhere in the local reference material.**

**Recommendation unchanged:** rewrite ADR-006's context to drop the unverifiable attribution and own the pattern as a sciagent-original design. The pattern is defensible on its own merits. If anyone wants to map it to the DeepMind tech report, that's a citation hunt, not a precedent claim.

## What does exist in the Anthropic skills ecosystem

[`anthropics/skills`](https://github.com/anthropics/skills) ships:
- `skills/skill-creator/` — interview-driven (covered above)
- Document skills: `pdf`, `docx`, `pptx`, `xlsx`
- A `spec/` directory with the Agent Skills spec
- A `template/` skill scaffold

None of these are biology-domain skill bundles. The financial-services-plugins, claude-plugins-official, and claude-code repos contain their own `skill-creator` mirrors.

[`agentskills.io`](https://agentskills.io/specification) is the standalone SKILL.md spec — covered in file 05. It is the official format both Anthropic, Pi, *and Google Antigravity science-skills* consume — confirmed by spot-checking `docs/.ref/science-skills/skills/*/SKILL.md` frontmatter, which matches the spec exactly (`name`, `description`, optional `>-` folding for multiline) and uses no `metadata.*` block at all. Anthropic and Google both follow the minimal-frontmatter style; sciagent's metadata-heavy style remains the outlier.

## Bottom line

- The `workflow-skill-creator` shipped in `docs/.ref/science-skills/` is **interview-driven, not trace-driven**. ADR-005's trace-ingesting variant is genuinely net-new. Reframe it as such — and cite the DeepMind skill as the closest precedent for the *four-phase gate model*, not for trace ingestion.
- The three-tier eval pattern in ADR-006 is **not present in the local Google reference repo**. Possibly in the linked DeepMind PDF; until cited, treat as sciagent-original.
- The benchmark numbers cited (49% → 93%, 41% → 61%) are **still unsourced from any locally-available material**. Drop them or cite the DeepMind PDF explicitly before they ship in any methods paper.
- ADR-001's reference to "the agentskills.io spec" is the only well-attributed claim in this neighborhood — see file 05. Google's science-skills repo confirms the spec is the de-facto cross-vendor convention.

Sources:
- [Agent Skills — Claude API docs](https://platform.claude.com/docs/en/agents-and-tools/agent-skills/overview)
- [anthropics/skills repo](https://github.com/anthropics/skills)
- [skill-creator SKILL.md](https://github.com/anthropics/skills/blob/main/skills/skill-creator/SKILL.md)
- [The Complete Guide to Building Skills for Claude (PDF)](https://resources.anthropic.com/hubfs/The-Complete-Guide-to-Building-Skill-for-Claude.pdf)
- `docs/.ref/science-skills/README.md`
- `docs/.ref/science-skills/plugin.json`
- `docs/.ref/science-skills/skills/workflow_skill_creator/SKILL.md`
- `docs/.ref/science-skills/skills/workflow_skill_creator/references/cli_script_template.py`
- `docs/.ref/science-skills/CONTRIBUTING.md`
- [DeepMind Science Skills technical report (PDF)](https://storage.googleapis.com/deepmind-media/papers/google_deepmind_science_skills_for_antigravity_towards_efficient_and_reliable_scientific_workflows.pdf) — referenced from README but not in local ref; possible source of unverified eval numbers
