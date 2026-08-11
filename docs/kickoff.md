# kickoff — sciagent extension grilling

> **Outcome note (superseded).** The decisions recorded here (§9) fed
> [`docs/implementation-kickoff.md`](implementation-kickoff.md), which planned
> an `inject --tag`/`eject`/`metadata.requires:`/`tags.yaml` extension. That
> plan did not ship as designed — a later "demolish the role layer" refactor
> went the opposite direction and removed `inject`/`eject`, the `metadata:`
> frontmatter block, and all `requires:`/tag mechanisms (see
> `docs/architecture.md`). This file's grilling record and the §7
> session-resume prompt are kept as audit trail; do not treat them as a
> description of the current system or as work to resume.

A working document for resuming the sciagent extension design review. Serves two purposes simultaneously:

1. **A guide for the human** (Anton) — what to read, in what order, what to decide, when to step back, what to hand-edit.
2. **A session-restoring prompt** — paste this file's contents (or its body — see §7) into a fresh Claude Code session to pick up exactly where the previous one left off.

The same content serves both because the instructions a future reader needs are the same instructions a future assistant needs.

---

## 1. Where things stand

The toolkit at `01_modules/SciAgent-toolkit/` has a proposed 10-ADR extension spec at `docs/proposals/sciagent-extension-design-spec.md`. Three research passes have produced eleven supporting files at `docs/proposals/ai-research/01-...md` through `11-architectural-backbone-and-trace-necessity.md`. No code has changed. The grilling is in progress; decisions are pending.

**Verified facts** (from research, not vibes):
- ERA paper exists — Aygün et al., Nature 2026; code at `github.com/google-research/era` (Apache-2.0). Algorithm is "Flat UCB Tree Search," not "UCB1." (File 06.)
- agentskills.io spec exists; `skills-ref` CLI does single-skill structural validation only — cross-skill graph checks (cycles, resolution) are sciagent's job. (File 05.)
- Pi coding agent has a fully documented hook/event system and `version: 3` JSONL session schema, but **no Task-tool / sub-agent primitive**. (File 02 + `docs/.ref/pi/`.)
- Google Antigravity science-skills repo ships `workflow_skill_creator` — but **interview-driven**, not trace-driven (four-phase gate model). (File 04 + `docs/.ref/science-skills/`.)
- "Antigravity Skills paper by Anthropic" does **not** exist — the spec's references to it are misattribution. ADR-005 is genuinely net-new design. (File 04.)
- Four community Pi extensions (`pi-actors`, `pi-multiagent`, `pi-subagents`, `pi-packages`) converge on a tiny shape: one `package.json` declaring `"pi".extensions`, one or two `pi.registerTool` calls, optional skills bundle, optional command, optional `session_shutdown` hook. None is a clean dependency. (File 11 §1.)

**Standing preferences** (in auto-memory, will load automatically next session):
- **Harness-agnostic preference** — prefer harness-agnostic designs over coupled, modulo engineering cost.
- **Minimal but no less** — start at upstream ecosystem minimum; each addition must justify against a *current concrete* use case. Grill where the seam crosses.

---

## 2. The decision taxonomy

Three kinds of decision live in this spec. Don't conflate them — each has a different *right answer process*. This is the load-bearing frame for the whole grilling session.

### [Arch] — Architectural
Hard to reverse. Cross-cutting. Shapes the system for everything downstream.
- **How to approach**: think slowly. Sketch alternatives. Quantify the cost of being wrong. Sleep on it.
- **Decision rule**: pick the option whose worst-case regret you can live with.
- **Examples in this spec**: harness-agnostic-vs-coupled deployment shape (file 11 §2), trace-substrate yes/no (file 11 §3), namespace metadata under `metadata.sciagent.*` (ADR-001), depth-2 stack cap (architecture spec).

### [Eng] — Engineering
Cost-quality tradeoffs in implementation. Reversible at refactor cost. Local impact.
- **How to approach**: estimate cost honestly, estimate value honestly, decide. If wrong, eat the refactor.
- **Decision rule**: pick the cheapest option that meets the bar.
- **Examples in this spec**: build own validator vs. wrap `skills-ref` (ADR-007), `uv` inside container (ADR-009), how much LOC for per-harness adapter (file 11 §2 option B).

### [Taste] — Taste
Multiple defensible answers; no objective winner. Workflow-dependent.
- **How to approach**: pick fast. Document *why*. Don't relitigate later.
- **Decision rule**: any defensible answer is fine; the cost of *not deciding* is higher than the cost of deciding wrong.
- **Examples in this spec**: which two-tier names (ADR-003: `orchestrator`/`atomic` vs. others), whether `scorer` is sub-agent or tool (file 11 §4a), big-bang vs. per-skill migration for ADR-001.

**Common failure mode**: treating an [Arch] question as [Taste] and shipping the first defensible answer. Treating a [Taste] question as [Arch] and spending three weeks deliberating. Tag the question before answering it.

---

## 3. The decision matrix (ADR-by-ADR)

> **Status (2026-05-24)**: all rows resolved. The grilling phase is complete. The table below is retained as an audit trail — strikethrough rows trace what was asked, how it was tagged, and where the resolution landed (§9 for decisions; `proposals/deffered.md` for deferrals). When the next design question arises, add a row here with its [Arch]/[Eng]/[Taste] tag and re-enter grilling mode.

Each ADR's outstanding question, tagged. See referenced files for full grilling material.

| ADR | Question | Tag | Where to grill |
|-----|----------|-----|----------------|
| 001 | ~~Big-bang codemod or per-skill migration?~~ | — | **DECIDED 2026-05-24** → §9 (big-bang) |
| 001 | ~~Keep volume of custom metadata, or trim toward upstream minimum?~~ | — | **DECIDED 2026-05-24** → §9 |
| 001 | ~~Tag vocabulary location/format~~ | — | **DECIDED 2026-05-24** → §9 (`tags.yaml` at toolkit root) |
| 001 | ~~Tag vocabulary seed list~~ | — | **DECIDED 2026-05-24** → §9 (10-tag seed) |
| 001 | ~~Tag vocabulary mutability~~ | — | **DECIDED 2026-05-24** → §9 (PR-only) |
| 001 | ~~`inject --tag <name>` × stack-depth cap~~ | — | **DECIDED 2026-05-24** → §9 (separate concept, cap stays at 2; revisit if uncomfortable) |
| 001 | ~~`inject --tag` ejection~~ | — | **DECIDED 2026-05-24** → §9 (symmetric `eject` verb, generalizes over both named-skill and `--tag` forms) |
| 002 | ~~Hard `requires` contract~~ | — | **DECIDED 2026-05-24** → §9 (warn-and-continue, visibly) |
| 002 | ~~`requires` warning surface format~~ | — | **DECIDED 2026-05-24** → §9 (STDERR end-of-activate summary block — simplest) |
| 002 | ~~Semver constraints on `requires`~~ | — | **DECIDED 2026-05-24** → §9 (plain names; one version of a skill lives in the repo at a time) |
| 003 | ~~Two-tier scope names~~ | — | **DECIDED 2026-05-24** → §9 (`concept` / `implementation`) |
| 004 | ~~All sub-questions~~ | — | **DEFERRED 2026-05-24** → `proposals/deffered.md` |
| 005 | ~~All sub-questions~~ | — | **DEFERRED 2026-05-24** → `proposals/deffered.md` |
| 006 | ~~All sub-questions~~ | — | **DEFERRED 2026-05-24** (with ADR-004) → `proposals/deffered.md` |
| 007 | ~~Build own validator or wrap `skills-ref`?~~ | — | **DISSOLVED 2026-05-24** → §9 (no subsystem; `sciagent validate` is a thin callable verb that `activate` invokes internally) |
| 007 | ~~Validator tag extension~~ | — | **DISSOLVED 2026-05-24** → §9 (tag-vocab check is one of the ~30 LOC of cheap checks inside `validate`) |
| 008 | ~~`lab-loop` template scope~~ | — | **DROPPED 2026-05-24** → `proposals/deffered.md` (motivation died with ADR-004 deferral) |
| 009 | ~~`uv` mandatory, optional, or container-fast-path?~~ | — | **DEFERRED 2026-05-24** → `proposals/deffered.md` (Anton works in Docker with pre-installed envs; not pressing) |
| 010 | ~~Cloud-lab deferred — anything to change in the deferred shape?~~ | — | **No change** — leave ADR-010 in its existing deferred state |

**Cross-cutting decisions** — all three earlier entries (trace/result substrate, stack-depth cap under scorable+lab-loop, out-of-process default) are dissolved by the 2026-05-24 deferrals of ADR-004/005/006. The only surviving stack-depth-cap question lives under ADR-001's `inject --tag × depth-cap` row.

---

## 4. ADR-004's sub-agent bisection (worked example)

> **Note (2026-05-24)**: ADR-004 is now deferred wholesale (see §9 + `proposals/deffered.md`). The bisection below is preserved as an illustrative example of the [Arch]/[Eng]/[Taste] taxonomy in action — useful pattern for grilling other ADRs. The underlying sub-agent decisions are themselves moot until ADR-004 reactivates.

This is a worked example of how the taxonomy bites a single ADR. Use it as a pattern when grilling the others.

The spec proposes three sub-agents inside `scorable`: `mutator`, `scorer`, `selector`. The agents-best-practices rule "*subagents only when decomposition improves measured results*" (`docs/.ref/agents-best-practices/SKILL.md:172`) bisects them differently:

- **`selector` → [Eng]**. It is deterministic Python. Spec justifies it as sub-agent "for tool-use parity" — architectural prettiness. Drop the sub-agent framing; use a tool. Cheapest option meets the bar. *Decide now.*
- **`scorer` → [Taste]**. Needs sandboxed code execution — structural argument, not measured-decomposition argument. Sub-agent or tool both defensible; choice is about retry semantics and error surfaces. *Decide fast, document why.*
- **`mutator` → [Arch] (deferred on borrowed evidence)**. The only case where the decomposition rule genuinely applies. ERA's paper measures decomposed-vs-monolithic on their benchmarks; decomposed wins. Accept borrowed evidence for v1, mark a re-A/B trigger ("revisit when tier-1 lands"). *Decide deliberately; mark the re-eval.*

This is what "minimal but no less" looks like in practice: three sub-agents collapse to one + two tools without losing any measured benefit. Apply the same lens elsewhere.

---

## 5. Reading order (for a coffee, sitting down)

If reading cold or after a break:

1. **`docs/architecture.md`** — current sciagent system. Re-grounds you in what exists.
2. **`docs/proposals/sciagent-extension-design-spec.md`** — the spec being grilled. Skim ADR titles + executive summary. Don't try to absorb all 10 ADRs at once.
3. **This file** — the decision taxonomy + matrix (§§2–3 above).
4. **`docs/proposals/ai-research/11-architectural-backbone-and-trace-necessity.md`** — the single most consequential research output. Two big calls (deployment shape, trace necessity) + the sub-agent bisection. *Most actionable file in the set.*
5. **`docs/proposals/ai-research/10-ADR-on-ADRs.md`** — the deep grilling matrix. Use as reference, not cover-to-cover.
6. **`docs/proposals/ai-research/04-...md`** — the misattribution finding. Recalibrates how to read ADR-005.
7. Drill into other files (01, 02, 03, 05–09) as specific ADRs come up in the grilling.

---

## 6. Step-back moments

Before each [Arch] decision, write down the answers to these out loud (notepad, comment, whatever):

1. **What's the cost of being wrong?** Hours of refactor? Public commitment? Loss of optionality?
2. **Is the seam load-bearing for your specific workflow?** Not "could this matter for a hypothetical user." For *your* concrete next 90 days of work, what breaks if you choose differently?
3. **Does this make harness-coupling more or less likely?** If more, what concrete value is it buying? (Per harness-agnostic preference.)
4. **What's the upstream minimum?** Could this decision be replaced with "do what agentskills.io / Pi / Claude Code already does"? If yes, pick that unless you have a current concrete reason not to.
5. **Is this actually one decision, or several stapled together?** Spec ADRs are sometimes bundles (ADR-005 was — split into 005a + 005b). Split before deciding.

For [Eng] decisions, skip 1 and run 2–5 fast. For [Taste], pick and move on.

---

## 7. Session-resume prompt

Paste the body below into a fresh Claude Code session to restore context. It assumes the model has access to the working directory `/workspaces/DC_hum_verse/01_modules/SciAgent-toolkit/`.

> ---
>
> Resuming the sciagent extension grilling session. Three preferences already in memory and will load automatically: harness-agnostic-preference, minimal-but-no-less, sciagent-extension-grilling project context. No need to re-establish those.
>
> **State of play**: the 10-ADR spec at `docs/proposals/sciagent-extension-design-spec.md` has been through three research passes. Outputs are in `docs/proposals/ai-research/01-...md` through `11-architectural-backbone-and-trace-necessity.md`. The decision taxonomy (architectural / engineering / taste) and per-ADR matrix live in `docs/kickoff.md` (this file).
>
> **Verified facts to anchor on**: ERA paper is real (Nature 2026, github.com/google-research/era). The "Antigravity Skills paper by Anthropic" does NOT exist (misattribution; the actual Google Antigravity science-skills repo at `docs/.ref/science-skills/` ships an interview-driven `workflow_skill_creator`, not trace-driven). Pi coding agent has a documented hook/event surface but no Task-tool / sub-agent primitive. agentskills.io's `skills-ref` validator does single-skill checks only — cross-skill graph is sciagent's job. Four community Pi extensions converge on a small backbone (see file 11 §1).
>
> **Two prior recommendations on the table**:
> 1. Ship sciagent as a harness-agnostic core + thin per-harness adapters (~150 LOC each). Do not depend on any single Pi extension. (File 11 §2.)
> 2. Defer JSONL trace recording (ADR-005b) until ADR-004 or ADR-006 also commits this cycle; ship ADR-005a (interview-driven distiller) first. (File 11 §3.)
>
> **ADR-004 has been bisected**: `selector` → tool (drop sub-agent framing, [Eng]); `scorer` → [Taste]; `mutator` → keep as sub-agent on borrowed ERA evidence, re-A/B when ADR-006 tier-1 lands ([Arch] deferred). (File 11 §4a.)
>
> **Today's task**: I am sitting down to make decisions on the queue in `docs/kickoff.md` §3. I will work through them tagged by category. For each:
> - [Arch] questions: I want you to push back, surface alternatives I might not have considered, and quantify cost-of-being-wrong honestly. Run §6's step-back checklist with me.
> - [Eng] questions: estimate the cheapest option that meets the bar; propose; defer to me.
> - [Taste] questions: name two or three defensible answers, pick the one I lean toward, document the *why* in one sentence.
>
> **Hard constraints**: no code changes; no new ADRs unless I explicitly ask; no edits to the original spec at `sciagent-extension-design-spec.md`. Updates go to files in `docs/proposals/ai-research/` or `docs/kickoff.md`. Voice in artifacts is impersonal; chat stays conversational.
>
> Confirm context is loaded and ask me which decision I want to start with — or wait for me to point at one.
>
> ---

---

## 8. After the next session

Once decisions land, the kickoff doc itself is the thing to update. The pattern:

- Move decided items out of §3 into a new §9 "Decided" with the resolution and a one-sentence why.
- Re-tag any question whose category turned out wrong (a question pre-tagged [Taste] that became [Arch] during grilling — record that lesson).
- If new questions surface, add them to §3 with a tag. Don't lose them.
- When §3 is empty, ADRs are ready for implementation planning — which is a different document, not this one.

This file should shrink over time, not grow. If it's growing, the grilling is producing new questions faster than it's closing old ones — that's a signal to step back to §6's checklist again.

> **2026-05-24**: §3 emptied. Implementation planning lives at [`docs/implementation-kickoff.md`](implementation-kickoff.md) — five-PR breakdown plus the same dual-use guide/session-resume prompt shape this file used.

---

## 9. Decided

Items that have moved from §3 to here. Each entry: the question, the resolution, and a one-sentence *why* — enough that future-you (or a future assistant) can challenge the call without reconstructing the conversation.

### 2026-05-24 — ADR-001 metadata volume → **A. Keep all four custom fields + add `tags`.**

Resolution: `metadata.sciagent.scope`, `requires`, `complementary-skills`, `contraindications`, and `tags` all live in skill frontmatter. Skills do not plug into sciagent's role/stack system by default — adding the frontmatter is the cost of admission.

**Why**: opinionated metadata *is* the architecture. The act of authoring frontmatter is the seam where a generic skill becomes a sciagent-aware skill. A generic LLM agent can do this tinkering for any human. The deviation from agentskills.io's minimum is intentional and load-bearing, not bloat.

[Arch] tag held; chose the more expressive option after grilling against "minimal but no less." The seam was named (every skill needs metadata to play along with the system) and that justified the addition.

### 2026-05-24 — Inject-by-tag mechanism → **A. Extend ADR-001 + the existing `inject` verb.**

Resolution: a new `tags` field on skill metadata (covered by the ADR-001 decision above), plus an extension of the existing `sciagent inject` verb with `--tag <name>` semantics. No new ADR-012; the work fits inside ADR-001's scope plus a small CLI extension.

**Why**: keeps the grilling surface contained. Tags are already metadata once ADR-001 lands; injection-by-tag is a natural inject-verb extension, not a new architectural concept. RPG framing: an active stack of two roles can receive a "buff" of skills sharing a tag — useful when no third role exists for the precise combination needed.

Follow-up questions surfaced and added to §3: vocabulary location, vocabulary seed, vocabulary mutability, stack-depth-cap interaction (the [Arch] one — is the buff a third stack frame or a separate concept?), and validator coverage of tag values.

### 2026-05-24 — ADR-005 → **DEFERRED wholesale.**

Resolution: full ADR-005 (both proposed 005a interview-driven distiller and 005b JSONL recorder) deferred. Replacement in scope now: a docs/codebase → skill(s) generator command, recorded as a future ADR placeholder. Full deferral context in `docs/proposals/deffered.md`.

**Why**: matches the user's current personal workflow; no second consumer for traces; the recorder is harness-coupled and violates the harness-agnostic preference without offsetting benefit. Revisit triggers documented in `deffered.md`.

### 2026-05-24 — ADR-001 tag vocabulary location → **A. Single `tags.yaml` at toolkit root.**

Resolution: tag vocabulary lives in a single `tags.yaml` at the toolkit root with `name`, `description`, `since` fields per entry. Diffable, greppable, single source of truth.

**Why**: minimal-but-no-less applied to the vocabulary itself — one file beats union-of-role-yamls (no SSoT) and beats per-skill ad-hoc tag invention (no governance). The cost of one extra file is trivial; the benefit (controllable, reviewable, evolvable as one artefact) is concrete.

### 2026-05-24 — `inject --tag <name>` mechanism scope → **B. Separate concept from "stack"; cap stays at 2.**

Resolution: `inject --tag` writes flat into `.claude/` and `.agents/` without participating in the role-stack model. The architectural cap of 2 (base + overlay) stays. Tag-injection is a workaround / buff mechanism, not a third stack frame.

**Why** (Anton verbatim): *"Right now, this independent, non-overlapping concept of injecting by the tag is a workaround and let's try to keep it like that for now, but we'll see."* Conservative call — preserves the depth-2 cap as load-bearing; if the workaround turns out to be insufficient or confusing, revisit then, not now.

Follow-up surfaced and added to §3: how to *eject* a tag-buff without disturbing the stack ([Eng]).

### 2026-05-24 — ADR-002 `requires` contract hardness → **B. Warn-and-continue, visibly.**

Resolution: activation does not block on missing `requires`. The user is loudly warned and prompted to investigate (manually or via an agent). The warning must be clearly visible and readable so the user can either fix the gap or judge it non-breaking and proceed.

**Why** (Anton): don't break the user's flow. If a missing dep is slop-level non-breaking for the current workflow, the user can continue. But the warning has to surface — silent failure is the actual worst case.

Follow-up surfaced and added to §3: warning surface format ([Eng]) — where it appears, whether machine-parseable so an agent can act on it.

### 2026-05-24 — ADR-003 two-tier scope names → **`concept` / `implementation`.**

Resolution: `metadata.sciagent.scope: concept` for skills that capture durable orientation / theoretical understanding / "step-back" perspective; `metadata.sciagent.scope: implementation` for skills that bundle practical snippets, scripts, references, specific swappable tools.

**Why**: captures the **durability axis** Anton articulated — concept-skills persist as taste/understanding grows across tooling changes; implementation-skills are expected to obsolete and get swapped. The word *implementation* naturally implies interchangeability (you can swap implementations of the same concept), which is a stronger semantic for swappability than alternatives like `tool` or `leaf` carried.

Rejected pairs (recorded so the call is not relitigated): `orchestrator`/`atomic` (opaque), `flow`/`leaf` (visual but semantically off), `composite`/`leaf` (technical but cold), `principle`/`practice` (practice undersells swappability), `concept`/`tool` (close — but lost to `implementation`'s framing).

### 2026-05-24 — ADR-001 migration approach → **A. Big-bang codemod.**

Resolution: one PR adds the new metadata fields (`scope`, `tags`, `requires`, `complementary-skills`, `contraindications`) to all existing skills at once via a mechanically-generated codemod. Empty values where not yet hand-curated; `scope` defaulted by best-guess heuristic and reviewed in the same PR.

**Why**: tight diff, mechanically generated, complete in one pass. Avoids the long-lived inconsistency between migrated and unmigrated skills that per-skill migration produces; avoids the "this skill works differently because it predates the convention" failure mode.

### 2026-05-24 — ADR-001 tag vocabulary mutability → **A. PR-only.**

Resolution: adding a tag to `tags.yaml` requires editing the file and shipping a PR / commit. No `sciagent tags add <name>` CLI verb. No auto-add on first use.

**Why**: the vocabulary is the small file that governs every skill's metadata. The cost of a slightly-friction-y add (open file, edit, commit) is negligible compared to the cost of a degraded vocabulary that accumulates typos, synonyms, and one-off tags via convenience verbs. Friction here is a feature.

### 2026-05-24 — ADR-001 tag vocabulary seed list → **DEFERRED to a separate quick pass.**

Resolution: only `trajectory` (or `pseudotime`) named so far; full day-one seed list to be sketched in a focused short session against Anton's actual project surface. Not blocking other decisions — the vocabulary file ships empty-or-near-empty at first.

Tracked in §3 as a pending sub-decision (not a `deffered.md` candidate — this is a sketching task, not a design-deferral).

### 2026-05-24 — ADR-001 `inject` ↔ `eject` verb pair → **Symmetric, generalized over both forms.**

Resolution: a new `sciagent eject` verb mirrors `sciagent inject` end-to-end:
- `sciagent inject <skill-name>` ↔ `sciagent eject <skill-name>` — named single skill.
- `sciagent inject --tag <name>` ↔ `sciagent eject --tag <name>` — tag-bundled group.

`eject` operates only on previously-injected items. Skills mounted as part of a role-stack layer are removed via `deactivate`, not `eject` — the two verbs cover orthogonal concepts (stack composition vs. ad-hoc injection).

**Why**: symmetry kept the user from asking "where did my X go?" The original Q12 only proposed eject for the `--tag` form; Anton flagged that the asymmetry would be confusing — eject should cover the named-skill form too, by the same rule. Clean generalization.

### 2026-05-24 — ADR-002 `requires` warning surface → **A. End-of-activate STDERR summary block.**

Resolution: missing `requires` produce a formatted, human-readable summary printed at the end of `activate` to STDERR. One place to look. Exit code remains 0 (warn-and-continue, per the earlier hardness decision). No JSON sidecar, no per-line inline noise.

**Why** (Anton): simplest. The other options (per-line + summary; STDERR + JSON sidecar) hypothesise an agent-consumer that doesn't exist today. Minimal-but-no-less.

### 2026-05-24 — ADR-002 semver constraints on `requires` → **A. Plain names.**

Resolution: `requires: [skill-a, skill-b]`. No `@^1.0`, no version-range syntax, no version-pinning grammar.

**Why** (Anton): only one version of each skill lives in the repo at a time. Updating a skill means editing the skill and bumping its `version:` for diff-tracking — old versions are not retained. A `requires` constraint over a single-version-at-a-time skill repo would be either trivially satisfied or impossible; the syntax buys nothing.

### 2026-05-24 — ADR-007 validator → **DISSOLVED. `sciagent validate` is a callable verb, not a subsystem.**

Resolution: a thin `sciagent validate` verb runs ~30 LOC of cheap checks over the skills/roles graph:
- `requires` resolution (each named dep exists as a skill).
- Cycle detection (DFS over `requires`).
- Tag-vocab compliance (each declared tag is in `tags.yaml`).
- *Optional*: shell-out to `skills-ref skill_name` per skill for per-skill structural validation, if `skills-ref` is installed; silently skipped otherwise.

`sciagent activate` calls `validate` internally before mounting. The standalone verb stays useful for debugging — run `sciagent validate` to see graph health without re-mounting.

**Why** (Anton verbatim): *"I think validate may be a separate word, but whenever you activate, activate just calls for validate internally so that we maintain validate separate from activation and it would be possibly available for like debugging purposes as a dedicated verb. … It's like three uh thirty lines of code, why not?"*

The earlier framing — "build own vs. wrap `skills-ref`" (Q15) and "validator tag extension" (Q16) — both dissolved here. Neither was a real architectural question; both presumed a validator subsystem that turns out not to exist. Recorded as a lesson in over-framing: tag every question with the *taxonomy* before grilling, not just the ADR number.

### 2026-05-24 — ADR-008 `lab-loop` → **DROPPED, folded into `proposals/deffered.md`.**

Resolution: ADR-008 dropped from the active spec. Its proposed `lab-loop` overlay was justified as a swap-target stub for the ADR-004 optimizer; with ADR-004 deferred, the justification is gone. Not killed forever — reactivation conditional on ADR-004 reactivation. Brief note appended to `deffered.md`.

### 2026-05-24 — ADR-009 `uv` → **DEFERRED, folded into `proposals/deffered.md`.**

Resolution: ADR-009 (`uv` pinning strategy for skill Python deps) deferred. Anton's current workflow is Docker-container-with-pre-installed-envs; venv management is not a current pain point.

**Why** (Anton): *"UV seems to be optional at the moment. I personally work mostly within Docker containers. I have everything pre-installed there and I'm not that concerned with managing virtual environments."*

Revisit trigger: sciagent gets adopted by a non-containerized workflow, or a collaborator needs reproducible env management without Docker.

### 2026-05-24 — ADR-010 cloud-lab deferred shape → **No change.**

Resolution: ADR-010 was already deferred in the original spec. The grilling question ("anything to change in the deferred shape?") resolves to: no. Leave it as-is. Removed from active grilling queue.

### 2026-05-24 — PR 2 hardness boundary (post-grilling clarification surfaced during implementation-kickoff review) → **(a) Keep resolver hard-fail on `requires`; soft-warn only for `complementary-skills`.**

Resolution: missing `requires` continue to hard-fail in `skill_deps.sh` (preserving the existing pre-mutation safety invariant in `activate.sh` Phase B). Soft warn-and-continue applies *only* to `complementary-skills` references that don't resolve. The STDERR end-of-activate summary block carries the complementary-skill warnings; `requires` failures still abort before any symlink lands.

**Why**: honors §9's spirit (no silent failures, user clearly sees what's off) without removing the existing pre-mutation safety invariant. `requires` is a *contract* — the skill literally cannot function without those deps; half-mounting would produce silent breakage. `complementary-skills` is a softer signal where warn is the right level.

Reading (b) — softening the resolver to literal warn-and-continue on `requires` — was rejected as it would remove the safety guarantee against half-mounted broken skills.

This refines (does not contradict) the earlier ADR-002 hardness decision in this same §9. The earlier decision said "warn-and-continue, visibly"; that phrasing turned out to elide the `requires` vs `complementary-skills` distinction. Reading (a) preserves both halves: hard contract enforcement on `requires`, soft visible warn on `complementary-skills`.

### 2026-05-24 — PR 1 codemod scope mapping → **Confirmed: `foundation` + `orchestrator` → `concept`; `atomic` → `implementation`.**

Resolution: the heuristic proposed in `implementation-kickoff.md` PR 1 is accepted as the codemod default. Hand-flips during code review for any individual skill where the mapping looks wrong on inspection.

**Why**: `foundation`/`orchestrator` skills (e.g. `scvi-framework`, `architecture-first-dev`) are durable concept-holders aligned with the `concept` semantic — they describe systems-level orientation that persists across tooling changes. `atomic` skills (e.g. `scvi-basic`, `pycistopic-atac-topic-modeling`) are practical tool wrappers aligned with `implementation` — they ship concrete swappable kit.

### 2026-05-24 — PR 1 body-line caps + folder relief + lint side items → **`concept` ≤ 500, `implementation` ≤ 350, folder relief blessed for both scopes; default-on-missing = `implementation`; lint CUTOFF bumped to 2026-05-24; new regression-block lint rejecting old scope names.**

Resolution: under the new `concept`/`implementation` split, `tests/test_skill_scope_lint.sh` enforces:
- `concept` skills ≤ **500** body lines (after frontmatter, fenced code excluded — the existing `body_loc_no_fences` counting rule is unchanged).
- `implementation` skills ≤ **350** body lines (same counting).
- **Folder-relief pattern blessed for both scopes**: a skill folder MAY contain reference files (scripts, snippets, data, templates) beyond `SKILL.md`, referenced by basename from the SKILL.md body. Not a new convention — codifies the existing one (12/62 skills already use it).
- Default `scope` when frontmatter field is missing: **`implementation`** (was: `atomic`). Tracks the migration mapping.
- Legacy WARN-vs-FAIL CUTOFF bumped 2026-05-21 → **2026-05-24** (the migration decision date) so newly-relabeled skills are held to the new policy strictly while pre-migration WARNs remain warnings.
- **New regression-block lint rule**: presence of `scope: atomic|orchestrator|foundation` in any non-`_TEMPLATE` skill is a FAIL — prevents accidental re-introduction of the old vocabulary.
- `architecture-treemap` (currently missing the `scope:` field) gets `implementation` assigned during the codemod, per Anton's confirmation.

**Why**: empirical scan (Opus agent, 2026-05-24) of all 62 skills:
- Concept population max is 381 body lines (`skill-creator`); 500 gives ~30% margin while tightening from the legacy 800 ceiling.
- Atomic (→ implementation) population max is 322 (`scrna-cxg-host`); 4 skills would breach a 300 cap (322/317/317/310). 350 catches all four with ~10% headroom.
- One outlier — `iterative-peak-merging` (317 body lines, **zero companion files**) — could not be relieved by folder-relief alone and would require a real refactor at 300. The 350 cap absorbs it without forcing a code change unrelated to the migration.
- Folder relief is already used by 12/62 skills (9 atomic, 2 foundation, 1 missing-scope). Spot-checks confirm SKILL.md bodies reliably reference companions by exact basename. Restricting the pattern to `implementation` only would force a refactor of `skill-creator` (381-line `concept` post-migration that uses the pattern heavily).

The numbers preserve the "minimal but no less" pressure (tighter than legacy 800) while not forcing refactors of skills that work today.

### 2026-05-24 — ADR-001 tag vocabulary seed → **10 tags, ship as `tags.yaml` content below.**

Resolution: day-one `tags.yaml` ships with the ten tags below. Comp-bio core (8) + `architecture` + `tooling` to seed the non-comp-bio direction. Naming follows the **shortest-unambiguous-form** rule (case-by-case): `de`, `qc`, `viz` shorten where the full form is noisy; `trajectory`, `integration`, `annotation`, `pathway`, `report`, `architecture`, `tooling` stay long because they already read fast.

Notable choices:
- `gsea` renamed to `pathway` — broader, covers enrichment + activity scoring + signatures + classic GSEA. Tag is the bundling concept; sub-methods live in skill descriptions.
- `docs` merged into `report` — same conceptual surface for the user; one tag avoids ambiguity at activation time.
- `debug`/`testing` not seeded; add via PR when a concrete skill needs the bucket. Per the mutability decision (PR-only), this is cheap.

**Proposed `tags.yaml` content** (records the decision; actual file gets created in the implementation phase, not now):

```yaml
# tags.yaml — sciagent skill tag vocabulary.
# Single source of truth. Adding a tag = edit + PR.
# See docs/proposals/sciagent-extension-design-spec.md ADR-001.

tags:
  - name: trajectory
    description: Time-trajectory inference (pseudotime, RNA velocity, dynamics).
    since: 2026-05-24

  - name: integration
    description: Batch correction and multi-sample harmonization (scVI, scANVI, harmony, MNN).
    since: 2026-05-24

  - name: annotation
    description: Cell-type assignment (marker-based, reference-projection, manual rescue).
    since: 2026-05-24

  - name: de
    description: Differential expression (limma, edgeR, DESeq2, scanpy.rank_genes_groups).
    since: 2026-05-24

  - name: pathway
    description: Pathway analysis — enrichment (fgsea, hypergeometric), activity scoring (GSVA, AUCell, decoupleR), signatures.
    since: 2026-05-24

  - name: qc
    description: Quality control (per-cell, per-sample, batch-level).
    since: 2026-05-24

  - name: viz
    description: Plotting, dimensionality-reduction display, figure composition.
    since: 2026-05-24

  - name: report
    description: Knit/render, publication-figure pipelines, captions, README/Methods generation, skill documentation.
    since: 2026-05-24

  - name: architecture
    description: Software architecture — design patterns, refactoring guidance, system structure.
    since: 2026-05-24

  - name: tooling
    description: General-purpose tool/utility building (CLIs, MCP servers, file-handling helpers, dev workflows).
    since: 2026-05-24
```

### 2026-05-24 — ADR-004 (optimizer overlay) + ADR-006 (benchmark harness) → **DEFERRED wholesale, as a pair.**

Resolution: both ADRs deferred together; sciagent ships as a fully opinionated personal context manager with no objective fitness signal. Full deferral context, revisit triggers, and reactivation surface in `docs/proposals/deffered.md`.

**Why**: no concrete second consumer for an optimizer; no clear benchmark tasks exist yet; N=1 self-evaluation is the wrong protocol; zero lock-in cost from waiting (FUTS is 154 LOC Apache-2.0, borrowable anytime). Anton: *"I'm not really sure how to benchmark which tasks to benchmark on. So for now we keep this opinionated."*

ERA-inspection pass (Opus agent, 2026-05-24) concurs: solo-researcher use case plus Anton's own "minimal but no less" rule jointly conclude defer. Path (c) of three options — no benchmark, fully opinionated.

The sub-agent bisection captured in §4 is preserved as illustrative of the taxonomy in action; the underlying decisions are moot under this deferral.

---

## 10. Decisions pending external input

Currently empty. The ADR-004 staged decision resolved on 2026-05-24 once the ERA-inspection pass landed (see §9). Populate this section when the next decision in §3 needs an external input before it can land.

---

*End of kickoff. Decisions live in §9; pending external inputs in §10. Walk in, sit down, pick one from §3.*
