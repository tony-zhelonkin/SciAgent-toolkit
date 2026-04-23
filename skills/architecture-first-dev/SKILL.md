---
name: architecture-first-dev
description: Architecture-first software development methodology. Use when planning non-trivial features, refactors, or redesigns where upfront design pays off. Encodes the 6-stage per-feature pipeline (map → review → synthesize → design → plan → implement) plus the optional 3-command meta-layer (meta-map → meta-design → meta-plan) for cross-feature work. Primes the assistant to stop at each gate, surface per-stage tips, and navigate the user from any state. Defer to docs/workflows/architect/00-quickstart.md as the cadence source of truth.
---

# Architecture-First Development

A discipline for non-trivial software changes: design before code, with composable expert review.

## For the assistant — navigation contract

**You are responsible for helping the user navigate this harness.** The user does not need to memorise stages, gates, or filenames; they ask, you answer. Your sources of truth, in priority order:

1. **`docs/workflows/architect/00-quickstart.md`** — cadence, decision tree, per-stage tips, prototypical meta-design walkthrough. Read this for any "what should I run next?" / "is meta-layer worth it?" / "why did X stop?" question.
2. **`docs/workflows/architect/02-behavior.md`** — exact per-stage behaviour, including the meta-stages (3a `/meta-map`, 3b `/meta-design`, 4b `/meta-plan`).
3. **`docs/workflows/architect/05-meta-approach.md`** — meta-layer reference; when, why, worked example.
4. **The command file** (`commands/<name>.md`) for the exact phase logic — when a command's behaviour is in question.

When the user lands mid-flow ("I have map.md and one review, now what?") inspect `docs/{slug}/` and `docs/_meta/` and recommend the next single command, with a one-line "why this and not the alternatives." Always offer the alternatives explicitly so the user can redirect.

**If the user feels disoriented ("where are we?", "what's the state?", "what's left?"), recommend `/status` first.** It produces a one-page portfolio-state summary in <30 seconds — much cheaper than answering those questions by reading individual artifacts. After `/status`, recommend the specific next command based on the surfaced state.

After every command finishes, surface the **Tips after this stage** entries from `00-quickstart.md` for that stage. Don't make the user dig for them.

When recommending a command, name the alternative paths (e.g., "/review --as all is overkill if you only need a stat lens — `/review --as stat` is cheaper"), and justify the recommendation in one sentence.

---

## When to use this methodology

**Use for:**
- New features or modules where the structure isn't obvious
- Refactors that touch ≥5 files or cross module boundaries
- Changes where scientific/biological reasoning matters as much as code structure
- Work that benefits from multiple expert lenses (statistical, ML, UX, wet-lab tractability)

**Skip for:**
- One-line fixes, typo corrections, trivial tweaks
- Throwaway scripts or exploratory data analysis
- Cases where design is already locked (just implementing an approved spec)

## The pipeline

```
/map <feature>                   — produces docs/{feature}/map.md
/review <feature> --as <spec>    — produces docs/{feature}/review/*.md (parallel)
/synthesize <feature>            — produces docs/{feature}/synthesis.md (if ≥2 reviewers)
/design <feature>                — produces docs/{feature}/design/*.md (gated by architect)
/plan <feature>                  — produces docs/{feature}/plan/phase-NN.md
/implement <feature> <phase>     — code, one phase at a time
```

### Stage 1: `/map` — codebase cartography

Runs the `mapper` agent (read-only). Produces a single artifact: `docs/{feature}/map.md` with blast radius (files, tests, configs touched), entry points, data flow with `file:line` references, existing patterns to reuse, and open questions. **This is the artifact all downstream stages read.** Never skip unless map.md already exists and is current.

### Stage 2: `/review` — composable expert panel

Runs 0-to-N reviewer agents in parallel:

| Flag | Runs |
|------|------|
| `--as all` | All 6 reviewers |
| `--as all --but stat,divergent` | All minus listed (exclusion syntax) |
| `--as bioinf,ml` | Just those |
| (no flag) | Asks user which to run |

Reviewers: `bioinf`, `wetlab`, `graphic`, `stat`, `divergent`, `ml`. Each reads `map.md` + feature scope and writes exactly one file to `docs/{feature}/review/<agent>.md`. **Reviewers never re-explore code** — if they need more, they write an "Open questions for mapper" section.

Skip this stage entirely for small features where a single perspective suffices. For medium features, invoke 1-2 reviewers. For large/cross-cutting features, invoke `--as all`.

**Two regimes, distinct epistemic purposes.** Round-1 (`/review --as ...`) is parallel-independent — reviewers are blind to each other; divergence is the signal the user is buying. Iterate (`/review --iterate --as ...`) is cross-informed — each reviewer sees every other reviewer's prior-round output plus the synthesis (incl. `§ User resolutions` where the user's locked direction lives); reconciliation is the point. Iterate reviewers are required to tag each prior position as UPDATED / STOOD BY / RETIRED so the position-shift is auditable. Prior rounds archive to `review/.history/round-N/` (never deleted). See MADR-H5 / ADR-019.

### Stage 3: `/synthesize` — consensus from reviews

Runs the `synth` agent. Reads all `review/*.md` files and writes `synthesis.md` identifying convergent concerns, preserved disagreements, and prioritized P0/P1/P2 recommendations. Conditional: only runs if ≥2 reviewers produced output (one reviewer's file is cheaper to read directly).

### Stage 4: `/design` — architecture views + gate

Main agent writes three design documents:
- `01-architecture.md` — structural view (components, boundaries); REQUIRED `## Delta view` subsection with an annotated Mermaid diagram (added nodes/edges green, modified yellow, removed red-dashed) so reviewers see *what changes* in one figure instead of mental-diffing prose against `map.md`.
- `02-behavior.md` — behavioral view (data flows, sequences)
- `03-decisions.md` — decision view (ADRs, trade-offs)

Gate-1 README now carries two coarse Mermaid diagrams (module sketch + data-flow sketch, ~5–10 lines each) — scope-level, non-binding; detailed diagrams still land in Phase 3b. This lets the user approve architectural scope with a picture, not prose alone.

Then the `architect` agent performs a consistency/completeness review → `design/review.md` with a Verdict field. Loop until READY.

### Stage 5: `/plan` — phased implementation plan

Decomposes into ordered phases (≤5 files per phase, testable each). Writes `plan/README.md` + `plan/phase-NN.md` files.

### Stage 6: `/implement` — one phase at a time

Executes a single phase, verifies against its checklist, reports completion, stops. User confirms before next phase.

## Meta-layer (optional, for ≥2 features in flight)

When two or more features are in flight at the same time and they touch overlapping surfaces, each feature's `/design` will make locally-optimal choices that can quietly contradict its siblings. The meta-layer is the optional tier above per-feature designs:

```
/meta-map       — produces docs/_meta/map.md (cross-feature touch-points + convergent concerns)
/meta-design    — produces docs/_meta/design.md (MADRs — inter-feature ADRs that per-feature designs inherit)
/meta-apply     — propagates accepted MADRs to every in-flight design/ in one parallel batch
                  (dispatches N feature-reviser sub-agents + N architect gates; one compound edit gate)
/meta-plan      — produces docs/_meta/plan.md (sequences per-feature phases across the portfolio)
```

When a per-feature `/design` runs and `docs/_meta/design.md` exists, it inherits every MADR. A local decision that conflicts with a MADR triggers a stop-and-ask gate — never silently overridden.

After `/meta-design` Gate 1 approval, `/meta-apply` is the primary propagation path: it replaces N serial `/design` re-runs with one parallel batch. Features whose MADR consequences need re-architecting (not just editing) are auto-flagged `MANUAL_REDESIGN_NEEDED` by their sub-agent and cleanly skipped — the user runs `/design <slug>` for those. See the "Two-path propagation" section below for the recommendation logic.

Skip the meta-layer when only one feature is in flight, or when in-flight features share zero touch-points. See [docs/workflows/architect/00-quickstart.md](../../docs/workflows/architect/00-quickstart.md) for the cadence at a glance and [05-meta-approach.md](../../docs/workflows/architect/05-meta-approach.md) for the full reference (when to invoke, worked example).

## Token-efficiency contract

**The single most important rule:** downstream stages MUST read upstream artifacts and MUST NOT re-grep the codebase. If map.md exists, reviewers use it. If design/ exists, architect reads it. This is how the pipeline stays affordable.

If upstream artifacts are stale or missing context, surface that in an "Open questions" section and have the user re-run `/map` — do not silently expand your own investigation.

## Human gates are non-negotiable — where judgment is required

Every stage that demands human judgment ends with a hand-off to the user. The assistant **waits** for approval before advancing. This is the feature, not the bug: the human is the interface on which architectural decisions happen.

**Gates exist where human judgment is required, not where execution is merely serialised.** See `docs/workflows/architect/03-decisions.md § ADR-018`. Two consequences to keep straight:

- **Judgment gates (preserved, always WAIT):** `/map` review, `/review` lens selection, `/synthesize` output, `/design` Phase 3a + Phase 4, architect `NEEDS_ITERATION`/`NEEDS_DISCUSSION`, `/verify` verdicts, `/meta-design` Gate 1, `/meta-apply` compound gate, `/meta-plan` sequencing, `/implement` stop-at-uncertainty.
- **Serialisation moments (relaxed, no WAIT):** `/plan` (plan README is stamped `status: APPROVED` at draft time; Phase 4 is a glance-check, silence-means-proceed), per-phase cadence in `/implement --auto` (end-to-end execution halts only on failure / uncertainty / completion).

The meta-architect's gates are MADRs, ADR content at `/design` Phase 4, `/meta-plan` cross-feature sequencing, and stop-at-uncertainty. They are NOT plan decomposition or per-phase cadence — those are DevOps-utility serialisations derivative of already-committed decisions.

**Scribe-on-latest — where binding chat decisions live.** Between commands the user routinely makes binding architectural decisions in chat ("actually, retire the tier badge entirely"). The harness does NOT have a separate chat-decisions file; instead, the assistant appends the decision to a `## User resolutions` / `## Binding chat decisions` / `## Post-READY amendments` section in the **most recent upstream artifact the next command will consume**. Every downstream command reads these sections and treats them as equal-priority to ADRs.

Selection table (keyed on file existence, deterministic — not a judgment call):

| State at time of chat | Append to | Section heading |
|---|---|---|
| `map.md` exists; no reviews | `map.md` | `## User resolutions` |
| Reviews exist; no synthesis | latest `review/*.md` (or new `review/user.md` if none is a natural home) | `## User resolutions` |
| `synthesis.md` exists | `synthesis.md` | `## User resolutions` |
| `design/03-decisions.md` exists (DRAFT or iteration loop) | `design/03-decisions.md` | `## Binding chat decisions` |
| `design/` READY, pre-`/plan` | `design/03-decisions.md` | `## Post-READY amendments` + re-dispatch architect |

**Cross-stage invalidation rule.** If a chat decision during a later stage logically invalidates an earlier frozen artifact (e.g., a mid-`/design` decision that overturns reviewer consensus), STOP and surface the conflict — same pattern as the MADR-conflict stop-gate in `/design` Phase 1:

```
This chat decision contradicts {frozen-artifact § section}.
Choose:
  1. Re-run /{upstream command} to revise the frozen artifact
  2. Override locally — I'll note the override in `## User resolutions` + cite the contradiction
  3. Cancel the chat decision
```

Never silently override. Never invent a new canonical sibling file (`chat-decisions.md`, `proposal.md`, etc.); the canonical file set is closed.

**Gate-framing rule: propose defaults, ask for exceptions.** When a gate surfaces Open Questions, present every OQ with a proposed default + one-line rationale, classified as `Domain? YES` (needs subject-matter call) or `Domain? NO` (architecture default the assistant is confident in). The user reads the defaults and flags only exceptions — cognitive work moves from *generation* ("what should k be?") to *discrimination* ("does single-median look right? yes/no"). Forcing a user to generate answers for architecture-level OQs cold produces under-informed decisions that are worse than a surfaced default. See `commands/design.md § Phase 4` for the pattern.

**Self-audit before presenting a gate**: if >50% of OQs carry `Domain? YES`, the design is asking too much of the user — pull the assistant's own position into ADRs (with TODO markers for validate-at-implement-time) and re-present. Most OQs should close under "approve all."

## Single-terminal portfolio workflow

When working across ≥2 features in a wave, **run the entire wave in one Claude Code session.** Features are topics inside a single long-running conversation, not separate terminal windows with independent Claude instances.

Why: every slash command opens a fresh context window, but within a single session the prompt cache stays warm, recently-read files are already in context, and cross-feature awareness (which MADRs apply, which features are APPROVED, which have open conflicts) lives in one place. Multi-terminal = N independent contexts = N times the re-orientation cost + zero in-session cross-feature awareness. The meta-layer was designed to operate on a coherent portfolio state; fragmenting the session fragments exactly what the meta-layer protects.

If you catch yourself switching between terminals to cross-check feature state, that's the anti-pattern. Consolidate into one session and use `/status` for fast re-orientation.

---

## Navigation cheat sheet (for the assistant)

When the user asks **"what should I run next?"**, inspect the project state and apply this table. Always justify the recommendation in one sentence and name the alternatives.

| Project state | Recommend | Why |
|---------------|-----------|-----|
| User feels disoriented or is cold-starting a session | `/status` | One-page portfolio snapshot in <30s — dispatches `status-reporter` (Sonnet, frontmatter+grep only); replaces 10-file orientation read |
| User hand-edited `design/*.md` and wants a fresh verdict without re-entering `/design` | `/architect <slug>` | Standalone gate — dispatches only the architect agent, not the full drafting flow |
| No `docs/{slug}/` for the feature | `/map <slug>` | Need a map before anything else |
| `map.md` exists, no `review/` | `/review <slug> --as <spec>` OR `/design <slug>` | Review for non-trivial; skip to design for simple-but-non-trivial. Round-1 is parallel-independent by design |
| `map.md` + 1 review file | Read the review directly | Synthesis is overkill for one reviewer |
| `map.md` + ≥2 review files, no `synthesis.md` | `/synthesize <slug>` | Collapses into one consensus doc |
| `synthesis.md` exists, user has read it and locked a design direction (scribed to `§ User resolutions`), wants reviewers to reconcile round-1 divergence under that constraint | `/review <slug> --iterate` | Cross-informed iterate round — each reviewer sees siblings' round-1 outputs + synthesis + locked direction; re-synthesize after (MADR-H5 / ADR-019) |
| `synthesis.md` exists, no `design/` | `/design <slug>` | Translate consensus into design |
| `design/` exists with `status: DRAFT` | Approve at Gate 2 OR redirect | The user is the gate |
| `design/review.md` Verdict = NEEDS ITERATION | Re-run `/design <slug>` | Loop until READY |
| `design/review.md` Verdict = READY, no `plan/` | `/plan <slug>` | Decompose into phases; plan is auto-stamped APPROVED (ADR-018) |
| `plan/` exists, user wants hand-paced | `/implement <slug> 1` | Default — one phase at a time, gate between |
| `plan/` exists, user trusts the plan and wants end-to-end | `/implement <slug> --auto` | Runs all remaining phases; stops only at failure / stop-at-uncertainty / completion |
| `/implement` phases have landed; want to check code matches plan + design | `/verify <slug>` | Mechanical drift check — catches missing files, unused ADRs, scope drift the human can't eyeball |
| User wants to see the current architectural picture without re-running `/design` | `/diagram <slug>` | Read-only utility — collects Mermaid from map + design into `diagrams.md` in <10s; no gate |
| `/verify` returned NEEDS REVIEW or DRIFT with >2 flags | `/review <slug> --as divergent,code-reviewer` (code-reviewer agent exists in toolkit; add to role if needed) | Mechanical drift was clean; behavioural quality needs a human/reviewer pass |
| ≥2 features have `map.md`, no `_meta/map.md` | Suggest `/meta-map` IF features share touch-points | Otherwise skip meta-layer |
| `_meta/map.md` exists with convergent concerns | `/meta-design` | MADRs resolve the concerns |
| `_meta/design.md` APPROVED, ≥1 feature has `design/` | `/meta-apply` (parallel propagation) | PRIMARY path — one command replaces N serial `/design` re-runs; features needing re-architecting auto-flag MANUAL |
| `/meta-apply` emitted `MANUAL_REDESIGN_NEEDED` for a feature | `/design <slug>` | Feature needs re-architecting beyond mechanical MADR application; sub-agent correctly refused to speculate |
| `/meta-apply` architect returned `NEEDS ITERATION` for a feature | `/design <slug> --iterate` | Loop to READY; `status:` already rolled back to DRAFT |
| `_meta/design.md` APPROVED, ALL features would be MANUAL (every feature needs re-architecting) | Serial `/design <feature>` cheap-first by conflict density | `/meta-apply` would write zero edits; skip straight to serial path |
| `_meta/design.md` APPROVED but `/meta-apply` not yet attempted, some features have no `design/` | `/meta-apply` first, then `/design` for first-draft features | `/meta-apply` skips features without `design/` cleanly; handle first-drafts serially after |
| `_meta/design.md` APPROVED, all per-feature designs aligned, ≥2 features have `plan/` | `/meta-plan` | Sequence the portfolio |

For the full decision tree and per-stage tips, defer to `docs/workflows/architect/00-quickstart.md`.

### Two-path propagation (after `/meta-design` Gate 1 approval)

When the user approves a meta-design and has multiple features whose `design/` already exists, **always recommend the primary path first, falling back to serial `/design` only for the cases that need it.**

**Primary path — `/meta-apply` (parallel batch):** use whenever ≥1 feature has `design/` and Phase 7 reports mostly `NEEDS_UPDATE` rather than structural `CONFLICT`. One command dispatches N `feature-reviser` sub-agents in parallel, applies every MADR's consequences mechanically across the batch, surfaces a single compound gate, then fans out architect verdicts in parallel. Features that genuinely need re-architecting emit `MANUAL_REDESIGN_NEEDED` — no speculative edits, no writes for those features. Typical saving: ~3× tokens, ~5× wall-clock, 1 compound gate vs ~2N per-feature gates. Because the batch lands in one turn, conflict-density ordering doesn't matter — every affected feature is updated together and any MADR revision surfaces once.

**Fallback path — serial `/design <feature>` cheap-first-by-conflict-density:** use when:

- `/meta-apply` flagged `MANUAL_REDESIGN_NEEDED` for a feature — run `/design <slug>` for that one, picking up from the sub-agent's explanation.
- The feature has no `design/` yet (first-draft; `feature-reviser` cannot draft from scratch).
- Architect returned `NEEDS_ITERATION` after `/meta-apply` — loop `/design <slug> --iterate` until READY.
- You want a hand-paced, fine-grained gate on the edits per feature (occasionally worth it for a single feature with many ADRs being touched).
- Most features would flag MANUAL anyway — running `/meta-apply` first would produce near-zero writes, so save the round-trip and go serial.

In the fallback path, order serial `/design` re-runs by descending MADR-vs-existing-ADR conflict count from `/meta-design` Phase 7. Most-conflicted feature first; zero-conflict features last.

**Why cheap-first ordering in the fallback**: a MADR-vs-ADR conflict is a stop-gate (`/design` Phase 1 STOPs and forces the user to revise the MADR, override locally, or cancel). Surfacing the conflict on the most-conflicted design first means any MADR revision happens *before* the less-conflicted features re-design and inherit the (now-revised) MADR. A zero-conflict feature is an "edit-and-cite pass" — cheap, but produces no signal that could change a MADR; it goes last. This ordering is *only* relevant in the serial fallback path; `/meta-apply` obviates it.

**Fallback-to-fallback heuristic — design weight (only when Phase 7 output is absent):** approximate by descending `design/` presence: features with full `01-architecture.md` + `02-behavior.md` + `03-decisions.md` first, then features with only `design/README.md`, then features with no `design/`. More existing ADRs ≠ more conflicts, but it's a usable proxy when no conflict report exists. Phase 7 supersedes this proxy as soon as it runs.

When recommending, state which path and which heuristic you used. **Always default to `/meta-apply` first** unless the feature pattern genuinely rules it out — running it and discovering MANUAL flags is cheaper than pre-emptively going serial.

**Worked example for the pathway-explorer wave** (`normalisation`, `umap`, `color-encoding`, `theme-bundles`, `biological-workflow`):
- With `/meta-apply` (primary): one compound batch. `color-encoding` likely flags `MANUAL_REDESIGN_NEEDED` (only README exists), the other four propagate in parallel, one compound gate, architect fans out. Follow-up: `/design color-encoding` to first-draft that one.
- Without `/meta-apply` (fallback only): serial order by descending conflict count — `normalisation` (≈5) → `umap` (≈3) → `color-encoding` (≈2) → `theme-bundles` (≈2) → `biological-workflow` (0). Five gates, five re-reads of `_meta/design.md`, five main-agent dispatches.

---

## After-command tips (surface these proactively)

After every command finishes, before returning control to the user, scan `00-quickstart.md` for the matching "Tips after `<command>` finishes" block and quote the 2–3 most relevant tips. Don't repeat tips the user has already acted on in this session — but DO surface them on first encounter per stage.

Example: after `/map` completes, surface "verify a handful of file:line citations" and "skim Existing patterns to reuse" before recommending `/review` or `/design`.

Special case — **after `/meta-apply`**: always surface the MANUAL_REDESIGN_NEEDED list first (if non-empty), the ±line-count outliers second (features with >200 lines changed may have been edited too aggressively for a mechanical pass), and architect `NEEDS_ITERATION` verdicts third. Then recommend `/design <slug>` per MANUAL-flagged feature and `/design <slug> --iterate` per NEEDS-ITERATION feature. Only recommend `/meta-plan` once every feature's verdict is `READY`.

---

## Anti-patterns to refuse

- **"Skip the gate"** / **"Run the whole pipeline"** — judgment gates are the feature; refuse and explain. Note this is distinct from `/implement --auto`, which is an ADR-018-sanctioned serialisation relaxation, not a gate skip; the judgment gates (failure, stop-at-uncertainty, completion, `/verify` verdicts) all still fire inside `--auto`.
- **"Asking the user to approve plan decomposition as a gate"** — the `/design` architect verdict is the real gate. Plan is derivative; re-gating at plan time re-litigates a decision already made. Surface the phase list as a glance-check, not an approval request.
- **"Waiting between every phase in `/implement --auto`"** — defeats the purpose of the flag. In `--auto` the inter-phase stop conditions are failure / uncertainty / completion; anything else is a halt that shouldn't exist.
- **"Re-grep the codebase from a downstream agent"** — only `mapper` reads code; downstream agents read `map.md`. If `map.md` is incomplete, re-run `/map`.
- **"Cross-inform reviewers on round-1"** — defeats the purpose of the parallel-independent regime; divergence is the signal. Cross-informing belongs to `/review --iterate` only, never to round-1.
- **"Run `/review --iterate` without a prior round"** — there is no "round-1 iterate"; the iterate regime requires something to reconcile against. The command rejects this explicitly.
- **"Silently rewrite a round-1 position in a round-2 review"** — iterate reviewers must tag every prior position UPDATED / STOOD BY / RETIRED. Silent rewrites destroy the reconciliation epistemics that justify running iterate in the first place.
- **"Overwrite a prior review or synthesis without archiving"** — MADR-H5 requires `.history/round-N/` archive before any replacement. If the archive slot already contains a file, abort — do not overwrite; a prior move was interrupted or the round-numbering is wrong.
- **"Edit `_meta/design.md` directly to resolve a conflict"** — only `/meta-design` writes there. Resolve via re-dispatch.
- **"Edit a per-feature design after architect verdict READY without re-dispatching architect"** — the verdict no longer applies; loop.
- **"Start re-design with the easiest / zero-conflict feature after /meta-design Gate 1 approval"** — violates cheap-first by conflict density; the zero-conflict feature is an edit-and-cite pass that produces no signal, while the most-conflicted feature is where MADR revisions need to happen first. Always go most-conflicted first; explain the rationale rather than silently letting the user start with the easy one.
- **"Recommend re-design order purely by design weight after Phase 7 has produced a conflict report"** — design weight is the *fallback* heuristic for when conflict counts are unknown; conflict density is primary. A heavy design that fully aligns with MADRs has zero conflicts and goes last.
- **"Run N serial `/design` re-runs after `/meta-design` approval when `/meta-apply` would work"** — 5× wall-clock and ~3× tokens for the same result. Use `/meta-apply` first; drop to serial `/design` only for features the sub-agent flagged `MANUAL_REDESIGN_NEEDED` or for features with no `design/` yet.
- **"Interpret `MANUAL_REDESIGN_NEEDED` as a failure or a defect"** — it's a correctness signal. The sub-agent refused to speculate on a feature whose MADR consequences can't be applied mechanically. Treat it as a routing decision ("this feature belongs in the `/design` lane, not `/meta-apply`") and run `/design <slug>` for that specific feature.
- **"Force `/meta-apply` on a feature that has no `design/` directory"** — the command already skips it cleanly with a recommendation; don't add flags or redirects to push it through. First-drafts belong to `/design`.
- **"Tell a `feature-reviser` to try again after `MANUAL_REDESIGN_NEEDED`"** — the sub-agent saw the same inputs twice; re-prompting won't change the answer. Escalate to `/design <slug>` where the main agent has the capability to re-architect.
- **"Ignore `NEEDS_ITERATION` from `/meta-apply`'s architect phase and proceed to `/plan`"** — the command already rolled `status:` back to `DRAFT` for that feature; `/plan` will refuse. Loop with `/design <slug> --iterate` first.
- **"Ship without `/verify` after all `/implement` phases land"** — the mechanical drift check is the bug class a human can't eyeball across thousands of lines. Missing files, unreached ADRs, and scope drift are the three failure modes `/verify` is purpose-built to surface. Running it takes <1 minute and turns "did we build what we said?" from a judgment call into a checked fact.
- **"Treat `/verify` as a quality audit"** — it is not. `/verify` is mechanical (anchors, files, mtimes). Behavioural correctness, architecture quality, and accidental complexity require human review or `/review <slug> --as divergent,code-reviewer`. `/verify` surfaces what it can mechanically check AND flags the rest as `NEEDS_HUMAN_REVIEW` — don't interpret CLEAN as "the code is good."
