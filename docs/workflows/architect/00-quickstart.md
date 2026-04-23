---
parent: ./README.md
view: how-to
audience: human + assistant
---

# 00 — Quickstart

Read this first. It's the operating manual for the harness: which sequence to run for what kind of work, when the meta-layer is worth invoking, and what to do at every gate.

If you're an assistant helping a user navigate the harness, this doc is the source of truth for cadence questions — defer to it over your training intuitions.

---

## Commands at a glance

**Per-feature pipeline** (one feature):

```
/map <slug>                       — codebase cartography → docs/{slug}/map.md
/review <slug> --as <spec>        — N reviewer lenses in parallel, ROUND-1 regime → review/*.md
/review <slug> --iterate --as <spec>  — cross-informed ITERATE round; reconciles after round-1 divergence
                                        → review/*.md (prior round archived to review/.history/round-N/)
/synthesize <slug>                — collapse reviews → synthesis.md (≥2 reviewers required;
                                        prior synthesis archived to synthesis/.history/round-N.md on re-run)
/design <slug>                    — staged drafting + architect gate → design/*.md
/architect <slug>                 — standalone architect gate (re-run after hand-edits; no /design re-entry)
/plan <slug>                      — phased decomposition → plan/phase-NN.md
/implement <slug> <N>             — one phase, hand-paced (default) → code
/implement <slug> --auto          — end-to-end: run all remaining phases, stops only at failure/uncertainty/completion
/verify <slug> [phases]           — mechanical drift check: plan + design vs code → verify.md
/diagram <slug>                   — render Mermaid from design/ into docs/{slug}/diagrams.md (read-only, <10s)
```

**Portfolio state** (any time):

```
/status [slug,...]                — one-page portfolio snapshot (<30s) → chat (no file writes)
```

**Meta-layer** (cross-feature, optional):

```
/meta-map                         — portfolio inventory → docs/_meta/map.md
/meta-design                      — inter-feature ADRs (MADRs) → _meta/design.md
/meta-apply                       — propagate ACCEPTED MADRs to every in-flight design/ in parallel
                                    (dispatches N feature-revisers + N architects; one compound gate)
/meta-plan                        — sequence portfolio phases → _meta/plan.md
```

Every command stops at a human gate. There is no autopilot — that is the feature.

---

## Decision tree: which sequence do I run?

```
Is this one-line, throwaway, or already-specified work?
├── YES → Don't use the harness. Just edit.
└── NO ↓

How many features are in flight at the same time?
├── ONE → Per-feature pipeline only:
│         /map → [/review → /synthesize] → /design → /plan → /implement
│         (review + synthesize optional for small features)
│
└── TWO+ ↓
    Do the in-flight features touch overlapping files / share concepts?
    ├── NO → Run per-feature pipelines independently. Skip meta-layer.
    │
    └── YES → Meta-layer wraps the per-feature work:
              [per-feature /map for each] → /meta-map →
              [optional per-feature /review + /synthesize] →
              /meta-design → [per-feature /design (inherits MADRs)] →
              /meta-plan → [per-feature /plan] → /implement (in portfolio order)
```

The "in flight" test is empirical: if you the architect feel "swirled" reading two synthesis docs and mentally cross-referencing them, that's the felt-sense signal that you needed the meta-layer two stages ago. See [05-meta-approach.md](./05-meta-approach.md) for the full rationale.

---

## Per-feature flow — the canonical run

The simplest case: one non-trivial feature, you want it well-designed before you write code.

### 1. `/map <slug> "<one-paragraph description>"`

Produces `docs/{slug}/map.md` — files touched, entry points, data flow, existing patterns to reuse, open questions.

**Tips after `/map` finishes:**
- Read the **Open questions** section first. These often need your input before any reviewer can do useful work.
- Skim the **Existing patterns to reuse** section — if mapper missed an obvious pattern, regenerate or hand-edit `map.md` before reviewers consume it. A wrong map poisons everything downstream.
- Verify a handful of `file:line` citations. If they're wrong, regenerate. Mapper makes mistakes.

### 2. `/review <slug> --as <spec>` (optional but usually worth it)

Six reviewer lenses available: `bioinf`, `wetlab`, `graphic`, `stat`, `divergent`, `ml`. Pick what's relevant.

```
/review <slug> --as all                          — full panel (large/cross-cutting features)
/review <slug> --as all --but stat,wetlab        — exclusion syntax
/review <slug> --as bioinf,ml                    — targeted subset
/review <slug>                                   — interactive; assistant asks
```

**Tips after `/review` finishes:**
- Reviewers ran in parallel — wall-clock is one reviewer, not N.
- Each reviewer wrote `docs/{slug}/review/{name}.md`. Skim each Headline and P0 list; that's where the high-leverage items are.
- A reviewer flagging "Open questions for mapper" means `map.md` is incomplete — re-run `/map` (or hand-edit) before continuing.
- **Two regimes, distinct purposes.** Round-1 (`/review --as ...`) is parallel-independent — reviewers are blind to each other. That's deliberate: divergence is the signal. Iterate (`/review --iterate --as ...`) is cross-informed — each reviewer sees the others' round-1 outputs, the synthesis, and your locked direction (via `synthesis.md § User resolutions`). Use iterate when you've read the round-1 synthesis, locked a design direction, and want reviewers to reconcile under that constraint. See MADR-H5 for details.
- **Iterate requires position-shift tagging.** Every iterate-round reviewer file opens with a `## Round N iterate summary` listing which prior positions were UPDATED / STOOD BY / RETIRED. This preserves the epistemics — round-2 that silently erases round-1 defeats the point.
- **Prior rounds are archived, not deleted.** `review/.history/round-N/<name>.md` and `synthesis/.history/round-N.md` keep the trail. Safe to iterate as many times as you need; nothing destroys.

### 3. `/synthesize <slug>` (only if ≥2 reviewers)

Collapses reviews into a consensus doc with Convergent concerns, Divergent opinions (preserved), Open questions, P0/P1/P2 recommendations.

**Tips after `/synthesize` finishes:**
- The **Headline** is the one sentence to remember. The **Divergent opinions** section is the second-most-important — synth deliberately does not flatten disagreement, and those disagreements are usually load-bearing.
- If everything is P0, synth has failed; ask it to re-prioritise.
- One reviewer? Skip `/synthesize` and read the single review directly. The harness will refuse anyway.

### 4. `/design <slug>`

Staged drafting:

- **Phase 3a → Gate 1:** README only (scope, acceptance criteria, ADR headlines, conditional-doc decisions). User approves scope.
- **Phase 3b → Gate 2:** Full bundle (`01-architecture.md`, `02-behavior.md`, `03-decisions.md` + conditional). User approves bundle.
- **Phase 5:** `architect` agent gates with Verdict (READY / NEEDS ITERATION / NEEDS DISCUSSION).

**Tips after `/design`:**
- **Phase 3a is the cheap correction point.** ~150 lines vs ~1500. Push back on a wrong scope here, not after the ADRs are drafted. Gate-1 now carries two coarse Mermaid diagrams (module sketch + data-flow sketch) — these are scope-level and non-binding; detailed diagrams land in Phase 3b.
- **Phase 3b's `01-architecture.md § Delta view`** is the annotated current-vs-proposed figure — reviewers count green nodes/edges to detect over-scoping. Revise here, not after architect.
- **Phase 4 (full bundle approval)** is for ADR-content review. If ADR-N's Decision is wrong, revise here — architect won't catch a wrong-but-internally-consistent design.
- **NEEDS ITERATION verdict** means the loop continues — the design has consistency or completeness gaps. It's not a failure; it's the harness working.
- **NEEDS DISCUSSION verdict** means architect saw something that needs your judgment, not a fix. Read `design/review.md` and reply.
- **READY verdict** unlocks `/plan`. Don't hand-edit design docs after READY without re-dispatching architect.

### 5. `/plan <slug>`

Decomposes design into ordered, testable phases.

**Tips after `/plan`:**
- Plan READMEs are stamped `status: APPROVED` at draft time (the architect verdict `READY` on the design is the real commitment; re-gating at plan time would re-litigate it). Phase 4's prompt is a *glance-check*, not a status gate — respond only if a phase boundary is obviously wrong.
- Each `phase-NN.md` lists Files-to-Create and Files-to-Modify. Skim and confirm the boundary makes sense.
- Phases follow data → logic → storage → API → integration → hardening. If they don't, push back *here* (cheap to re-plan; expensive to unwind after implementation).
- A phase touching >5 files is a smell — ask the assistant to split it.

### 6. `/implement <slug> <N>` or `/implement <slug> --auto`

Default — `<slug> <N>` — executes exactly one phase, runs tests, stops. Hand-paced; you confirm before the next phase.

Opt-in — `<slug> --auto` — runs all remaining phases end-to-end, stopping only at (a) phase failure, (b) stop-at-uncertainty, (c) all phases complete. Per-phase output is a one-liner (`Phase 3 ✅ — 4 files, 87 lines`); full Phase-4 verification block only fires when the loop stops. Use this when you trust the plan and don't want to hand-pace 12 phases.

**Tips after `/implement`:**
- The meta-architect's gates are at MADRs, `/design` Phase 4 ADRs, `/meta-plan`, and stop-at-uncertainty — not at plan decomposition or per-phase cadence. `--auto` is the mode that matches that lens.
- Verify the phase checklist after each `/implement`. In `--auto` mode, the checklist runs on every phase automatically and the loop halts on any failure.
- If `/implement` halted with "stop at uncertainty", read what it asked — don't override blindly. This rule applies equally in `--auto`.
- Phase failures often mean the design assumed something the codebase doesn't actually have. Go back to `/map` (often) or `/design` (sometimes), not deeper into `/implement`.
- After an `--auto` completion, the follow-up is always `/verify <slug>` — the mechanical drift check catches what silent per-phase runs could hide.

### 7. `/verify <slug> [phases]` — mechanical drift check

Compare `plan/phase-NN.md` (Files-to-Create + Files-to-Modify) and `design/03-decisions.md` (ADR Decisions) against the actual code. Produces `docs/{slug}/verify.md` with a Verdict: `CLEAN | INCOMPLETE | DRIFT | NEEDS REVIEW`. Main-agent only; no sub-agent. Targets <1 minute.

**Tips after `/verify`:**
- **INCOMPLETE** = implementation didn't finish. Named phases still have `create` files missing or `modify` files un-touched. Continue `/implement <slug> N` for each missing phase.
- **DRIFT** = either files were modified that weren't in the plan (scope creep), or at least one ADR's grep-anchor is missing from code (the design said X but the code doesn't mention X). Hand-read the flagged items; decide whether to backfill the plan, revert the change, or promote a drift to a real ADR.
- **NEEDS REVIEW** = mechanical checks are clean but >30% of ADRs are behavioural-only (can't be verified by grep). The design is heavy on judgment; dispatch `/review <slug> --as divergent,code-reviewer` or hand-read the flagged ADRs.
- **CLEAN** ≠ "the code is good." It means every mechanical check passed. Code quality, architecture soundness, and accidental complexity still need human or reviewer-panel eyes. `/verify` is a floor, not a ceiling.
- Run `/verify <slug> 1-3` to scope to specific phase ranges during a long implement wave — useful for early-warning drift detection without waiting for all phases to land.

### `/status` — portfolio state snapshot (any time)

Produces a one-page table summarising every feature's stage, architect verdict, plan status, and any suspected cross-feature conflicts. Reads artifacts only (frontmatter + grep); writes nothing. Target wall-clock: <30 seconds.

**Tips after `/status`:**
- Use before any command when you're unsure where the portfolio is, after a context switch, or at the start of a new session. It's cheaper than opening 10 files.
- The "Outstanding" column names the one blocker per feature — that's usually what to act on next.
- "Open cross-feature conflicts" is best-effort, heuristic-grade. A flagged conflict warrants opening the cited `03-decisions.md` sections yourself, not automatic action.
- `/status` does not gate or write. If you want durable state, commit the output to `docs/_meta/STATUS.md` yourself.

### `/diagram <slug>` — render current architectural picture

Collects Mermaid blocks from `map.md` and `design/*.md` into one consolidated `docs/{slug}/diagrams.md`. Read-only utility; no gate; no sub-agent; target wall-clock <10s. Parallel to `/status` (portfolio-level) but scoped to one feature.

**Tips after `/diagram`:**
- Use mid-design-conversation when you want to see "the shape as it stands" without triggering `/design` re-entry.
- The "Missing" section tells you which diagrams *haven't* been drafted yet — often a useful to-do list before the next gate.
- `diagrams.md` is auto-generated. Hand-edits are overwritten on the next run; fix the source design docs instead.

### `/architect <slug>` — standalone architect gate

Re-runs the architect agent on `docs/{slug}/design/*.md` without re-entering `/design`. Use after you've hand-edited a design doc, after `/meta-design` revised a MADR this feature inherits, or after a peer feature's APPROVED design changed. Produces a fresh `design/review.md` with a Verdict.

**Tips after `/architect`:**
- Unlike `/design` Phase 6, this command does NOT auto-loop on NEEDS_ITERATION. It's a single-shot gate; you iterate by re-invoking.
- If design frontmatter says `DRAFT` and the verdict is READY, the command asks before flipping to APPROVED — keep intent explicit.
- Use this instead of re-running `/design` when the change is small and already in the docs. `/design` re-drafts; `/architect` just re-gates.

---

## Meta flow — when ≥2 features are in flight

This is the prototypical sequence when you have two or more features whose designs would otherwise drift apart.

### Prototypical case — wave of 5 features

Suppose pathway-explorer has five features in flight: `normalisation`, `umap`, `color-encoding`, `theme-bundles`, `biological-workflow`. They share files (`similarity.py`, `data_loader.py`, `main.py`, `html_generator.py`) and concepts (what "score" means, what "in-set gene" means, etc.).

Without the meta-layer, each `/design` will pick its own answer to "what is a score?" — and the answers will silently contradict each other. The architect catches the contradictions in their head while reading 5 synthesis docs (the swirl).

With the meta-layer, here's the cadence:

```
# Bootstrap — produce per-feature maps
/map normalisation
/map umap
/map color-encoding
/map theme-bundles
/map biological-workflow

# Aggregate
/meta-map normalisation,umap,color-encoding,theme-bundles,biological-workflow
   → docs/_meta/map.md surfaces:
       - Shared touch-points (similarity.py touched by 4 features, etc.)
       - Convergent concerns (the same defect surfaced by ≥2 features)
       - Inter-feature dependencies
       - OQ-M-* questions

# Per-feature review (where useful)
/review theme-bundles --as bioinf,ml,divergent
/synthesize theme-bundles
# … repeat for features that need a panel

# Inter-feature ADRs — the heart of the meta-layer
/meta-design
   → docs/_meta/design.md drafts MADR-001..N covering every convergent concern
   → Gate 1: you approve / revise / reject MADRs
   → Status flips PROPOSED → ACCEPTED on approval

# Propagation — ONE command, parallel sub-agents, one compound gate (PRIMARY path)
/meta-apply
   → Phase 2: main agent allocates deferred-ID blocks per feature (serialised, race-free)
   → Phase 3: dispatches N feature-revisers in parallel (one message, N Agent calls)
       - Each reads _meta/design.md + its own docs/{slug}/ only
       - Each applies MADR consequences mechanically (edit-and-cite; no re-architecting)
       - Each appends its feature-local deferrals to _meta/deferred.md via its exclusive ID block
       - Features needing re-architecting emit MANUAL_REDESIGN_NEEDED (no writes) — clean skip
   → Phase 4: compound diff report (per-feature table + highlights + MANUAL flags quoted)
   → Gate 1: approve all / reject some / redo one / drill into a feature's diff
   → Phase 6: on approval, status DRAFT → APPROVED; dispatch architect subagents in parallel
   → Phase 7: collect per-feature verdicts; NEEDS_ITERATION rolls status back to DRAFT

# Fallback — manual per-feature for MANUAL flags, first-drafts, or NEEDS_ITERATION loops
/design color-encoding         # MANUAL: only README — needs first-draft of 03-decisions.md
/design biological-workflow    # architect verdict was NEEDS_ITERATION after /meta-apply — loop
# Any /design run in this fallback position still:
#   → Phase 1 reads _meta/design.md
#   → If a local decision conflicts with a MADR, /design STOPS and asks:
#       (a) re-run /meta-design to revise the MADR
#       (b) override locally (recorded in 03-decisions.md § Meta overrides)
#       (c) cancel

# Cross-feature sequencing
/meta-plan
   → docs/_meta/plan.md sequences phases:
       - P0: normalisation phase-01 (introduces within_method_z; unblocks 3 downstream)
       - P1: umap phase-01 + normalisation phase-02 (parallel, no shared files)
       - P2: theme-bundles phase-01 + color-encoding phase-01
       - P3: biological-workflow phase-01
   → Collision check flags any same-file overlap as a blocker

# Per-feature plan (each honours its portfolio slot)
/plan normalisation
/plan umap
…

# Implement in portfolio order
/implement normalisation 1   # P0
/implement normalisation 2   # P1
/implement umap 1            # P1 (parallel-safe with normalisation 2)
…
```

**Tips after `/meta-map`:**
- The **Convergent concerns** section is the priority list for `/meta-design`. If it's empty, you don't need the meta-layer for this wave.
- The **Shared touch-points** matrix tells you which files will get hit by multiple features. Mentally rehearse the merge order.
- **OQ-M-N** items are cross-feature questions that nobody can answer alone. Resolve them before `/meta-design`.

**Tips after `/meta-design`:**
- Every MADR has a "Consequences per feature" line — that's what each feature's `/design` must honour.
- After Gate 1 approval, Phase 6.5 appends time-boxed items from `Deliberate meta-level rejections § Time-boxed sub-items`, deferred `OQ-M-N` items, and future-round MADR Consequences to `_meta/deferred.md` (rows tagged `Source feature: _meta`). Verify the row list looks right before continuing — `/meta-plan` reads them as Parking items.
- Then the assistant lists "downstream impact" — features whose existing design conflicts with a new MADR, with a **conflict count per feature**. Two propagation paths are available; the primary path makes conflict ordering mostly moot.
- **Primary path: `/meta-apply`** (parallel batch). Dispatches N `feature-reviser` sub-agents in parallel; each applies its feature's MADR consequences mechanically. One compound gate on the batch of edits; architect gate fans out per feature. Features whose MADR consequences would need re-architecting emit `MANUAL_REDESIGN_NEEDED` and are cleanly skipped — no speculative writes. Typical saving versus serial `/design`: ~3× tokens, ~5× wall-clock, 1 compound gate versus ~2N per-feature gates. Because the compound gate surfaces all edits in one turn, conflict-density ordering no longer matters: every affected feature is updated together, and MADR revisions surfaced by any MANUAL flag happen once, before any re-propagation.
- **Fallback path: serial `/design` cheap-first-by-conflict-density.** Use when:
  - `/meta-apply` flagged `MANUAL_REDESIGN_NEEDED` for a feature — run `/design <slug>` for that feature, picking up from the MANUAL's explanation.
  - The feature has no `design/` yet (first-draft; `feature-reviser` cannot draft from scratch).
  - Architect returned `NEEDS ITERATION` after `/meta-apply` — loop `/design <slug> --iterate` until READY.
  - You want a hand-paced, fine-grained gate on the edits instead of a batch compound review.

  In the fallback, order serial `/design` re-runs by descending MADR-vs-existing-ADR conflict count (from `/meta-design` Phase 7). Most-conflicted feature first; zero-conflict features last.
- **Why cheap-first ordering in the fallback**: a MADR-vs-ADR conflict is a stop-gate — `/design` Phase 1 STOPs and forces you to revise the MADR, override locally, or cancel. Surfacing the conflict on the most-conflicted design first means any MADR revision happens *before* the less-conflicted features re-design and inherit the (now-revised) MADR. A feature with zero conflicts is an "edit-and-cite pass" — cheap, but produces no signal that could change a MADR; it goes last. This ordering is *only* relevant in the fallback path; `/meta-apply` obviates it by batching.
- **Fallback-to-fallback heuristic (pre-Phase-7, conflict counts unknown):** approximate by descending **design weight** — features with full `design/03-decisions.md` first, then features with only `design/README.md`, then features with no `design/`. More existing ADRs ≠ more conflicts, but it's a proxy when you don't yet have the conflict report. Phase 7 supersedes this proxy as soon as it runs.
- **Worked example (pathway-explorer wave)**: Phase 7 surfaced ~5 conflicts in `normalisation`, ~3 in `umap`, ~2 each in `color-encoding` + `theme-bundles` (both had only README), ~0 in `biological-workflow`. With `/meta-apply`: one compound batch — `color-encoding` likely flagged MANUAL (first-draft), the other four propagated in parallel, one compound gate, architect fans out. Without `/meta-apply` (the old flow): serial order `normalisation → umap → color-encoding → theme-bundles → biological-workflow`, five gates, five re-reads of `_meta/design.md`.
- "Deliberate meta-level rejections" splits into permanent rejections (stay in `_meta/design.md`) and time-boxed sub-items (also appended to `_meta/deferred.md`). Permanent rejections without a time horizon stay only in design.md — don't overload deferred.md with non-actionable noise.

**Tips after `/meta-apply`:**
- Read the **MANUAL_REDESIGN_NEEDED** list first — those are the features where a sub-agent refused to speculate, and the user's cue to run `/design <slug>` manually for each. MANUAL is a routing decision, not a failure.
- Skim the compound diff table for ±line counts per feature. A feature with >200 lines changed probably warranted `/design` rather than `/meta-apply` — consider re-running that feature via `/design` for a careful hand pass if the edits feel architecturally aggressive.
- Architect `NEEDS ITERATION` on a feature means the MADR application produced an internal inconsistency. The command has already rolled that feature's `status:` back to `DRAFT`. Re-run `/design <slug> --iterate` to loop architect until READY.
- `_meta/deferred.md` has N new rows appended, one contiguous ID block per sub-agent. Verify no block was exhausted (sub-agents were told to emit MANUAL rather than steal IDs, but check for gaps or overlaps).
- Frontmatter: after Phase 7, `READY` features have `status: APPROVED`, `NEEDS_*` features are `DRAFT`, MANUAL features are unchanged from pre-run. `/plan <slug>` requires `APPROVED`.

**Tips after `/meta-plan`:**
- Collisions are blockers, not warnings. Resolve before `/implement`.
- The **Parking items** section is `_meta/deferred.md` projected forward — items the architect parked across features.
- Re-run `/meta-plan` whenever a per-feature `/plan` changes the phase boundaries.

---

## Re-running meta-* — five triggers

1. **A new feature enters the wave** → `/meta-map` (and probably `/meta-design`)
2. **A per-feature `/design` flags a MADR conflict** → `/meta-design` with the user's resolution
3. **Implementation reveals a wrong assumption in the meta-plan** → `/meta-plan`
4. **A new portfolio wave starts** → fresh `/meta-map`; archive the old `_meta/` to `docs/_meta/.archive/{date}/` (manual)
5. **`/meta-apply` emitted `MANUAL_REDESIGN_NEEDED` for one or more features** → NOT a meta re-run. Run `/design <slug>` per flagged feature. If the flag pattern is systemic (the same MADR keeps getting flagged across features), that's a meta re-run signal — re-run `/meta-design` to revise *that specific MADR*.

---

## Common gotchas

| Symptom | Cause | Fix |
|---------|-------|-----|
| `/synthesize` refuses to run | <2 reviewer files | Read the single review directly, or run `/review` for another lens |
| Reviewer's output references wrong file:line | Mapper had stale or incomplete `map.md` | Re-run `/map` (or hand-edit), then re-dispatch the reviewer |
| `/design` Gate 2 surfaces a wrong scope decision | You missed it at Gate 1 | Cheap to revise now; very expensive to revise after the architect verdict |
| Architect verdict is `NEEDS ITERATION` repeatedly | Persistent consistency gaps in design | Read `design/review.md` line by line; fix specifically what it lists |
| `/plan` output ignores `_meta/plan.md` slot | `_meta/plan.md` was added after `/plan` ran | Re-run `/plan` |
| `/design` STOPs with "MADR-N conflicts" | A MADR in `_meta/design.md` requires X; this feature wants Y | Choose: re-run `/meta-design`, override locally, or cancel |
| `/meta-plan` flags collisions in every phase | Two features touch the same file in their phase-01s | Sequence one feature's phase-01 to a later portfolio phase, or split the file change |
| Meta-architect references source code | Wrong tool — meta-architect should only read `docs/` | Re-dispatch and tell it to use artifacts only |
| Per-feature `/design` doesn't read `_meta/design.md` | `_meta/` is in the wrong place | Must live at `docs/_meta/`, sibling to per-feature dirs |
| `/verify` reports INCOMPLETE but you think you finished | Plan's Files-to-Modify list includes a file that got renamed/split during implementation | Either backfill the plan with the new file names, or accept the flag as a known-stale row — don't `touch` the file just to silence the mtime check (that hides real state) |
| `/verify` reports DRIFT on a helper you added intentionally | Small helper got added outside the plan's scope | Either backfill a phase (honest: captures what actually happened) or accept the flag as incidental; don't revert just to silence the warning |
| `/verify` reports NEEDS REVIEW with most ADRs behavioural | Design is heavy on judgment (biology, UX, statistics) that grep can't verify | Dispatch `/review <slug> --as divergent,code-reviewer` or pair with a domain expert; this is expected for domain-heavy features |
| `/verify` CLEAN but code feels off | Mechanical pass is a floor, not a ceiling | CLEAN only means "every anchor present"; it never claims correctness or quality. Pair with a reviewer panel pass on the code |

---

## When NOT to use the harness

- One-line fixes, typo corrections, trivial tweaks
- Throwaway scripts or exploratory data analysis
- Cases where the design is already locked (just implementing an approved spec) — go straight to code

## When NOT to use the meta-layer

- Only one feature in flight
- Features have zero shared touch-points (e.g., a backend refactor and an unrelated UI tweak)
- Quick prototypes where cross-feature sequencing doesn't matter

The cost of skipping meta when it's not needed is zero. The cost of skipping it when it IS needed is the swirl — concerns that re-surface at integration time, when fixes are 10× more expensive.

---

## What to ask the assistant for, by stage

The assistant has a system-prompt-loaded skill (`architecture-first-dev`) that knows this cadence. When in doubt, you can ask:

- **"What stage am I at?"** — Assistant inspects `docs/{slug}/` and `docs/_meta/` and reports.
- **"What should I run next?"** — Assistant reads the latest artifact and recommends the next command.
- **"Is the meta-layer worth invoking here?"** — Assistant checks how many features have `map.md`, whether they share touch-points, and recommends.
- **"Why did `/design` stop?"** — Assistant explains MADR conflicts, Gate failures, or architect verdicts.
- **"What did I park, and where?"** — Assistant reads `_meta/deferred.md` and surfaces relevant items.
- **"What's blocking `/plan` for feature X?"** — Assistant reads `design/review.md` Verdict and reports.
- **"Show me the architecture as it stands"** — `/diagram <slug>`. Collects every Mermaid block from the feature's artifacts into one file.

The harness is meant to be navigated, not memorised. Ask.

---

## See also

- [README.md](./README.md) — index + activation
- [01-architecture.md](./01-architecture.md) — what files exist, what binds them
- [02-behavior.md](./02-behavior.md) — per-stage behaviour with sequence diagrams
- [03-decisions.md](./03-decisions.md) — why the harness has this shape
- [04-extending.md](./04-extending.md) — adding new reviewers / commands / roles
- [05-meta-approach.md](./05-meta-approach.md) — the meta-layer in depth, including the worked pathway-explorer example
