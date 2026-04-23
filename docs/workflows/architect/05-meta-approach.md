---
parent: ./README.md
view: how-to
---

# Meta-approach — designing across interoperating features

The per-feature pipeline (`/map → /review → /synthesize → /design → /plan → /implement`) works well for one feature at a time. When two or more features are in flight simultaneously and they touch overlapping surfaces, each feature's `/design` makes locally-optimal choices that can quietly contradict its siblings. The meta-layer is the optional tier above per-feature designs that holds the *continuities* across them.

---

## Why this layer exists

Empirical pattern (observed running the harness on five interoperating pathway-explorer features — `normalisation`, `umap`, `color-encoding`, `theme-bundles`, `biological-workflow`):

The same underlying defect surfaces independently across multiple features' reviews — same concern refracted through different lenses. Examples:

- **Score incommensurability** flagged by `normalisation`, `theme-bundles/bioinf`, `color-encoding/graphic+ml`
- **Leading-edge vs full-set gene membership** flagged by `biological-workflow`, `theme-bundles`
- **Shared UMAP canvas implying false metric equivalence** flagged by `color-encoding/ml`, `theme-bundles/ml`, `umap`
- **Static-HTML payload budget** flagged by `biological-workflow/graphic+divergent`
- **Null-model absence for "cluster = convergence" claims** flagged by `theme-bundles/ml`, `biological-workflow/divergent`
- **Column rename `signed_sig → within_method_z`** blocking downstream consumers in `normalisation`, `theme-bundles`

These are not five feature-local critiques — they are one system-level critique refracted through five lenses. The architect ends up doing the refraction-inversion manually, by reading five independent synthesis docs and cross-referencing in their head. That's "the swirl."

The meta-layer formalises the cross-reference. Per-feature designs stay local; the meta-layer makes the inter-feature agreements explicit as MADRs (meta-ADRs) that each per-feature design *inherits*.

### Pirsig framing

Per-feature reviewers see **chunks** (one feature at a time, one lens at a time). The meta-architect sees **continuities** (what holds across features). The harness's job is to make traversal between the two cheap. Without the meta-layer, the architect is the only place continuities live — which is the swirl.

### Frederik framing

- **Reviewers (bioinf, ml, stat, graphic, wetlab, divergent)** = engineers. Deep specialists, one lens.
- **Pipeline agents (mapper, synth, architect)** = per-feature discipline.
- **User** = architect of the architects.
- **meta-architect** = the generalist agent that reads across features, drafts integrative ADRs, and sequences cross-feature implementation. The bridge.

---

## When to invoke

**Use the meta-layer when:**

- ≥2 features are in flight at the same time
- The features touch overlapping files (the shared touch-point matrix has off-diagonal entries)
- OR the architect feels "swirled" reading the per-feature syntheses — that's the felt-sense signal
- A feature's `/design` would otherwise make a decision (e.g., "we'll use within-method z-scores") that another feature would also have to make independently

**Skip the meta-layer when:**

- Only one feature is in flight
- Features have zero shared touch-points (e.g., a backend refactor and an unrelated UI tweak)
- Quick prototypes where sequencing across features is not load-bearing

The cost is low when the meta-layer is unnecessary: no `_meta/` directory means per-feature `/design` and `/plan` behave exactly as before.

---

## Team of agents (with meta-layer)

Unchanged:

- 6 reviewers (`bioinf`, `wetlab`, `graphic`, `stat`, `divergent`, `ml`) — per-feature lenses
- 3 pipeline agents (`mapper`, `synth`, `architect`) — per-feature discipline

Added:

- `meta-architect` — pure orchestrator. Reads across features. Writes only `docs/_meta/*`. Never reads source code; never writes inside per-feature directories. Never gates anything; the user gates instead.

- `feature-reviser` — narrow-scope propagator. Dispatched in parallel by `/meta-apply`, one per feature. Applies accepted MADRs from `_meta/design.md` to exactly one feature's existing `design/` tree (edit-and-cite only; never re-architects). Appends feature-local deferrals to `_meta/deferred.md` using a pre-assigned ID block. Emits `MANUAL_REDESIGN_NEEDED` when a MADR's consequence can't be applied mechanically (first-drafts, structural re-thinking) — those features fall back to serial `/design`.

Total: 11 agents, 10 commands.

---

## Artifact layout

```
docs/
├── _meta/                         (meta-layer — sibling of features)
│   ├── map.md                     produced by /meta-map
│   ├── design.md                  produced by /meta-design
│   ├── plan.md                    produced by /meta-plan
│   └── deferred.md                appended by per-feature /design Phase 3b AND /meta-design Phase 6.5 AND /meta-apply feature-reviser sub-agents (using pre-allocated ID blocks)
│
├── normalisation/                 (per-feature)
│   ├── map.md
│   ├── review/
│   ├── synthesis.md
│   ├── design/
│   └── plan/
├── umap/                          (per-feature)
│   └── ...
└── theme-bundles/                 (per-feature)
    └── ...
```

The underscore prefix on `_meta/` signals "framework-level, not a feature." It sorts to the top in alphabetical listings, which is also useful.

`_meta/` lives in the project's `docs/` tree, NOT inside the toolkit. It's project-specific runtime state, like the per-feature directories.

---

## Commands and the canonical cadence

```
A. Bootstrap (once per wave)
   1. /map f1                                          (per feature)
   2. /map f2
   ...
   3. /meta-map                                        (aggregate all)

B. Review (optional, per feature)
   4. /review fX --as <list>                           (per feature, as needed)
   5. /synthesize fX                                   (per feature, conditional on ≥2 reviewers)

C. Design envelope
   6. /meta-design                                     (produces MADRs — inter-feature ADRs)

D. Per-feature design propagation — PARALLEL (primary path after /meta-design Gate 1)
   7. /meta-apply                                     (propagates every accepted MADR to every in-flight design/
                                                       in one parallel batch; compound diff gate; architect gate
                                                       fans out in parallel; features needing re-architecting
                                                       auto-flag MANUAL and fall back to step 7' below)

D'. Per-feature design — SERIAL (fallback path, used for MANUAL flags + first-drafts + deep reworks)
   7'. /design fX                                      (per feature; inherits MADRs)
       → if /design surfaces a MADR conflict, goto step 6 (re-run /meta-design with user resolution)
       → if /meta-apply flagged MANUAL_REDESIGN_NEEDED, this is the fallback

E. Plan envelope
   8. /meta-plan                                       (sequences phases across features)

F. Per-feature plan
   9. /plan fX                                         (per feature; honours portfolio slot)

G. Implementation (in portfolio order)
   10. /implement fX N                                 (per feature, one phase at a time)

H. Close-out (optional)
   11. /meta-map (regenerate)                          (snapshot of shipped state)
```

Every meta-stage is human-gated, just like every per-feature stage. There is no "run the whole pipeline" macro for the meta-layer either.

---

## Propagation phase — the team-parallel step

After `/meta-design` Gate 1 approves the MADRs, those decisions have to land in per-feature designs. There are two paths, and they compose rather than compete.

### Team-parallel (primary) — `/meta-apply`

`/meta-apply` dispatches N `feature-reviser` sub-agents in parallel, one per feature in scope that already has a `design/` directory. Each sub-agent:

- Reads `_meta/design.md` + its feature's `design/` only (never the codebase; never other features).
- Applies MADR consequences mechanically — edit-and-cite the existing ADRs, add an "Inherited meta-design constraints" block (pattern at `docs/biological-workflow/design/03-decisions.md § Inherited meta-design constraints`), never re-architects.
- Appends feature-local time-boxed deferrals to `_meta/deferred.md` using a **pre-assigned ID block** (20 IDs; the main agent allocates blocks in `/meta-apply` Phase 2 before parallel dispatch, which keeps the parallel writes race-free without any cross-agent coordination).
- Returns a diff summary OR a `MANUAL_REDESIGN_NEEDED` flag with a one-paragraph explanation of why the MADR can't be applied mechanically.

Main agent collects the per-feature diff summaries, surfaces a single compound gate (options: approve all / reject some / redo one / drill into a specific diff), and on approval dispatches `architect` subagents in parallel — one per approved feature — whose per-feature verdicts (`READY` / `NEEDS ITERATION` / `NEEDS DISCUSSION`) are collected and surfaced together. Frontmatter status flips `DRAFT → APPROVED` at the compound gate; `NEEDS` verdicts roll back to `DRAFT` automatically.

Cost shape for N = 5 features:

| Measure | Serial N × `/design` | `/meta-apply` |
|---|---|---|
| Tokens (dominant: redundant upstream re-reads) | ~300–400K | ~100–150K |
| Wall-clock (dominant: serial main-agent dispatch) | N × one-`/design` latency | ~max(per-feature) latency |
| User gate turns | ~2N (scope + bundle per feature) | 1 compound + 0 (architect verdicts are surfaced, not gated) |

### Serial (fallback) — `/design <feature>`

Use `/design` directly for:

- Features flagged `MANUAL_REDESIGN_NEEDED` by `/meta-apply` (sub-agent correctly refused to speculate; feature needs a real re-think, not an edit)
- Features with no `design/` yet (first-drafts are outside `feature-reviser`'s scope)
- Features where `architect` returned `NEEDS ITERATION` or `NEEDS DISCUSSION` after `/meta-apply` — loop with `/design <slug>` until verdict is `READY`
- Any feature where you want a hand-paced, fine-grained gate on the edits rather than a batch compound review

### Structural analogy

`/meta-apply` is to the propagation step what `/review --as all` is to the review step: narrow-scope sub-agents dispatched in parallel via a single message, with one compound gate. `meta-architect` owns authoring (`/meta-map`, `/meta-design`, `/meta-plan`); `feature-reviser` owns propagation; `architect` owns gating. Each agent's scope stays narrow; the command layer orchestrates parallelism.

### The MANUAL flag is a correctness signal

A `feature-reviser` that emits `MANUAL_REDESIGN_NEEDED` is behaving exactly right — it saw a MADR whose consequence couldn't be applied mechanically, and refused to write a speculative revision the user would then have to undo. Treat MANUAL as a routing decision ("this feature is in the `/design` lane, not the `/meta-apply` lane"), not a failure. Anti-pattern: interpreting MANUAL as a defect and re-prompting the sub-agent to try anyway.

---

## When to re-run meta

Five triggers:

1. **A new feature enters the wave** → re-run `/meta-map` and optionally `/meta-design`. The new feature's concerns may collapse into existing MADRs or surface fresh ones.
2. **A per-feature `/design` surfaces a MADR conflict** → re-run `/meta-design` with the user's resolution (either revising the MADR or relaxing it).
3. **Implementation reveals a wrong assumption in the meta-plan** → re-run `/meta-plan` with the new constraint (e.g., a P0 feature turns out to depend on a P1 feature).
4. **A new portfolio wave starts** → fresh `/meta-map`. Old `_meta/` artifacts move to `docs/_meta/.archive/{date}/` (manual; not automated).
5. **`/meta-apply` emitted `MANUAL_REDESIGN_NEEDED` for one or more features** → NOT a re-run of meta. Run `/design <slug>` per flagged feature. If the flag pattern is systemic (e.g., the same MADR keeps getting flagged across features), re-run `/meta-design` to revise *that specific MADR* — the portfolio decision may be too aggressive.

---

## Inheritance rules

- A per-feature `/design` reads `docs/_meta/design.md` (if present) during Phase 1.
- Every MADR in `_meta/design.md` is binding on the features named in its "Features affected" line.
- If the per-feature design needs to make a decision that contradicts a MADR, the `/design` command STOPS and forces the user to choose:
  1. Re-run `/meta-design` to revise the MADR
  2. Override locally and document the override under "Meta overrides" in the feature's `03-decisions.md`
  3. Cancel the feature decision

The override path exists deliberately — sometimes a feature has a reason to deviate, and forcing the meta-design to flex for one feature would make the MADR weaker for the others. But the override must be explicit and visible.

---

## Deferred items

`_meta/deferred.md` is the portfolio's roadmap of "we said we'd come back to this." Each row carries provenance (which doc, what was deferred, what blocks reactivation) so an architect can re-prioritise across the portfolio without re-reading every source archive.

Two commands append, sharing one schema (defined in `commands/design.md` Phase 3b):

**Per-feature `/design` Phase 3b** — when it writes its "Deliberate departures from synthesis/audit" section, it splits the items:
- **Permanent rejections** — stay in the feature's `03-decisions.md`. Not portfolio-relevant.
- **Time-boxed deferrals** — ALSO appended to `_meta/deferred.md` with `Source feature: <slug>`.

**`/meta-design` Phase 6.5** — after Gate 1 approval, it walks `_meta/design.md` and appends:
- Items in `Deliberate meta-level rejections § Time-boxed sub-items` (the meta-architect agent populates this subsection per its output template).
- `OQ-M-N` items whose resolution line says `still-open — deferred to ...` or `parked-in-deferred`.
- MADR Consequences referencing future rounds (e.g., "no further keys in Round 2").

Meta-level rows use `Source feature: _meta` and `Source doc: _meta/design.md § <section>`. Permanent meta-level rejections (no time horizon, no monitoring trigger) stay only in `_meta/design.md` — they are not roadmap items.

`/meta-plan` consumes `deferred.md` for its "Parking items" section regardless of which command appended each row.

---

## Worked example — pathway-explorer's 5 features

**Wave**: `normalisation`, `umap`, `color-encoding`, `theme-bundles`, `biological-workflow`. All have `map.md` and `synthesis.md`; some have `design/`.

### After `/meta-map`

`docs/_meta/map.md` would capture:

- **Features in scope**: 5 (status varies)
- **Shared touch-points**: `pathway_explorer/similarity.py` (touched by 4); `data_loader.py:151` (touched by 3); `main.py` (touched by 3); `html_generator.py` (touched by 4)
- **Convergent concerns**: 6 (the list at the top of this doc)
- **Inter-feature dependencies**: `theme-bundles → normalisation (within_method_z)`, `biological-workflow → umap (lasso primitive)`, etc.

### After `/meta-design`

`docs/_meta/design.md` would carry MADRs such as:

- **MADR-001**: Within-method z-score is the canonical commensurable score (resolves score incommensurability; affects `normalisation`, `theme-bundles`, `color-encoding`)
- **MADR-002**: Per-family colour scales are mandatory wherever multiple entity families share an axis (resolves shared-canvas implication; affects `color-encoding`, `theme-bundles`, `umap`)
- **MADR-003**: Leading-edge gene set is the source of truth for in-set membership; full-set is a separate, opt-in field (resolves LE/full-set ambiguity; affects `theme-bundles`, `biological-workflow`)
- **MADR-004**: Static-HTML payload budget is 5 MB compressed; >5 MB triggers progressive loading (resolves payload-blast-radius; affects `biological-workflow`, `theme-bundles`)
- **MADR-005**: Cluster-as-convergence claims require a permutation null model; absence forces softer language (resolves null-model gap; affects `theme-bundles`, `biological-workflow`)
- **MADR-006**: Schema rename `signed_sig → within_method_z` is portfolio-coordinated (consequence of MADR-001; affects `normalisation`, `theme-bundles`)

### After `/meta-plan`

`docs/_meta/plan.md` would sequence:

- **P0**: `normalisation phase-01` (introduce `within_method_z` column — unblocks 3 downstream features)
- **P1**: `umap phase-01` + `normalisation phase-02` in parallel (no shared files)
- **P2**: `theme-bundles phase-01` + `color-encoding phase-01` (both unblocked by P0 + P1)
- **P3**: `biological-workflow phase-01` (lasso + full-set gene support)
- **P4+**: Round-2 work (per-feature phase-02s as they're drafted)

The collision check would warn that P2 has both `theme-bundles phase-01` and `color-encoding phase-01` touching `html_generator.py` — likely on different functions, but flagged for human resolution.

### Without the meta-layer

Each per-feature `/design` would have made its own choice on:

- What "score" means (z-score? signed sig? raw NES?) — likely 5 different answers
- Whether full-set or LE is the primary gene reference — likely 2-3 different answers
- Whether colour scale is global or per-family — likely contradictions between `color-encoding` and `theme-bundles`
- Whether implementation order matters — likely no coordination, leading to KeyErrors at integration

The architect would resolve these in their head while reading 5 separate `synthesis.md` files. With 5 features × 6 reviewers, that's 30 documents to cross-reference. The meta-layer reduces this to 3 documents (`_meta/{map, design, plan}.md`).

---

## Cost

Per portfolio wave of N features:

- Per-feature: `N × (map + maybe-review + maybe-synthesis + design + plan)` invocations. Same as without meta-layer.
- Meta: `+ 1 (meta-map) + 1 (meta-design) + 0..1 (meta-plan)` invocations.

Overhead: 2-3 invocations per wave, regardless of N. For N = 5: ~3 meta vs ~15 per-feature → ~20% overhead. For N = 1: 0% overhead (don't invoke meta).

The amortised cost drops as N grows. The continuity cost (architect's mental refraction-inversion) is what you're saving — that scales O(N²) without the meta-layer.

---

## Limitations

What the meta-layer intentionally does NOT do:

- **Does not auto-resolve conflicts.** Every MADR is human-approved at Gate 1. Every collision in `_meta/plan.md` is flagged, never silently sequenced.
- **Does not re-run reviewers.** Specialists stay per-feature. There is no "meta-review."
- **Does not dispatch new architect reviews.** No meta-architect gate in v1 (see [03-decisions.md § ADR-016](./03-decisions.md)). The user is the gate.
- **Does not create features.** Only consolidates existing ones.
- **Does not modify per-feature artifacts.** All edits to feature design docs go through `/design {feature}` — the meta-layer can only add MADRs that the next `/design` invocation will inherit.

---

## See also

- [01-architecture.md](./01-architecture.md) — agent catalog now includes `meta-architect`
- [02-behavior.md](./02-behavior.md) — per-stage behaviour, including new meta-stages 3a/3b/3c
- [03-decisions.md](./03-decisions.md) — ADR-014/015/016 cover meta-layer design choices
- [04-extending.md](./04-extending.md) — how to add a new meta-command if needed
