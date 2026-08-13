---
doc: architectural-slicing-audit
date: 2026-05-21
status: PROPOSED
author: user + assistant
scope: SciAgent-toolkit/{commands,agents,skills/architecture-first-dev}/ + docs/workflows/architect/
motivation: A pathway-explorer audit ran the "thin architectural slice + synthesis" pattern by hand — six sequential single-concern slice agents, then a synthesis agent applying cohesion / connascence / decoupling lenses. The pattern proved generalizable. This proposal packages it as a first-class workflow inside the architect harness so any future audit, on any module, follows the same shape.
---

# Architectural slicing audit — toolkit extension

## 1. Shape of the extension — recommendation: command suite under `/audit`, NOT a new architect-mentor stage

The audit pattern has three primitives — namespace setup, per-slice tracing, cross-slice synthesis. They map cleanly to three commands:

```
/audit <name> [--lenses <list>] [--slices <list>]   — bootstrap: namespace + README + slice catalog
/audit-slice <name> <slug> [--depends-on <slice-ids>] — run one slice agent (sequential by default)
/audit-synthesize <name> [--overlays <list>]        — synthesis across all completed slices
```

Why a command suite and not the alternatives:

- **Single slash command** (`/audit ...` with internal phases) would replicate the bad pattern of `/design` Phase-3a/3b/4/5/6 — one command, six gates, no resumption. Audits have N slices (variable N), each potentially deferred or hand-edited between runs. A single command can't naturally express "run slice 4 now, slice 5 tomorrow, re-synthesize when slice 7 lands." The existing harness already separates `/map`, `/review`, `/synthesize`, `/design` for exactly this reason — each is a resumable, idempotent stage.
- **Single subagent with internal multi-step flow** would hide gates from the main agent. The user's primary value-add at gate-time is rejecting / re-mapping / re-scoping a slice before it consumes a context window. Burying the loop inside a subagent removes that lever.
- **New stage inside the architect-mentor pipeline** (e.g., a 7th stage between `/verify` and the next feature) was tempting but is wrong: audits are not per-feature, they are per-module / per-subsystem, and they typically run *without* an upstream `map.md`. Forcing them into the `docs/{feature}/` namespace overloads the slug. They are siblings to `/design` and `/meta-design`, not their continuation.

The suite mirrors the proven shape of `/map → /review → /synthesize` while replacing the parallel reviewer panel with a sequential slice panel. The new suite shares zero stages with the per-feature pipeline by design — co-tenancy is solely at the architect-mentor skill level (see §8).

## 2. Slice template

Every slice produces exactly one file at `docs/_meta/architecture-audit/<date>-<name>/slices/NN_<slug>.md` with the following structure. Required sections are unmarked; **[OPTIONAL]** is honest about what may be empty for a small slice.

```markdown
---
date: YYYY-MM-DD
audit: <name>
slice: NN_<slug>
slice-axis: <one-of-catalog OR free-form>
depends-on: [NN_<slug>, ...]    # other slices this one referenced
lenses-applied: [cohesion, connascence, decoupling, <overlays>]
files-touched: <count>
branching-points: <count>
---

# Slice NN — <human title>

## Concern (one paragraph)
What this slice traces and *why this axis was chosen* — not what's found yet.

## Trace path
Linear walk through the code that implements this concern, in execution / data-flow order.
Every step cites `file:line`. Mermaid `flowchart` encouraged.

## Branching points
| # | file:line | Branch condition | Branches | Reasonable? (Y/N/?) | Note |
|---|-----------|------------------|----------|---------------------|------|

A "branching point" is any place the trace can take more than one path: `if/elif`, dispatch
table, polymorphic call, feature flag, optional config key. Count this and put it in
frontmatter (`branching-points`). The branching-factor metric is `Σ(branches − 1)` over the
table — a slice with 4 binary branches has factor 4; one 5-way dispatch has factor 4. Use it
as a smell-check, not a hard threshold.

## Cohesion read
Apply the cohesion ladder (coincidental → logical → temporal → procedural → communicational →
sequential → functional). For each module the trace passes through, name where on the ladder
it sits and cite the evidence. Functional cohesion is the only one that doesn't need a follow-up
question.

## Connascence findings
Apply Page-Jones taxonomy. Static forms (Name, Type, Convention, Algorithm, Position) and
dynamic forms (Execution-order, Timing, Value, Identity) — each finding tagged with form,
strength (1=weakest=Name, 9=strongest=Identity), locality (same-function / same-class /
same-module / cross-module / cross-package), and degree (how many sites).

| Form | Locality | Degree | Sites | Strength | Note |

The rule of thumb stays Page-Jones's: stronger connascence is acceptable only at higher
locality. Cross-module Connascence-of-Position (strength 5) is a finding; same-function
Connascence-of-Algorithm (strength 6) usually isn't.

## Decoupling proposals (optional — empty if no actionable suggestion)
Each proposal: current shape → proposed shape → which connascence form it weakens. **[OPTIONAL]**

## Extension-readiness implications [OPTIONAL — only if `extension-readiness` is in lenses]
What this slice tells us about future extension along the audit's strategic axis.

## Open questions for synthesis
Things only the synthesis agent (with cross-slice visibility) can resolve.

## Files touched
Bullet list; one line per file. The frontmatter `files-touched` count is the length of this list.
```

**Required vs optional matrix.** Required: frontmatter, Concern, Trace path, Branching points (table may be empty but the section must exist), Cohesion read, Connascence findings, Open questions, Files touched. Optional: Decoupling proposals, Extension-readiness implications, any user-requested overlay sections.

**Why quantify branching.** A slice with zero branching is rare and worth flagging (often means the code is so linear that the slice axis was probably wrong). A slice with very high branching factor signals either (a) an actual decision-density hotspot worth a refactor proposal, or (b) the slice axis is too broad and should be split. Frontmatter counts make this skimmable across all slices without re-reading.

## 3. Synthesis template

`docs/_meta/architecture-audit/<date>-<name>/synthesis.md`:

```markdown
---
date: YYYY-MM-DD
audit: <name>
slices: [NN_<slug>, ...]
lenses-primary: [cohesion, connascence, decoupling]
lenses-overlay: [<as configured>]
---

# Synthesis: <audit name>

## Headline
1–2 sentences. The single most surprising / load-bearing finding across all slices.

## Cross-slice convergence
Findings flagged by ≥2 slices. Naming each slice in parens, same convention as `synth.md` uses
for reviewers.

## Cross-slice divergence
Where slices contradict each other (e.g., slice-3 calls a module functionally cohesive that
slice-5 calls communicational). State both, do not flatten.

## Connascence atlas
A single table aggregating every connascence finding across slices, sorted by (strength desc,
locality asc, degree desc). The first 5 rows are the priority targets.

## Cohesion atlas
Per-module cohesion verdicts aggregated across slices, with the worst-cohesion module ranked
first.

## Decoupling roadmap
Prioritised proposals:

### Must-fix (P0)
Findings where connascence strength exceeds locality budget. Each entry: finding (slice ref),
proposal, estimated blast radius (files touched), risk class.

### Should-fix (P1)
Same shape, lower urgency.

### Watch (P2)
Smells that aren't actionable yet but should be re-checked on the next audit.

## Strategic overlays
Per overlay lens (extension-readiness, MVP-vs-deep, performance, security, …), one subsection
each. Empty subsection if the overlay was configured but the slices surfaced nothing.

## Open audit questions
Cross-slice questions still unresolved after synthesis. These bound the next audit's slice set.

## Files most-touched
Hot files across slices (deduplicated, ranked by appearance count). Helps the reader judge
where future refactors will collide.
```

**Ordering rationale.** Headline first (the reader may stop after section 1). Convergence before divergence (convergence is decision-grade; divergence requires judgment). Atlases before roadmap (the roadmap *cites* atlas rows, so the reader needs the atlas in working memory). Strategic overlays last — they are the user's overlays, not the primary lens output, and a reader skimming for "what must change" finds the decoupling roadmap above them.

**Prioritization rule.** P0 means the connascence strength of a finding exceeds its locality budget under Page-Jones's principle. P1 means strength = budget (acceptable, watch for drift). P2 means smell-only. If every finding ranks P0, the synth has failed — same rule as `synth.md` for reviewer synthesis. Forced re-prioritization.

## 4. Namespace conventions

Default: `docs/_meta/architecture-audit/<YYYY-MM-DD>-<name>/`. Living inside `_meta/` matches the existing meta-layer namespace (`_meta/map.md`, `_meta/design.md`); audits are cross-feature analyses and belong there. The date prefix supports multi-audit histories naturally:

```
docs/_meta/architecture-audit/
├── README.md                                 # audit index (see below)
├── 2026-05-12-pathway-explorer/
│   ├── README.md                             # scope, slices, status table
│   ├── slices/
│   │   ├── 01_data-lineage.md
│   │   ├── 02_identity-badges.md
│   │   └── ...
│   └── synthesis.md
└── 2026-08-30-coresh-cache/
    └── ...
```

**Audit index** (`docs/_meta/architecture-audit/README.md`) is a one-table-per-row index of every audit run, with date, name, slice count, P0 count, and a one-line headline. `/audit <name>` maintains this file (append-on-bootstrap, update-on-synthesize). The index is the answer to "what audits has anyone run on this codebase?"

## 5. Slice catalog — ship a default, allow free-form

A named catalog of common slice axes ships with the toolkit at `docs/workflows/architect/audit-slices.md`. Proposed default catalog (the names are deliberately code-shape-driven, not domain-driven):

| Slice axis | Traces | Typical first finding |
|---|---|---|
| `data-lineage` | A datum's birth → mutation → death | Connascence-of-Type across boundary layers |
| `event-flow` | Event production → routing → consumption | Connascence-of-Execution-order |
| `state-machine` | All transitions of a stateful entity | Connascence-of-Algorithm replicated across handlers |
| `visual-channels` | (UI) data → visual encoding (color/size/shape) | Connascence-of-Meaning between data layer and view |
| `lifecycle` | An object's creation → use → teardown | Connascence-of-Timing |
| `persistence` | Read/write boundaries with storage | Connascence-of-Convention with on-disk format |
| `ipc-rpc` | Inter-process / network boundaries | Connascence-of-Type + serialisation drift |
| `extension-points` | Plugin / hook / strategy surfaces | Connascence-of-Name pollution |
| `filter-logic` | Predicate construction → application → display | Connascence-of-Position in argument lists |
| `entry-points` | Module's public API surface | Cohesion downgrades at the boundary |
| `error-paths` | Error production → propagation → presentation | Procedural cohesion masquerading as functional |
| `config-resolution` | Config source → merge → consumption | Connascence-of-Convention + Value |

`/audit-slice <name> <slug>` accepts a catalog name OR a free-form slug. Catalog names auto-fill the slice's Concern section with a starter template; free-form slugs require the user to write the Concern by hand. The catalog is a recommendation set, not a closed list — new entries are added by PR.

## 6. Sequencing — sequential default with explicit parallel override

Default: strictly sequential. Each `/audit-slice` invocation runs one slice agent, gates, and writes one file. The user runs them in order.

Override: `/audit-slice <name> --batch <slug1>,<slug2>,<slug3>` runs the listed slices in parallel as one dispatch (same idiom as `/review --as ...`).

**Independence detection** (user-overridable, not automatic). Two slices are *eligible* for parallel batching when:

1. Neither lists the other in its `depends-on` frontmatter.
2. Their planned file overlap (slice-axis × module) does not exceed 30% of either's predicted blast radius.

The harness will not auto-batch on this rule — too many false positives, and the user's reading rhythm matters. Instead, after each `/audit-slice` completes, the assistant surfaces a one-line "candidates for parallel batching" suggestion based on the rule, and the user opts in via `--batch`. The default conservative-sequential behavior matches the pathway-explorer audit's actual cadence.

## 7. Synthesis frame configuration

Lenses are configured at `/audit-synthesize` time (not at bootstrap), so the same set of slice docs can be re-synthesized under different lens overlays without re-slicing.

```
/audit-synthesize <name>                                   # defaults
/audit-synthesize <name> --overlays extension-readiness    # add one overlay
/audit-synthesize <name> --overlays performance,security   # add multiple
/audit-synthesize <name> --lenses cohesion,connascence     # narrow primaries
```

**Defaults.** Primary lenses = `cohesion, connascence, decoupling` (always present, cannot be narrowed below `connascence`). Optional overlays:

- `extension-readiness` — what each finding means for the module's ability to absorb new use-cases
- `mvp-vs-deep-vertical` — tension between minimal-viable surface and deep functional verticals
- `performance` — connascence findings tagged with hot-path implications
- `security` — connascence findings tagged with trust-boundary implications
- `testability` — cohesion findings tagged with mock surface area

Overlays are additive — each adds one subsection to the synthesis `## Strategic overlays` block. They never replace the primary lenses; they re-read the primary findings through an additional axis. This composition matches how `/review` reviewers stack on top of `map.md` without replacing it.

Each overlay is a small agent definition (`agents/audit-overlay-<name>.md`) following the same shape as the reviewer panel agents (e.g., `agents/divergent.md`). Adding a new overlay = one new agent file + one entry in the synth's prompt template.

## 8. Integration with architect-mentor

**Recommendation: SIBLING workflow, not a new stage.**

The architect-mentor pipeline (`map → review → synthesize → design → plan → implement → verify`) is *forward-looking* — it produces design artifacts for changes that will happen. The audit workflow is *backward-looking* — it characterises code that already exists, without a change in flight. Cramming audits into the per-feature pipeline would conflate two distinct epistemic activities. They share the architect-mentor skill's vocabulary (slices ≈ reviewers, synthesis ≈ synthesis, gates ≈ gates) but they answer different questions.

Concretely:

- The audit suite lives at `commands/audit*.md`, `agents/slice-tracer.md`, `agents/audit-synth.md`.
- The skill `architecture-first-dev/SKILL.md` gains a short "Audit workflow" section pointing at `docs/workflows/architect/audit/00-quickstart.md` (new file, parallel to the existing `00-quickstart.md`). The skill stays one document; the audit subtree is internally self-contained.
- The state-table in the skill gets two new rows: "user wants to characterise an existing module" → `/audit <name>`; "audit slices done, want the cross-cutting picture" → `/audit-synthesize <name>`.
- No changes to `/map`, `/review`, `/synthesize`, `/design`, `/plan`, `/implement`, `/verify`, or any meta-* command. Co-tenancy is purely additive.

The cleanest signal that this is sibling and not stage: an audit can output a P0 finding that *triggers* a new architect-mentor feature (`/map <fix-slug>` etc.). Audits feed into the per-feature pipeline, they do not extend it.

## 9. Deferred to a later iteration

- **Architect-gate equivalent for synthesis.** The per-feature pipeline runs `architect` as a verdict-issuing agent after `/design`. The audit synthesis has no such gate yet — slice quality is the user's read. Add a `audit-architect` agent issuing READY/NEEDS-ITERATION on the synthesis in v2.
- **Iterative slices.** No `/audit-slice --iterate` analogue to `/review --iterate` in v1. A slice is one-shot; if findings are weak, the user re-runs the same slice. v2 could add cross-informed iteration once slice quality stabilises in real use.
- **Cross-audit diffing.** Two audits of the same module a year apart should produce a diff: which findings closed, which regressed, which are new. The namespace supports it (`<date>-<name>/`); the tooling is v2.
- **Slice templates per language ecosystem.** A Python-heavy codebase's `event-flow` slice differs structurally from a TypeScript front-end's. v2 ships ecosystem-specific Concern templates inside the slice catalog.

---

## References

- Page-Jones, *What Every Programmer Should Know About Object-Oriented Design* (Dorset House, 1995) — connascence taxonomy (Name, Type, Convention, Algorithm, Position, Execution-order, Timing, Value, Identity), strength ordering, locality budget rule.
- Yourdon & Constantine, *Structured Design* (Yourdon Press, 1979) — cohesion ladder used in §2.
- `skills/architecture-first-dev/references/review.md` — parallel-dispatch idiom borrowed for `--batch`, two-regimes distinction borrowed for sequential-default/parallel-override.
- `commands/synthesize.md` and `agents/synth.md` — synthesis output structure (Headline / Convergent / Divergent / P0/P1/P2) borrowed wholesale.
- `docs/workflows/architect/00-quickstart.md` — navigation cadence and "tips after stage" pattern to be replicated in `audit/00-quickstart.md`.
- `docs/workflows/architect/proposals/2026-04-22-iterative-review-MADR.md` — `.history/` archive idiom for round-based artifacts; the audit suite reuses the same archive shape if `/audit-synthesize` is re-run after new slices land.
