# Plan — Craft as Architecture: making the "mind-palace" mindset a first-class artifact

**Date:** 2026-06-04 · **Author:** Opus 4.8 (investigation + synthesis) · **Status:** narrative report, for discussion (no changes made)
**Scope:** `SciAgent-toolkit` (primary) with one boundary note for `scbio-docker`
**Evidence base:** live audit of `DC_mouse_cancer` (a real analysis project at `@v3.1.0`), the toolkit at `v3.1.0-1-gf0f6c9d`, and the scbio-docker ↔ toolkit boundary docs

---

## 0. Why this document exists

The prompt was not "fix a bug." It was an observation worth taking seriously: *the
environment itself set the working tone.* An agent dropped into `DC_mouse_cancer`
behaved tidily — numbered stages, config-not-hardcoding, captioned artifacts,
decision logs — **not because it was told to in the moment, but because the repository
it woke up in already encoded that disposition.** The structure was the prompt.

That is a real and somewhat rare phenomenon, and it is worth naming before we refine
anything: **in an agent-operated repository, the layout, the docs, and the context
files are not documentation _of_ the work — they are the _interface_ through which the
work is done.** They are agent-facing artistry first, human-facing second. The
"mind-palace" framing is apt: the rooms, their labels, and the order you walk them
through are load-bearing cognition, not decoration.

This report answers the user's standing question —

> *should I have any general skills or AGENTS.md-managed blocks that define the
> mindset, the code style, the architectural preferences/characteristics?*

— and the answer, with evidence, is: **yes, but the mindset already exists and already
works; the problem is that it is implicit, triplicated, and drifting across three
uncoordinated homes.** The refinement is not "write more rules." It is *consolidate the
craft into one single-sourced, toolkit-owned, always-on artifact, and stop letting each
project's copy drift.* The sections below show where the mindset lives today, where it
leaks, and four concrete proposals.

---

## 1. Executive summary

1. **The mindset is real and it is working.** `DC_mouse_cancer` adheres to the scaffold
   with unusual fidelity — config-driven, compute/viz split, captioned artifacts,
   dated decision logs, cross-language consistency (R mirrors Python). This is the
   environment doing its job. (§2)

2. **But the craft philosophy lives in three places that do not agree.** The
   per-project `AGENTS.md` body states **4 "Critical rules"**; `docs/guidelines/`
   states a different **"Five Rules"**; and a third copy is *emergent* in the scaffold
   itself (`new.sh`, the templates). These are three statements of one philosophy, and
   they have already diverged. (§3)

3. **The named "single source of truth" is not actually wired in.** `docs/guidelines/`
   calls itself the SSOT, but the analysis `AGENTS.md` template **never references it** —
   and the guidelines are bulk-RNA-seq-shaped (`1.x.*.R`, limma-voom, DESeq2) while real
   projects run scRNA/scanpy/scVI with `NN_<slug>.py` naming. The SSOT describes a
   different pipeline than the one agents actually execute. (§3.3)

4. **The mindset is frozen per-project, so it cannot be improved centrally.** The craft
   rules live in the *hand-editable body* of each project's `AGENTS.md`, above the
   managed block. Fix a rule in the toolkit and existing projects never receive it —
   the opposite of how the `SCIAGENT:ROLES` managed block already works. (§3.4)

5. **There is no explicit aesthetic / craft statement anywhere.** The tidiness is
   inferred from structure and a handful of operational rules. The *disposition* — "we
   are building a place we will be proud to return to, not producing output for its own
   sake" — is nowhere written. It works by emergence, which means it is fragile to
   model changes, role changes, and template edits. (§3.5)

6. **A skill is the wrong vehicle for ethos; a managed block is the right one.**
   (Direct answer to the question.) Skills are load-on-demand, tool-scoped procedures;
   ethos must be always-on. Ethos belongs in the always-loaded `AGENTS.md`, single-sourced
   from the toolkit via a managed block — *not* in `skills/`. (§4.1)

7. **One boundary-doc drift surfaced as a side effect.** `scbio-docker/CLAUDE.md` still
   describes `setup-ai.sh` emitting `.mcp.json`, `GEMINI.md`, and `context.md`. The real
   `v3.1.0` harness emits a lean `CLAUDE.md → @AGENTS.md` shim with a managed block and
   *none* of those. Flag for reconciliation. (§5)

**The thesis in one line:** the environment already teaches the disposition — so make
the disposition a *maintained, single-sourced, toolkit-owned artifact* instead of an
emergent accident frozen into per-project copies.

---

## 2. What is working — the environment as a control surface

The evidence that the scaffold guides behavior is concrete, not aspirational.

| Signal | Evidence |
|--------|----------|
| **Config, not hardcoding** | Every script imports `config.py`; `analysis_config.yaml` (288 lines) is the sole home for thresholds, marker panels, palettes. No naked literals in callers. |
| **Compute/viz split** | Every stage is two files: `NN_<slug>.py` (no matplotlib import) + `NN_<slug>_viz.py` (reads checkpoint, plots). Enforced uniformly across ~10 scripts. |
| **Decision traceability** | `docs/_internal/reasoning/2026-06-03_object-build-plan.md` logs D1–D10 as *choice + rationale*, not just outcomes. Handoffs in `docs/_internal/sessions/` frame scope, results, and caveats. |
| **Artifact provenance** | Each `03_results/<stage>/README.md` captions every file as a one-sentence *finding* + a `Script | Function | Config | Input` table. |
| **Cross-language parity** | `04_convert_to_seurat.R` mirrors the Python style — section banners, config loading, defensive dimension assertions. The disposition crosses the language boundary. |
| **Commit hygiene** | `<scope>: <action>` subjects that link config edits to script edits ("config: add Lyz2 to tme_panel … TME dotplot only"; "viz: idempotent purge_figures per prefix"). |
| **Disciplined stopping** | Project delivered clustered objects and *stopped before annotation* — exactly as scoped — deferring cell-typing to the injected `mllmcelltype-consensus-annotation` skill. |

The lesson: **legibility is a control surface.** The agent did not need to be reminded to
be tidy in each turn; the rooms were already labelled. This is the asset we want to
protect and make robust — not re-derive per project.

---

## 3. Where the mindset leaks — three homes, no single source

The craft philosophy is stated, in whole or in part, in **three uncoordinated places**.
They overlap, they disagree in detail, and only one of them is actually loaded into the
agent's context.

### 3.1 Home A — the per-project `AGENTS.md` body (always loaded, hand-frozen)

`DC_mouse_cancer/AGENTS.md:22-29` states **four** Critical rules:

```
1. Config, not hardcoding.
2. Normalize, then visualize.
3. Read-only input data.
4. Cache expensive operations.
```

Plus inline conventions for figure standards (`:53-68`), artifact captions (`:70-76`),
the docs namespace + one-way reference rule (`:80-97`). This is the text the agent
*actually reads* (via `CLAUDE.md → @AGENTS.md`). It is good — and it is a hand-curated
**subset** of the canonical list in Home B, copied verbatim from the template at scaffold
time, never updated since.

### 3.2 Home B — `docs/guidelines/` (named SSOT, not loaded)

`docs/guidelines/core_architecture.md:10-26` states **five** rules:

```
1. Normalize Once, Visualize Many
2. Single Source of Truth
3. Checkpoint Everything Expensive
4. Master Tables as Bridges
5. Separation of Concerns
```

The guidelines README declares itself *"the single source of truth for analysis
patterns, coding conventions, and workflow standards"* and says *"`CLAUDE.md` (via
`@AGENTS.md`) references these guidelines rather than duplicating content."* **It does
not.** `DC_mouse_cancer/AGENTS.md` contains zero references to `docs/guidelines/`. The
SSOT is an island.

### 3.3 The two lists do not reconcile — and Home B is domain-mismatched

Mapping the two statements of "the same philosophy":

| Guidelines "Five Rules" (Home B) | Project "Critical rules" (Home A) | Status |
|----------------------------------|-----------------------------------|--------|
| Normalize Once, Visualize Many | Normalize, then visualize | ~same |
| Single Source of Truth | Config, not hardcoding | ~same, narrower wording |
| Checkpoint Everything Expensive | Cache expensive operations | ~same |
| Master Tables as Bridges | — | **dropped in A** |
| Separation of Concerns | — | **implicit in A (compute/viz)** |
| — | Read-only input data | **added in A, absent from B** |

Worse, Home B is written for **bulk RNA-seq**: its phase table is `1.x.*.R` / `2.x.*.R` /
`3.x.*.py`, its data flow is "limma-voom/DESeq2," its code examples are `run_gsea()` /
`DE_FDR_CUTOFF`. The real project is **scRNA**: scanpy/scVI, `NN_<slug>.py` two-digit
stages, AnnData checkpoints. The named SSOT describes a *different pipeline* than the one
agents run — so even if it were wired in, it would partly misdirect.

### 3.4 Home A is frozen per-project — improvements cannot propagate

The craft rules sit in the **hand-editable region** of `AGENTS.md`, *above* the managed
block (`<!-- BEGIN SCIAGENT:ROLES … -->` at `:105`). The managed block already proves the
right pattern: the toolkit *owns* it, stamps it, hashes it for drift detection, and can
re-render it on re-activation. But the *craft rules* enjoy none of that — they are a
copied template body. Improve rule #2 in the toolkit and every existing project keeps its
stale copy. The most important text in the file is the *least* maintainable.

### 3.5 No explicit disposition — only operational rules

Every rule above is *operational* ("put config here," "cache that"). Nowhere is the
*disposition* stated: why we work this way, what "done well" feels like, that the
repository is a place we return to and should be proud of. The tidiness is therefore
**emergent** — a happy consequence of good structure plus good operational rules — rather
than a stated value an agent can reason from when it hits a case the rules don't cover.
Emergent virtue is fragile: it survives only as long as the structure is untouched and the
model infers the rest. The user's instinct to *name* it is correct.

---

## 4. Proposals (for discussion — nothing implemented)

Four proposals, ordered by leverage. P1 is the heart of the answer.

### P1 — Promote the craft ethos to a toolkit-owned managed block

**Answers:** "should I have AGENTS.md-managed blocks that define the mindset?" → **yes.**

Introduce a second managed block, e.g. `<!-- BEGIN SCIAGENT:CRAFT vN hash=… -->`,
rendered by the toolkit alongside `SCIAGENT:ROLES`, carrying the **single canonical
statement** of disposition + the operational rules. Because it is managed (like the roles
block in `block.sh`), it is:

- **single-sourced** — defined once in the toolkit, stamped into every project;
- **updatable** — `sciagent activate`/a future `sciagent sync` re-renders it, so a rule
  improvement reaches all projects;
- **drift-detectable** — the existing SHA1 mechanism flags hand-edits.

Shape (terse — depth lives in guidelines, see P2):

```
<!-- BEGIN SCIAGENT:CRAFT v1 hash=… -->
## Craft

We are building a place we will return to — legible to the next agent and to a
human reading six months from now. Tidiness here is not decoration; it is the
interface through which the work is done. When a rule below does not cover a
case, choose the option a careful colleague would be glad to inherit.

1. Config, not hardcoding.
2. Normalize once; visualize many (compute never plots; viz never computes).
3. Checkpoint anything > ~1 min.
4. Read-only input; outputs are captioned with provenance.
5. Decisions are logged with rationale, not just outcomes.

Depth & rationale: docs/guidelines/ (single source of truth).
<!-- END SCIAGENT:CRAFT v1 hash=… -->
```

**Budget caveat:** `AGENTS.md` has a "never-grow-this-file (~150 lines)" rule (`:96`) and
the roles block alone is ~60 lines. Keep CRAFT ≤ ~20 lines and push detail to guidelines.

### P2 — Reconcile to one canonical rule-list and actually wire the SSOT

- Merge Home A's 4 and Home B's 5 into **one** canonical list (the union above is 6).
  The managed CRAFT block carries the terse list; `docs/guidelines/` carries the rationale.
- Make the wiring real: the CRAFT block (and the `AGENTS.md` template) must *reference*
  `docs/guidelines/` so the "SSOT" claim becomes true, not aspirational.
- **De-bulk-ify the guidelines.** `core_architecture.md` and `code_style.md` are
  bulk-RNA-seq-shaped while analysis projects are scRNA. Either (a) generalize the phase
  model to be assay-agnostic (`NN_<slug>` two-digit stages, AnnData/DGEList both first
  class), or (b) split guidelines per assay and let the role select which to reference.
  Today the SSOT documents a pipeline nobody runs.

### P3 — A "house style" the `code-reviewer` agent can cite, not a skill

**Answers:** "should I have *general skills* that define code style?" → **no, not skills.**

Skills are load-on-demand, tool-scoped procedures (`anndata`, `scanpy`, `scvi-basic`).
Ethos and house-style must be *always-on*; loading them on demand defeats the purpose, and
a skill the model has to *choose* to consult is a skill it will skip under pressure.
Instead:

- Keep **disposition** in the always-loaded CRAFT block (P1).
- Keep **house-style depth** in `docs/guidelines/code_style.md` (already exists, 526 lines).
- Give the existing `code-reviewer` and `doc-curator` agents an explicit instruction to
  *review against* `docs/guidelines/code_style.md`. That turns the style doc from passive
  prose into an enforced gate, without inventing a new skill or a new always-loaded file.

### P4 — Make the disposition legible to *humans* too (one public artifact)

The internal craft is invisible to a collaborator opening the public `docs/`. A single
short `docs/README.md`-level statement ("how this project is organized and why") — derived
from the same canonical source — closes the loop the user described: the environment that
disciplines the agent should also welcome the human. Low effort, high signalling value,
and it reinforces that these artifacts are *interfaces*, not paperwork.

---

## 5. Boundary note (scbio-docker) — stale harness description

Surfaced incidentally and worth a separate small fix (lives in scbio-docker, not here):

- `scbio-docker/CLAUDE.md` (Project Scaffolding + AI Tools sections) and
  `docs/repo-structure.md` describe `setup-ai.sh` creating `.mcp.json`, `GEMINI.md`, and
  `context.md`, with a richer MCP server roster.
- The real `v3.1.0` harness in `DC_mouse_cancer` is **lean**: `CLAUDE.md` is a one-line
  `@AGENTS.md` shim; `AGENTS.md` carries the managed `SCIAGENT:ROLES` block; there is **no**
  `.mcp.json`, `GEMINI.md`, or `context.md`. Scientific framing moved to
  `docs/_internal/scientific-context.md`.

Recommendation: reconcile the scbio-docker boundary docs to the lean reality (or, if MCP
is still intended, document it as optional). This is a doc-truth fix, not an architecture
change — but it matters precisely *because* these docs are the interface a new agent reads
first.

---

## 6. Open questions for our discussion

1. **One block or two?** Fold disposition into the existing `SCIAGENT:ROLES` block, or a
   separate `SCIAGENT:CRAFT` block? (Separate is cleaner to version and reason about;
   costs ~15 lines of the AGENTS.md budget.)
2. **Propagation verb.** Do we want a `sciagent sync` that re-stamps managed blocks into an
   existing project, so craft improvements reach past projects — or is re-`activate`
   enough? (Today there is no clean "update the frozen body" path.)
3. **Guidelines: generalize vs split.** Make one assay-agnostic guideline set, or per-assay
   sets selected by role? Bulk vs scRNA vs scATAC genuinely differ in phase shape.
4. **How much voice in a managed block?** The CRAFT preamble is where "mind-palace"
   disposition would actually live in agent context. How much *tone* do we want there
   vs. terse rules? (This is the one place tone is load-bearing, not ornamental.)
5. **Enforcement appetite.** Is wiring `code-reviewer`/`doc-curator` to gate against
   `code_style.md` (P3) the right amount of enforcement, or do we want a pre-commit/CI
   check on caption completeness and config-not-hardcoding (raised by the project audit)?

---

## 7. One-paragraph recommendation

The environment is already a teacher — that is the finding, and it is a good one. The
work now is to stop letting the lesson live in three drifting copies and one emergent
accident. Name the disposition, state it once, let the toolkit own it as a managed block
(P1), reconcile the rule-lists and wire the SSOT it claims to be (P2), enforce house-style
through the agents we already have rather than a new skill (P3), and let the human see what
the agent sees (P4). Skills stay what they are — load-on-demand procedures. The mindset
becomes what it should always have been: a single, maintained, always-on artifact — the
first room you enter in the palace, the one that tells you how to walk through the rest.
