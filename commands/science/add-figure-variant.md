# Add Figure Variant

Add a new figure family with its own grep-isolable namespace, strict compute→viz split, dual print+screen variants, and a gated PDF-open review. Opus designs the mini-plan and runs the review gate; Sonnet implements compute and viz; the acceptance gate runs `figure-audit` and verifies namespace isolation before the mandatory `captions` cleanup pass.

## Phase 0: Parse arguments

| Param | Shape | Default | Meaning |
|---|---|---|---|
| `slug` | `$ARGUMENTS[0]` (positional, **required**) | — | Figure family slug (e.g. `genotype-umap`). Used as the base name for scripts, artifacts, and the research dir. |
| `stage-id` | `$ARGUMENTS[1]` (positional, **required**) | — | Target stage identifier (e.g. `03_clustering`). Resolves the output root to `03_results/<stage-id>/`. |
| `--namespace-token` | flag | derived from `slug` by replacing `-` with `_` (e.g. `genotype_umap`) | A grep-isolable token that MUST prefix every new identifier, filename, and variable in the compute and viz scripts. Used in the namespace-isolation gate in Phase 5. |
| `--n-implementers` | flag | `1` | Number of parallel Sonnet implementers. `2` = dispatch compute and viz concurrently once both briefs are ready; `1` = sequential (default, safer for inter-script dependency). |

Flag parsing is order-independent. If either positional argument is missing, reject:

```
Usage: /add-figure-variant <slug> <stage-id> [--namespace-token <token>] [--n-implementers <int>]
```

Resolve and announce the run config:

```
/add-figure-variant {slug}  stage={stage-id}
  Research dir:      docs/_internal/research/{date}-{slug}/
  Results root:      03_results/{stage-id}/
  Namespace token:   {namespace-token}
  Implementers:      {n-implementers} (Sonnet)
  Review tier:       Opus
```

**STOP condition — unknown stage.** If `03_results/<stage-id>/` does not exist AND `<stage-id>` does not appear as a declared stage id in `02_analysis/config/analysis_config.yaml:stages`, do NOT create it silently. Stop:

```
Stage {stage-id} is not declared in analysis_config.yaml:stages and 03_results/{stage-id}/ does not exist.
Declare the stage first (add it to analysis_config.yaml:stages) or pass an existing stage-id.
```

Stop.

## Phase 1: Evidence → research (persist, not chat-only)

Create `docs/_internal/research/{date}-{slug}/` if it does not exist. Do NOT proceed if the directory cannot be created.

Capture the rationale and data contract for the new figure family into:

```
docs/_internal/research/{date}-{slug}/01_rationale.md
```

This file MUST exist on disk before Phase 2 begins. The file records:

```markdown
# Research: {slug}

**Date:** {date}  ·  **Stage:** {stage-id}  ·  **Namespace token:** {namespace-token}

## Scientific rationale
Why this figure family is needed; what biological question it addresses; what existing figures do not show.

## Data contract
- **Compute inputs:** what data the compute script reads (file paths, object names, columns used).
- **Checkpoint outputs:** what the compute script writes (`.rds`/`.h5ad` checkpoint path + schema).
- **Viz inputs:** the checkpoint(s) the viz script reads.
- **Figure outputs:** the artifact stems and sub-layout (`_overview/` or `by_contrast/<c>/`).

## Genotype / facet structure
List the facet levels or grouping variables the figure family iterates over (e.g. genotypes, contrasts, time-points). Empty if the figure is a single-panel overview.

## Acceptance criteria (human-readable)
What the figure must show, at what statistical threshold, for the family to be considered correct.
```

**STOP condition — rationale absent.** If `01_rationale.md` is not on disk after this phase, stop:

```
docs/_internal/research/{date}-{slug}/01_rationale.md was not written.
A chat-only rationale cannot anchor the mini-plan. Re-run or write the rationale manually.
```

Stop.

## Phase 2: Mini-plan (model tier: **Opus**)

Dispatch **one Opus planner**. The planner reads:

1. `docs/_internal/research/{date}-{slug}/01_rationale.md` — the data contract and acceptance criteria.
2. `02_analysis/config/analysis_config.yaml` — stage ids, `figures:` block (geometry, font floors, sub-layout names).
3. `skills/figure-style/SKILL.md` + `lib/figure-style/figure_helpers.{R,py}` — the figure-style contract (helper functions, anti-patterns, dual-variant semantics).
4. Existing `02_analysis/scripts/` file listing — to determine the next available `NN` script index.

The planner writes `docs/_internal/research/{date}-{slug}/02_miniplan.md`:

```markdown
# Mini-plan: {slug}

**Date:** {date}  ·  **Planner:** Opus  ·  **Namespace token:** {namespace-token}

## Namespace declaration
- **Token:** `{namespace-token}` — must prefix every new identifier (variables, function names,
  checkpoint filenames, figure stems). MUST NOT match any token already used in sibling scripts
  (verified by grep in Phase 5).
- **Script names:** `02_analysis/scripts/NN_{slug}_compute.{R,py}` and `NN_{slug}_viz.{R,py}`
  where NN is the next available index.
- **Checkpoint path:** `03_results/objects/{namespace-token}_checkpoint.rds` (or `.h5ad`).
- **Figure stems:** `{namespace-token}_<panel>` under `03_results/{stage-id}/figures/<sub-layout>/`.

## COMPUTE outputs
- What the compute script writes (checkpoint/table path, schema, row count expected).
- Must NOT contain any ggplot/matplotlib call. Must NOT write to `03_results/figures/`.

## VIZ outputs
- Figure stems and sub-layout (`_overview/` or `by_contrast/<c>/`).
- Dual variants per stem: `<stem>.print.pdf` (vector PDF, column geometry) + `<stem>.screen.png`
  (raster PNG, screen geometry) — both produced by `save_figure(variant="both")` or `save_overview()`.
- Same-stem source table: `tables/<sub-layout>/<stem>.csv`.
- README caption: path-qualified `## figures/<sub-layout>/<stem>.screen.png` section with
  non-empty `**How to read:**` block — written atomically by `save_overview()`.

## Acceptance gate
State verbatim: what artifacts must exist, what panels must be non-empty, what the
namespace-isolation grep must confirm, and what figure-audit criteria the family must pass.
```

**STOP condition — mini-plan absent.** If `02_miniplan.md` is not on disk after the Opus planner returns, stop:

```
docs/_internal/research/{date}-{slug}/02_miniplan.md was not written.
The compute and viz implementers cannot proceed without a persisted mini-plan.
```

Stop. Do NOT pass the plan inline to Phase 3.

## Phase 3: Compute first (model tier: **Sonnet**)

Dispatch **one Sonnet implementer** with the mini-plan as context. The implementer creates:

```
02_analysis/scripts/NN_{slug}_compute.{R,py}
```

where `NN` is the index declared in `02_miniplan.md`.

**Compute discipline (non-negotiable):**

- The compute script reads input data and writes a checkpoint / table to the path declared in the mini-plan.
- **NEVER calls `ggplot()`, `plt.subplots()`, `ggsave()`, `plt.savefig()`, `save_figure()`, or `save_overview()`.** Compute scripts do not plot. Any plotting call in a compute script is a hard violation of the figure-style contract.
- Every new object, variable, and function name in the compute script MUST begin with the `{namespace-token}` prefix (e.g., `genotype_umap_cells`, `genotype_umap_compute()`). No unprefixed new identifiers.
- The script must be runnable stand-alone from the project root: `Rscript 02_analysis/scripts/NN_{slug}_compute.R` or `python 02_analysis/scripts/NN_{slug}_compute.py`.

**After the compute script is written — VERIFY the checkpoint on disk before proceeding to viz:**

```
ls -lh 03_results/objects/{namespace-token}_checkpoint.{rds,h5ad}
```

**STOP condition — checkpoint missing.** If the checkpoint file is absent or empty after the compute script runs, STOP:

```
Compute checkpoint 03_results/objects/{namespace-token}_checkpoint.{rds,h5ad} is absent or empty.
Viz cannot proceed without a verified checkpoint. Fix the compute script and re-run it first.
```

Stop. Do NOT start the viz implementer.

## Phase 4: Viz (model tier: **Sonnet**)

Only after the checkpoint is verified on disk, dispatch **one Sonnet implementer** with the mini-plan and the confirmed checkpoint path. The implementer creates:

```
02_analysis/scripts/NN_{slug}_viz.{R,py}
```

**Viz discipline (non-negotiable):**

- Import the per-project figure-style shim:
  - R: `source("02_analysis/helpers/figure_style.R")` — exposes `FIG_CFG`, `project_theme()`, `save_figure()`, `save_overview()`, `contrast_path()`, `overview_path()`.
  - Python: `from helpers.figure_style import set_paper_style, save_overview, FIG_CFG` then `set_paper_style(config=FIG_CFG)`.
- Call `save_overview()` or `save_figure(variant="both")` for every figure stem — **never** `ggsave()` / `plt.savefig()` directly. The `variant="both"` parameter produces both variants from one call:
  - `<stem>.print.pdf` — vector PDF, column geometry (`width_column` × `height_column`), print font tier (`base_size_column`). Illustrator-editable text (`pdf.fonttype=42`).
  - `<stem>.screen.png` — raster PNG at `dpi`, screen geometry (`width` × `height`), screen font tier (`base_size`).
- `save_overview()` writes three things atomically: the dual-variant figure files, the same-stem source table (`tables/<sub-layout>/<stem>.csv`), and the path-qualified README caption with `**How to read:**`. This is the only sanctioned path for a figure + table + caption.
- **No inline `theme()` / `element_text(size=<num>)` / `ggsave(width=<literal>)` / raw hex color strings.** All style decisions go through `project_theme(config=FIG_CFG)` (R) or `set_paper_style(config=FIG_CFG)` (Python). Inline overrides break the dual-variant font floors silently.
- Every new identifier MUST begin with `{namespace-token}`.
- Cap categorical axes to `FIG_CFG$figures$top_n` (R) / `FIG_CFG["figures"]["top_n"]` (Python) before plotting.
- Use `direction_cue(value)` for signed labeling — never a bare `*` or raw colored dot.

After the viz script runs, verify outputs:

```
ls -lh 03_results/{stage-id}/figures/<sub-layout>/{namespace-token}_*.{pdf,png}
ls -lh 03_results/{stage-id}/tables/<sub-layout>/{namespace-token}_*.csv
grep -n "## figures/" 03_results/{stage-id}/README.md | grep {namespace-token}
```

**STOP condition — viz artifacts absent.** If any of the following are missing or empty after the viz script runs, STOP:

```
One or more viz artifacts are absent or empty:
  - <stem>.print.pdf missing            → dual-variant contract violated
  - <stem>.screen.png missing           → dual-variant contract violated
  - <stem>.csv missing                  → source-table adjacency violated
  - README.md caption section missing   → README-adjacency violated
Fix the viz script. Do not proceed to review.
```

Stop.

## Phase 5: Review — the gate (model tier: **Opus**)

Dispatch **one Opus reviewer**. The reviewer MUST perform all four checks. A review that lives only in chat is a failure; the reviewer MUST write its verdict to:

```
docs/_internal/reasoning/{date}_{slug}_review.md
```

### (a) Open the produced PDFs — confirm non-empty panels

Open each `<stem>.print.pdf` under `03_results/{stage-id}/figures/`. Confirm that:

- Every panel contains visible data (not blank axes).
- Text is present and not overlapping with panel borders.
- The print variant uses column geometry.

If any PDF is blank or unreadable, FAIL the review and stop.

### (b) Run `figure-audit` on the stage

Invoke the `figure-audit` subagent on `03_results/{stage-id}/`:

```
figure-audit  stage={stage-id}
```

`figure-audit` runs the D-theme legibility checklist (all ten criteria: dual-scale legibility, font-size floors, no truncated labels, top-N capping, line/point weight, unambiguous glyphs, row/column clustering, residualized channel, both variants present, path-qualified caption with How-to-read). The reviewer reads the verdict table and confirms:

- All new `{namespace-token}_*` figures have a `figure-audit` verdict of **PASS** on all ten criteria.
- Any **FAIL** is surfaced with the criterion letter and the exact fix required before the family is accepted.

### (c) Namespace isolation — grep gate

Run verbatim:

```bash
grep -rn "{namespace-token}" 02_analysis/scripts/ 03_results/{stage-id}/
```

**The namespace-isolation gate passes if and only if every match falls within:**

- `02_analysis/scripts/NN_{slug}_compute.{R,py}` — the compute script
- `02_analysis/scripts/NN_{slug}_viz.{R,py}` — the viz script
- `03_results/{stage-id}/figures/<sub-layout>/{namespace-token}_*` — figure artifacts
- `03_results/{stage-id}/tables/<sub-layout>/{namespace-token}_*` — source tables
- `03_results/{stage-id}/README.md` — captions section

**If the token appears in any other file** (a sibling script, a master table column not declared in the mini-plan, a different stage's README, an unrelated config), the review FAILS with:

```
Namespace isolation FAILED. Token "{namespace-token}" found outside intended scope:
  {file}:{line}: {matching line}
The figure family is leaking into sibling scripts/artifacts. Fix before accepting.
```

This gate is non-negotiable — it is the primary mechanism that makes the figure family isolatable, auditable, and safely removable.

### (d) Caption completeness

For every new figure stem, confirm `03_results/{stage-id}/README.md` contains:

- A `## figures/<sub-layout>/<stem>.screen.png` heading (path-qualified, exact path).
- A non-empty `**How to read:**` block covering glyph semantics, sign convention, and claim tier.

A figure with a missing or empty `**How to read:**` is an **automatic FAIL** on this gate.

### Review verdict

The reviewer writes the verdict to `docs/_internal/reasoning/{date}_{slug}_review.md`:

```markdown
# Review verdict: {slug}

**Date:** {date}  ·  **Reviewer:** Opus  ·  **Stage:** {stage-id}

## (a) PDF panels   — {PASS|FAIL}
{evidence: file paths opened, panel descriptions or failure details}

## (b) figure-audit — {PASS|FAIL}
{paste the figure-audit verdict table for {namespace-token}_* figures}

## (c) Namespace isolation — {PASS|FAIL}
Command run:
  grep -rn "{namespace-token}" 02_analysis/scripts/ 03_results/{stage-id}/
Matches found: {count}
Files in scope: {list}
Leaks outside scope: {none | list of offending file:line}

## (d) Caption completeness — {PASS|FAIL}
{list of stems checked; caption path; How-to-read present yes/no}

## Overall: {PASS|FAIL}
{1–2 sentence summary}
```

**STOP condition — review fails.** If any of the four checks is **FAIL**, STOP:

```
Review FAILED on check ({a|b|c|d}).
See docs/_internal/reasoning/{date}_{slug}_review.md for details.
Fix the indicated issue and re-run the review phase.
```

Do NOT proceed to Phase 6 until all four checks are **PASS**.

## Phase 6: Mandatory `captions` cleanup pass

After the review passes, dispatch the **`captions`** agent on `03_results/{stage-id}/` as a mandatory cleanup pass (README-adjacency rule):

```
captions  stage={stage-id}
```

The `captions` agent sweeps every file under `03_results/{stage-id}/figures/` and verifies the sibling README contains a path-qualified caption section with a non-empty `**How to read:**` block. Any gap is filled by `captions`. This pass is mandatory even when `save_overview()` was used correctly — `save_overview()` writes the caption at figure-creation time, but a manual edit or a re-run may have left the README in an inconsistent state.

After `captions` returns, verify:

```bash
grep -c "## figures/" 03_results/{stage-id}/README.md
```

The count must equal the number of figure stems produced. If any caption is still missing, re-run `captions` with the specific stem as context.

## Output block

When all phases and checks are green:

```
/add-figure-variant {slug}  stage={stage-id}  complete.

Research dir:    docs/_internal/research/{date}-{slug}/
  01_rationale.md   — data contract + scientific rationale
  02_miniplan.md    — Opus mini-plan (namespace decl, compute/viz outputs, acceptance gate)

Scripts:
  02_analysis/scripts/NN_{slug}_compute.{R,py}   — COMPUTE only (no plots)
  02_analysis/scripts/NN_{slug}_viz.{R,py}        — VIZ only (reads checkpoint, calls save_overview)

Artifacts produced:
  03_results/{stage-id}/figures/<sub-layout>/{namespace-token}_*.print.pdf  (print variant)
  03_results/{stage-id}/figures/<sub-layout>/{namespace-token}_*.screen.png (screen variant)
  03_results/{stage-id}/tables/<sub-layout>/{namespace-token}_*.csv         (source tables)
  03_results/{stage-id}/README.md                                            (captions, updated)

Review:
  docs/_internal/reasoning/{date}_{slug}_review.md   (Opus verdict — all 4 checks PASS)

Namespace isolation: grep "{namespace-token}" — confined to intended scripts/artifacts ✓
Dual variants:       both .print.pdf + .screen.png present for every stem ✓
Caption completeness: all stems have path-qualified README section + How-to-read ✓
```

If stopped at a failed phase, report the failing phase, the missing or offending artifact path, and the recommended fix loop.

## Rules

1. **Compute never plots; viz never computes.** The compute script writes only checkpoints/tables. The viz script reads only the checkpoint and calls `save_overview()`/`save_figure(variant="both")`. Any crossover is a hard violation.
2. **Checkpoint verified before viz starts.** The `ls -lh` checkpoint check after Phase 3 is mandatory. A viz that starts without a verified checkpoint is not reproducible.
3. **Dual variants from one call.** Every figure stem must have both `<stem>.print.pdf` and `<stem>.screen.png`. These are produced by `save_figure(variant="both")` or `save_overview()` — never by two separate save calls with different parameters.
4. **Namespace token prefixes everything new.** Every new identifier in both scripts, every checkpoint filename, every figure/table stem uses `{namespace-token}` as a prefix. The grep gate in Phase 5(c) is the enforcement mechanism.
5. **Namespace isolation is non-negotiable.** If the token appears in sibling scripts or unrelated artifacts, the review fails. A leaking namespace corrupts the audit trail and makes the family non-removable.
6. **`save_overview()` is the only sanctioned figure+table+caption path.** It writes three things atomically. `save_figure()` alone is permitted only when the source table and caption are written separately in the same script call. Never ship a figure without its same-stem CSV and its `## figures/...` README section.
7. **No inline style overrides.** No `ggsave(width=<literal>)`, `element_text(size=<num>)`, `theme(...)` blocks, or raw hex color strings in the viz script. All style decisions go through `project_theme(config=FIG_CFG)` / `set_paper_style(config=FIG_CFG)`. Inline overrides silently break the font floors on the print variant.
8. **Persist every decision before proceeding.** Research rationale → `01_rationale.md`. Mini-plan → `02_miniplan.md`. Review verdict → `{date}_{slug}_review.md`. A decision with no trace is non-reproducible.
9. **Model tiering is explicit.** Mini-planner = **Opus**. Compute implementer = **Sonnet**. Viz implementer = **Sonnet**. Reviewer = **Opus**. `captions` cleanup = `captions` agent (Sonnet). State the tier at every dispatch.
10. **Mandatory `captions` pass always runs.** Even when `save_overview()` wrote captions at creation time, run the `captions` agent as the final cleanup pass — it is the backstop for README-adjacency.
