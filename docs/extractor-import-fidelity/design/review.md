---
date: 2026-05-26
feature: extractor-import-fidelity
kind: review
gate: pre-/plan
reviewer: architect
depends_on: [decisions.md, map.md]
repo_root: /workspaces/DC_hum_verse/01_modules/SciAgent-toolkit
---

# Architecture Review: extractor-import-fidelity

**Reviewed:** `docs/extractor-import-fidelity/design/decisions.md`, `docs/extractor-import-fidelity/map.md`
**Verified against live source:** `extract_components.py`, `validate_components.py`, and the
witness snapshot `pathway-explorer/docs/_meta/architecture-audit/2026-05-25-postrefactor/components.json`.

## Verdict

**PASS-WITH-NITS**

The load-bearing decision (D-1 resolver) is correct by construction for every
relative-import form, including the multi-dot/bare forms with no live witness;
the regular-module control path is provably byte-identical. The gate hardening
(D-3) is sound and would have correctly flagged the masked `run-metadata`
defect without false-positiving on the existing snapshot. D-4 correctly
identifies exactly one redundant audit-asserted edge and the phase ordering is
sound. No BLOCK-level issue: the golden oracle is NOT circular as specified
(the expected edge set is hand-authored from import semantics in the design
itself; the extractor run only confirms), and no resolver form is mishandled.
The nits are (N-1) a weakly-phrased oracle-finalization instruction that should
be tightened to forbid authoring-from-the-code-under-test, and (N-2) one
Phase-C sequencing dependency (channel-owners) the design should call out as a
verification gate.

---

## D-1 — Resolver correctness  [PASS]

**Anchors verified against source.** Every `file:line` the design pins is
correct in the live file:

- `module_name_for_file` strips trailing `__init__` at
  `extract_components.py:468-469` (verified).
- `collect_imports` signature at `:493`; `package_parts = source_module.split(".")[:-1] if source_module else []`
  at `:510` (verified — design Edit 2 targets the exact line).
- `resolve_import_from_base` level-walk
  `base_parts = package_parts[: len(package_parts) - (node.level - 1)]` at `:533`
  (verified verbatim).
- Sole caller `imports = collect_imports(absolute, source_module)` at `:585`,
  with `relative` in scope from the loop at `:582` (verified). Grep confirms
  **exactly one** caller — `extract_components.py:585`, def at `:493`, no others.
  Design Edit 3 is accurate.
- Self-import guard `if target_file == relative: continue` at `:595-596`
  (verified — the design's "no spurious self-edges" claim holds).

**Helper correctness — every form walked by hand.** Let `P` = containing
package parts. The helper returns `parts` for `*/__init__.py`, else `parts[:-1]`.

| Form | Trace | Result |
|------|-------|--------|
| regular sibling `beta.py: from .alpha` | `parts=['pkg']`, not init → `['pkg']`; level1 → base `pkg.alpha` | ✓ resolves `pkg/alpha.py` |
| init sibling `api/__init__.py: from .channels` | `parts=['pathway_explorer','api']`, init → keep all; level1 → base `…api.channels` | ✓ (was collapsing to root) |
| **multi-dot** `sub/__init__.py: from ..alpha` | P=`['pkg','sub']`; level2 → `P[:len(P)-1]=['pkg']`; +alpha → `pkg.alpha` | ✓ resolves up to `pkg/alpha.py` |
| **bare** `pkg/__init__.py: from . import beta` | P=`['pkg']`; level1, module None → base `pkg`; candidates `pkg.beta`,`pkg` → `pkg.beta` | ✓ resolves `pkg/beta.py` |
| absolute `from pkg.mod import Y` | level0 → `node.module` untouched (`:530-531`) | ✓ never touches package_parts |

The design's claim that **the L533 arithmetic is already correct for `level>=1`
once fed a correct `P`** is true: for `level==1` it keeps all of `P`; each
additional level drops one ancestor (`P[:-(level-1)]`). The `level==3` table row
(`P[:-2]`) equals `package_parts[:len-2]` for `len(P)>=2`. No change to
`:530-536` is needed — confirmed.

**(a) `/__init__.py` suffix robustness — interrogated the edge cases the gate
demanded:**

- **0-dot top-level regular module** (`foo.py` at root): `source_module='foo'`,
  `parts=['foo']`, not init → `['foo'][:-1]=[]`. A relative import in a
  top-level module is a Python error anyway; `[]` reproduces prior behaviour.
  No regression.
- **Top-level `__init__.py`**: `module_name_for_file` strips `__init__` →
  `best_module` becomes `""` and `return best_module or None` yields **None**
  (`:480`). The helper's `if not source_module: return []` guard fires FIRST, so
  the explicit `relative_path == "__init__.py"` branch is **dead code** — never
  reached because `source_module` is None. Harmless, but the design's helper
  carries an unreachable clause. Cosmetic only.
- **Namespace packages (PEP 420, no `__init__.py`)**: no init file exists to
  misclassify; constituent modules are regular `.py` handled by `parts[:-1]`,
  unchanged from today. Not a regression. The design folds this under
  "layout-agnostic by construction" rather than naming it — acceptable.
- **Module name derived from a non-matching root**: `module_name_for_file`
  returns None → helper returns `[]`. Same as today.

No relative-import form is mishandled by the helper. **Not BLOCK-worthy.**

**(c) Regular-module path provably identical.** Old: `source_module.split(".")[:-1]`.
New: `parts = source_module.split(".")`, then `parts[:-1]`. Byte-identical
expression. The `html_generator.py` 8/8 control is preserved by construction —
the helper only diverges on the `*/__init__.py` branch. ✓

**Implementation contract concreteness:** sufficient. An implementer holding
this doc can apply Edits 1-3 without re-deriving (signature change, exact line
replacement, single caller update, grep-to-verify already done — confirmed 1
caller). The invariant is stated in the docstring.

---

## D-2 — Robustness harness  [PASS-WITH-NIT N-1]

**Coverage of every form — verified.** The fixture tree (decisions.md L247-264)
exercises: regular sibling import (`beta.py: from .alpha`), init facade
re-exports (`pkg/__init__.py`), nested subpackage init (`sub/__init__.py`),
**multi-dot** (`sub/__init__.py: from ..alpha`), **bare** (`pkg/__init__.py:
from . import beta`), and a **test→production** edge for D-3. This covers the
two forms with NO live witness (multi-dot, bare-dot) — the exact gap map.md open
question (d) flagged. R3 is genuinely retired by a test, not an assertion.

**Expected edge set internally consistent — spot-checked two non-obvious pairs
against the corrected resolver semantics:**

- `src/pkg/__init__.py -> src/pkg/sub/__init__.py` (from `from .sub import S`):
  base `pkg.sub`, candidates `pkg.sub.S`/`pkg.sub`; `resolve_to_repo_module`
  trims to `pkg.sub` = the subpackage init. ✓
- `src/pkg/__init__.py -> src/pkg/beta.py` (from `from . import beta`):
  candidates `pkg.beta`/`pkg`; `pkg.beta` in index. ✓

The seven-pair oracle (L290-296) matches the corrected resolver's output.

**Circularity assessment — the gate's BLOCK trigger — NOT triggered.** The
prompt asks whether "golden finalized at implement time by running the corrected
extractor once" is a circular oracle. Reading the design fully:

- The expected edges ARE hand-specified from Python import semantics **in the
  design itself** (decisions.md L290-296), independent of any extractor run.
- R5 (L346-350) already mandates: hand-verify each pair against fixture source
  AND confirm a **pre-fix** run fails the same oracle.

So the trustworthy path (author from semantics first; the run only confirms) is
present. This is **not** an untrustworthy oracle and is **not** a BLOCK.

**N-1 (nit).** The lead phrasing at L299-302 ("Exact set to be finalized at
implement time by running the corrected extractor once and hand-verifying") puts
the extractor run *before* hand-verification in the sentence order, which reads
as "author the oracle from the code under test." This contradicts the safer R5
procedure. Tighten the normative instruction to: **"The expected edge set is
hand-specified from Python import semantics BEFORE any extractor run (it already
is, L290-296). The corrected extractor run only CONFIRMS the hand-authored set;
it never AUTHORS it. A pre-fix run MUST fail the same hand-authored oracle."**
This makes the oracle's provenance unambiguous and the regression proof airtight.

**Other D-2 anchors verified:** `skills/architecture-treemap/tests/` does not
exist (map.md L76 — confirmed: no such dir). `extract_components.py` is a script,
not a package — subprocess-invocation strategy is appropriate.

---

## D-3 — `--strict` gate hardening  [PASS]

**Anchors verified in `validate_components.py`:**

- `validate_referential_integrity` at `:374`; connectivity loop
  `logical_nodes_with_edges` at `:421-428` (verified verbatim — no filter on
  edge source kind, exactly as the design and map describe).
- `logical_meta` built at `:431-435`, **after** the connectivity loop — the
  design's Edit 1 note "move its construction above the connectivity loop" is
  necessary and correct.
- Per-node exemption loop at `:437-455`; `is_removable`/`is_leaf` at `:441-442`;
  the `lc_id in logical_nodes_with_edges` test the design replaces is at `:438`
  (verified). Soft-warning vs error split at `:443-455`.
- `--strict` exits 1 only on `ref_issues` at `:536-537` (verified). Since
  `issues` is the first tuple element returned (`:457`) and becomes `ref_issues`,
  appending `production-orphan` to `issues` makes `--strict` fail with **no
  change to main** — the design's claim is correct.

**Would it have flagged the real defect? — verified against the snapshot.**
`run-metadata` (logical node, `:705`) has `classification: seam` (`:713`), NOT
removable, NOT leaf. Its ONLY logical edge is `test-suite -> run-metadata`
(`:4495-4503`, rolled up from `tests/test_run_metadata_contract.py`). Today it is
in `logical_nodes_with_edges` (as a `to`) and non-exempt → **passes the current
gate silently**, exactly the blind spot. Under the hardened gate the only edge is
test-sourced → `production-orphan` → `--strict` fails. Correct, by direct
evidence.

**`is_test_source_node` misclassification check — interrogated:**

- `test-suite` physical_files (`:923-997+`) are ALL under `tests/` (incl.
  `tests/__init__.py`, `tests/conftest.py`, `tests/fixtures/*.txt`,
  `.html`, `.json`). The predicate `all(p.startswith("tests/") or "/tests/" in p
  or basename.startswith("test_"))` returns True via the `tests/` clause for
  every path (the non-`test_`-prefixed fixtures are caught by the prefix clause).
  ✓ The `all()` quantifier means a production node with a single `test_`-named
  file is NOT a test node — the design's R6 mitigation is sound.
- A production node with only test IMPORTERS but which is itself production: the
  predicate keys on the node's OWN `physical_files`, not on who imports it, so a
  production node is never misclassified as a test source. ✓
- Mixed node (some test, some production files): `all()` is False → treated as
  production source → its edges count toward production connectivity. Correct
  (conservative — won't suppress a real production edge).

**False-positive sweep over the existing snapshot — no spurious flags:**

- Core nodes (`data-loader`, `html-generator`, `similarity`, etc.) have many
  production edges → production-connected.
- `package-init` (`:901`, classification core) has inbound production edges at
  `:4281` and `:4362` (the mis-resolved `api/__init__.py` and
  `rendering/__init__.py` collapsing to root). After D-1 those vanish but it
  gains real OUTBOUND sibling edges → stays connected. No flip to orphan.
- Seam nodes `plugin-payload`, `projection-method`, `dashboard-sink`,
  `public-api-surface` all have production edges (verified at `:4000-4002`,
  `:4046`, `:4118`, `:4190-4200`, `:4262`, `:4271`) → connected.

**The one node that needs Phase-C attention — see N-2.** `channel-owners`
(`:636`) has EXACTLY ONE edge in the whole manifest: the audit-asserted
`public-api-surface -> channel-owners` (`:3972-3973`) — which is precisely the
edge D-4 drops. It receives no test edge either. The design's logic chain
handles this correctly (D-1 restores the STATIC `api/__init__.py:39 .channels`
edge, which rolls up to `public-api-surface -> channel-owners` static; D-4 then
drops the redundant audit-asserted twin), so net connectivity is preserved. But
this is a non-obvious dependency: if regeneration drops the audit edge while the
rollup fails to produce the static logical edge, channel-owners would flip to a
production-orphan and (correctly) fail the hardened gate. Flagged as N-2.

**R7 / phase ordering:** D-3 lands before Phase C re-validates the corrected
repo; by then D-1 has restored the facade edges. The design does not spuriously
fail the corrected repo. Sound.

---

## D-4 — Backward-compat + completeness  [PASS]

**Exactly one redundant edge — verified against the snapshot.** There are
exactly THREE `audit-asserted` edges in the witness manifest:

| Edge | type | evidence | Verdict |
|------|------|----------|---------|
| `public-api-surface -> channel-owners` (`:3972`) | `direct-call` | `api/__init__.py:40 -> from .channels …` (`:3977`) | **DROP** — a real static import edge the resolver bug lost; D-1 restores it. Correct. |
| `link-outs -> metadata-assembler` (`:3981`) | `background-knowledge` | `link_outs.py:81-84 … prose only` (`:3986`) | **RETAIN** — not a static import; no Python edge exists. Correct. |
| `link-outs -> render-core` (`:3991`) | `shared-state` | JS `selection.py` + `payload.py` browser state (`:3996`) | **RETAIN** — not a Python import (browser state). Correct. |

The design's D-4 (decisions.md L557-564) identifies exactly the right one to
drop and the right two to keep, distinguished by edge `type` (only the
`direct-call` one is a structural import). **Confirmed correct.**

**Minor line-number note (non-load-bearing).** The committed audit edge's
evidence string says `api/__init__.py:40 -> from .channels` but the actual source
(`api/__init__.py`) has `.channels` on **line 39** and `.payload` on line 40
(verified by reading the file). This is a typo in the *human-authored audit edge*,
faithfully quoted by the design. It does not affect the drop decision (D-4 keys
on endpoints + type, not the line number), and the map's own trace (map.md L142)
correctly uses line 39. No action required; noted for completeness.

**All six map open questions resolved or routed:**

- (a) regression strategy → D-2 (both layers). ✓
- (b) `--strict` blind-spot scope → in scope per user decision 1, designed in D-3. ✓
- (c) backward-compat → D-4 (new dated snapshot; old immutable). ✓ Matches the
  fixed user decision.
- (d) multi-dot/bare → D-1 (correct by construction) + D-2 fixture enforcement. ✓
- (e) non-`__init__` / non-src layouts → confirmed NO effect for src; non-src
  routed to **DEF-3** with src-layout default. ✓
- (f) loud-fail on unresolved base → routed to **DEF-1** with keep-trim default +
  optional diagnostic. ✓

All deferreds carry explicit recommended defaults (DEF-1/2/3). None blocks.

**Phase ordering rationale sound.** A before C (C re-runs corrected extractor);
B before C (C re-validates under hardened `--strict`, which must pass for the
corrected repo). A and B independent except the shared D-2 fixture. Verified
consistent with the channel-owners / run-metadata dependencies above.

**Monolith constraint honored.** All D-1 edits sit within the self-contained
import subsystem `extract_components.py:439-618`; one new pure helper + a
signature extension + one call-site update. No restructure. D-3 edits confined to
`validate_components.py`. Surgical, as required.

---

## Internal consistency

- Line anchors in decisions.md match the live source AND the map's anchors
  (cross-checked L468-469, L493, L510, L528-536, L533, L585, L595-596 in D-1;
  L374, L421-428, L431-435, L437-455, L536-537 in D-3). No drift.
- Component/edge names are consistent across map and design
  (`run-metadata`/`run_metadata.py`, `public-api-surface`, `channel-owners`).
- Each decision's implementation contract is concrete enough to execute without
  re-deriving (exact edits, exact insert points, grep-verified caller count).

---

## What must change before /plan (nits — none block the gate)

1. **N-1 (D-2 oracle phrasing).** Rewrite decisions.md L299-302 so the normative
   procedure is unambiguous: the expected edge set is hand-specified from Python
   import semantics BEFORE any extractor run (it already is at L290-296); the
   corrected extractor run only CONFIRMS, never AUTHORS, the oracle; a pre-fix run
   MUST fail the hand-authored set. This removes the only sentence that reads as
   authoring-the-oracle-from-the-code-under-test.

2. **N-2 (D-3/D-4 Phase-C gate).** Add an explicit Phase-C verification step:
   after dropping the audit-asserted `public-api-surface -> channel-owners` edge,
   confirm the regenerated snapshot contains the rolled-up STATIC
   `public-api-surface -> channel-owners` edge BEFORE committing, so
   `channel-owners` (whose sole edge today is the dropped audit edge) does not
   flip to a `production-orphan`. This is the one node whose connectivity hinges
   on the D-1-restore and D-4-drop being atomic in the same regeneration.

3. **N-3 (optional, cosmetic).** The helper's `relative_path == "__init__.py"`
   clause (decisions.md L143) is unreachable — a top-level `__init__.py` yields
   `source_module is None`, so the `if not source_module: return []` guard fires
   first. Either drop the clause or note it as defensive. No behavioural impact.

These are clarifications, not corrections to the design's logic. The design is
ready to proceed to /plan once N-1 and N-2 are folded in (N-3 optional).
