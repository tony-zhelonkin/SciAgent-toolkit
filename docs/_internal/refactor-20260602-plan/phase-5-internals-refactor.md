# Phase 5: SciAgent-toolkit Internals Refactor

## Summary

This phase cleans up the brittle, redundant internal logic inside `lib/sciagent/`. 
Phases 1–4 move project-structure and scaffolding responsibilities out of scbio-docker and into SciAgent-toolkit; that work makes `sciagent` the single authority for project-level context management, which means its internals now carry more load and must be sound. The four worst offenders are: (a) `stack.sh:stack_walk()`, which open-codes the same insertion-order / last-wins / shadow-tracking algorithm three times for SKILL/AGENT/COMMAND; (b) the `inject.sh` ↔ `collisions.sh` contract, where `inject` writes symlinks and rewrites the managed block without ever consulting collision detection, so collisions surface post-hoc in `validate`/`status` instead of at the point of action; (c) `inject.sh:_inject_rewrite_block()`, which is a thin wrapper that is fine in principle but sits on top of a managed-block I/O surface that must be the exclusive property of `block.sh`; and (d) cross-cutting concerns — the per-verb module load graph in `bin/sciagent` and the inconsistent `return`/`exit`/silent-failure error conventions.

The refactor is behavior-preserving for the happy path (a stack walk produces the same tagged table; an inject produces the same symlinks and block) but adds one new user-facing behavior: inject now performs an inline collision check before it writes. All changes are internal to `lib/sciagent/` and `bin/sciagent`; no role/skill/agent/command content changes. Breaking changes are acceptable (single user, no downstream consumers), so we optimize for clean modularity over compatibility shims.

## Current State

All paths below are relative to `toolkits/SciAgent-toolkit/`.

- `lib/sciagent/stack.sh` (8773 bytes). Two public functions:
  - `stack_walk()` — ~95 lines. Holds six associative arrays (`_sw_skills`, `_sw_agents`, `_sw_commands` plus three `*_shadows` mirrors) and three order arrays (`_sw_skill_order`, `_sw_agent_order`, `_sw_command_order`), then runs a `case "$kind"` with three near-identical `SKILL)`/`AGENT)`/`COMMAND)` arms. Each arm is the same eight-line "if already present, append previous provider to the shadow CSV; else push to the order array; then record this role as the current provider" logic, differing only in which array triple it touches. `OUTPUT_STYLE` is a separate fourth, simpler arm (single-value, last-wins, one shadow slot). Emission is four near-identical `for n in "${_sw_*_order[@]}"` loops.
  - `render_block_body()` — consumes `stack_walk` output via `IFS=$'\t' read` and rebuilds a *second* set of six arrays (`_rb_*`) to render the Markdown block body. This duplicates the array-rehydration shape but is a legitimate consumer, not a target for 5.1.
- `lib/sciagent/inject.sh` (15661 bytes; the prompt's snapshot was ~411 lines). Key functions: `cmd_inject`, `_inject_detect_kind`, `_inject_named_entry`, `_inject_by_tag`, `_inject_entry_is_stack_mounted`, `_inject_rewrite_block`.
  - `_inject_named_entry()` resolves a canonical source, checks stack-mount, checks injected-idempotency, creates dual symlinks (`.claude/` + `.agents/`), appends to the manifest, optionally synthesizes `_injected` overlay, then calls `_inject_rewrite_block`. **It never calls anything from `collisions.sh`.**
  - `_inject_rewrite_block()` reads the manifest stack, collects injected skill names, calls `render_block_body` (from `stack.sh`), then calls `block_write AGENTS.md "$body"` and `manifest_update_block_hash "$(block_stored_hash AGENTS.md)"` (both from `block.sh`). So inject already *delegates* managed-block I/O to `block.sh` — see 5.3 for the nuance.
- `lib/sciagent/collisions.sh` (5321 bytes). Public `collisions_enumerate()` walks all four canonical namespaces (skills/agents/commands/roles) and emits `<name>\t<kind-csv>` for every basename present in ≥2 namespaces. Pure read-only enumeration of the **toolkit**, with no notion of the **active stack** or of "what an inject would add." Consumers today: `validate.sh` (soft-warn) and `status.sh` (Notes section).
- `lib/sciagent/block.sh` (5328 bytes). The managed-block authority: `block_read`, `block_stored_hash`, `block_hash_check`, `block_write`, `block_remove`, plus private `_sha1`/`_count_lines`. Owns the marker constants `_BLOCK_BEGIN_PREFIX` / `_BLOCK_END` and the hash-canonicalization invariant.
- `lib/sciagent/new.sh` (3655 bytes). `cmd_new` with `_new_project` (copies only 3 templates — flagged as too thin in Phase 1/2, out of scope here), `_new_role`, `_new_skill`, `_new_agent`.
- `bin/sciagent` (4045 bytes, 137 lines). Lazy per-verb sourcing. Observed graph:
  - `activate` / `deactivate`: block, symlinks, roles, skill_deps, **collisions**, validate, stack, claude_settings, activate (+deactivate).
  - `inject`: block, symlinks, roles, skill_deps, stack, inject. **No collisions, no validate.**
  - `eject`: block, symlinks, roles, skill_deps, stack, inject, eject.
  - `status`: block, symlinks, roles, skill_deps, stack, **collisions**, status.
  - `list`: roles, skill_deps, status.
  - `validate`: skill_deps, **collisions**, validate.
  - `new`: new.
- Error-handling conventions are mixed: library functions mostly `return 1`; `bin/sciagent` top-level uses `exit 1`; several spots swallow failures (`|| true`, `2>/dev/null`, `manifest_*` calls whose return is unchecked). `set -u` is on in `bin/sciagent`; **`set -e` / `set -o pipefail` are not**, and library functions deliberately rely on `return`-code propagation rather than `errexit`.

## Changes

### 5.1 Extract `_sw_record()` from `stack_walk()`

**Files:** `lib/sciagent/stack.sh`

**What:** Replace the three duplicated `SKILL)`/`AGENT)`/`COMMAND)` arms in `stack_walk()` with three calls to a single shared helper, `_sw_record()`, that encapsulates the "insertion-order + last-wins + append-prev-to-shadow-CSV" pattern. `OUTPUT_STYLE` keeps its own small arm (different shape: scalar, not a map). Emission stays four explicit blocks (they differ only in the leading tag string, but keeping them explicit is clearer than a generic loop over a map-of-maps in bash).

**Why:** The user called the current code "terrible nested boolean logic of if-else arrow code." The three arms are byte-for-byte identical modulo the array-variable prefix, so any fix or audit currently has to be made in triplicate. One helper means one place to reason about the shadow-CSV semantics.

**How:**

1. Define the helper. Because bash cannot pass associative arrays by value, `_sw_record` takes the array **names** and uses namerefs (`local -n`, available in bash 4.3+, which the codebase already assumes via `declare -A`). Signature:

   ```bash
   # _sw_record <map-name> <shadow-map-name> <order-array-name> <name> <role>
   # Records <name>=<role> into the provider map with last-wins semantics.
   # On a repeat name, the *previous* provider is appended to the shadow CSV
   # (preserving the existing comma-joined accumulation order). On a first
   # sighting, <name> is pushed onto the order array.
   _sw_record() {
       local -n _map="$1"
       local -n _shadow="$2"
       local -n _order="$3"
       local name="$4" role="$5"

       if [[ -n "${_map[$name]:-}" ]]; then
           local prev="${_map[$name]}"
           if [[ -n "${_shadow[$name]:-}" ]]; then
               _shadow[$name]="${_shadow[$name]},$prev"
           else
               _shadow[$name]="$prev"
           fi
       else
           _order+=("$name")
       fi
       _map[$name]="$role"
   }
   ```

2. Rewrite the dispatch loop body in `stack_walk()` to delegate:

   ```bash
   case "$kind" in
       SKILL)
           _sw_record _sw_skills   _sw_skill_shadows   _sw_skill_order   "$name" "$role" ;;
       AGENT)
           _sw_record _sw_agents   _sw_agent_shadows   _sw_agent_order   "$name" "$role" ;;
       COMMAND)
           _sw_record _sw_commands _sw_command_shadows _sw_command_order "$name" "$role" ;;
       OUTPUT_STYLE)
           [[ -n "$_sw_style" ]] && _sw_style_shadow="$_sw_style_role"
           _sw_style="$name"; _sw_style_role="$role" ;;
   esac
   ```

3. Leave the declarations of the six arrays + three order arrays and the four emission loops unchanged. The emission loops already read cleanly; do not over-refactor them.
4. Nameref hygiene: prefix the helper's locals with `_` (`_map`, `_shadow`, `_order`) and use distinct names from the caller's arrays so a nameref can never alias one of the helper's own locals (a known bash nameref footgun). The caller arrays are `_sw_*`; the helper locals are unprefixed-`_map` etc., so there is no collision.
5. Add/extend a unit test in `tests/` that activates a base+overlay pair where the overlay shadows a base skill, an agent, and a command simultaneously, and asserts the `stack_walk` TSV has the correct `<shadowed-roles-csv>` column for all three kinds. This locks the behavior the helper must preserve (multi-shadow CSV accumulation order in particular).

### 5.2 Inject performs an inline collision check before writing

**Files:** `lib/sciagent/inject.sh`, `lib/sciagent/collisions.sh`, `bin/sciagent`

**What:** Before `_inject_named_entry` writes any symlink or touches the manifest, compute whether the entry being injected would create (or sit on) a cross-namespace name collision, and surface it inline. Add a small, stack-aware predicate to `collisions.sh` so `inject` does not have to re-derive the four-namespace enumeration. Wire `collisions.sh` into the `inject` (and `eject`) verb in `bin/sciagent`.

**Why:** Today a user runs `sciagent inject foo`, it succeeds silently, and only a later `sciagent status` or `sciagent validate` reveals that `foo` now collides (e.g. exists as both a skill and a command). The error manifests far from its cause — user-hostile. Inject is the action that *introduces* the load-bearing ambiguity, so inject is where the warning belongs. (Note: `_inject_detect_kind` already hard-fails on *toolkit-side* ambiguity for the bare `inject <name>` form; this change is different — it concerns the *cross-namespace collision* that becomes load-bearing in the project's mounted trees, including the explicit-flag forms that bypass `_inject_detect_kind`.)

**How:**

1. Add a focused helper to `collisions.sh` that answers the single question inject needs, reusing `collisions_enumerate` rather than re-walking:

   ```bash
   # collisions_for_name <name>
   # Emits the kind-csv (e.g. "skill,command") iff <name> collides across
   # ≥2 namespaces in the toolkit; empty output + return 1 otherwise.
   collisions_for_name() {
       local name="$1" line nm csv
       while IFS=$'\t' read -r nm csv; do
           if [[ "$nm" == "$name" ]]; then
               printf '%s\n' "$csv"
               return 0
           fi
       done < <(collisions_enumerate)
       return 1
   }
   ```

2. In `_inject_named_entry`, after the source is resolved and the kind is known, but **before** the symlink/`mkdir`/`manifest_append_inject` block, run the check:

   ```bash
   local _coll_csv
   if _coll_csv=$(collisions_for_name "$name"); then
       _inject_handle_collision "$name" "$kind" "$_coll_csv" || return $?
   fi
   ```

   where `_inject_handle_collision` applies the policy chosen in ADR-5.1 (warn-and-proceed by default; `--force` suppresses; optional strict mode fails). Keep this in a single helper so the policy lives in one place.

3. Wire the module in `bin/sciagent`: add `. "$LIB/collisions.sh"` to the `inject` arm (and the `eject` arm, which sources `inject.sh` and may need the same guard on re-mount paths). This is the load-graph fix for 5.4 as it pertains to inject.
4. Ordering guarantee: the check is **read-only** and runs before the first side-effecting line (`mkdir -p`). On a fail policy, nothing has been written, so no rollback is needed. On a warn policy, the warning is emitted, then writes proceed exactly as before. Re-run idempotency is unaffected because the stack-mount and injected-idempotency guards still run first (a second `inject foo` short-circuits before reaching the collision check, so the warning is not re-emitted on no-op re-injects — desirable).
5. Allowlist integration: `collisions.sh`'s header references `tests/collision-allowlist.txt` as the source of truth for approved "family overlaps." `_inject_handle_collision` must consult that allowlist and downgrade an allowlisted collision to either silence or a single informational line (not a warning). Add a `collisions_is_allowlisted <name> <csv>` reader to `collisions.sh` so both inject and validate share one allowlist parser.

### 5.3 `block.sh` is the sole managed-block authority; inject is a pure caller

**Files:** `lib/sciagent/inject.sh`, `lib/sciagent/block.sh` (audit only), `lib/sciagent/eject.sh` (audit only)

**What:** Confirm and enforce the invariant that *all* managed-block byte manipulation (marker constants, BEGIN/END detection, hash canonicalization, write/remove) lives only in `block.sh`. Audit `inject.sh` (and `eject.sh`) for any logic that re-implements, rather than calls, that surface. The prompt flagged "YAML block-rewriting logic that duplicates `block.sh`."

**Why:** Single source of truth for the managed block prevents hash-canonicalization drift (the `block_write` invariant about trailing-newline-before-hashing is subtle and must exist in exactly one place) and prevents two code paths writing subtly different marker bytes.

**How:**

1. Audit finding (from the snapshot): `_inject_rewrite_block()` does **not** re-implement block I/O — it correctly calls `render_block_body` (body generation, owned by `stack.sh`) then `block_write` + `block_stored_hash` + `manifest_update_block_hash` (all owned by `block.sh`). So the *primary* duplication the prompt worried about is **not present** in `inject.sh`'s block path. The real overlap is conceptual: `_inject_rewrite_block` and the analogous "rebuild the block after a stack change" routine in `activate.sh` (and the teardown routine in `eject.sh`/`deactivate.sh`) each open-code the same 4-step recipe (read stack → collect injected names → `render_block_body` → `block_write` + refresh hash).
2. Extract that recipe into one shared function so inject/activate/eject call it instead of each maintaining a copy:

   ```bash
   # block_render_and_write <agents-file>
   # Canonical "recompute the managed block from current manifest + role YAMLs
   # and write it" routine. Lives in stack.sh (it needs render_block_body) but
   # delegates ALL byte I/O to block.sh. Returns block.sh's write status.
   block_render_and_write() {
       local file="${1:-AGENTS.md}"
       local stack base overlay
       stack=$(manifest_stack)
       base=$(printf '%s\n' "$stack" | awk '{print $1}')
       overlay=$(printf '%s\n' "$stack" | awk '{print $2}')

       local -a injected=()
       local ov nm _via kind
       while IFS='|' read -r ov nm _via kind; do
           [[ -n "$ov" ]] || continue
           [[ "${kind:-skill}" == "skill" ]] || continue
           injected+=("$nm")
       done < <(manifest_injected)

       local body
       body=$(render_block_body "$base" "$overlay" "${injected[@]+"${injected[@]}"}")
       block_write "$file" "$body" || return $?
       manifest_update_block_hash "$(block_stored_hash "$file")"
   }
   ```

3. Replace the body of `_inject_rewrite_block()` with a one-liner `block_render_and_write AGENTS.md` (or delete `_inject_rewrite_block` entirely and call `block_render_and_write` at the call site). Do the equivalent substitution in `activate.sh`/`eject.sh` wherever they currently open-code the same recipe.
4. Enforce the marker-constant boundary: grep the tree for `_BLOCK_BEGIN_PREFIX`, `_BLOCK_END`, `BEGIN SCIAGENT:ROLES`, and `sha1sum` outside `block.sh`. Any hit outside `block.sh` is a layering violation; route it through a `block.sh` accessor. (Expected result: zero hits today, but add this as a guard test so it stays zero.)
5. Keep `render_block_body` in `stack.sh` (it is body *content* generation, not block *framing*). The clean seam is: `stack.sh` owns "what the body says," `block.sh` owns "how the body is framed and hashed in the file." `block_render_and_write` is the thin orchestrator that joins them and is the only thing inject/activate/eject call.

### 5.4 Tighten the per-verb module load graph

**Files:** `bin/sciagent`

**What:** Audit the lazy per-verb `source` lists for missing, redundant, and ordering-sensitive entries; make each verb load exactly its transitive closure, no more, no less.

**Why:** The graph is currently maintained by hand and has at least one correctness gap 
(`inject` lacks `collisions.sh`, fixed in 5.2) and some apparent over-loading. 
A wrong graph means either a runtime "command not found" (under-load) or wasted source time and accidental coupling (over-load).

**How:**

1. Correctness gaps to close:
   - `inject`: add `collisions.sh` (required by 5.2). After 5.3, `inject` calls `block_render_and_write` which lives in `stack.sh` (already sourced) and uses `render_block_body` (also `stack.sh`) — graph stays consistent.
   - `eject`: add `collisions.sh` if 5.2's re-mount guard touches it; otherwise leave.
2. Over-load to investigate (verify with a call-graph grep before removing — removing a needed source is a hard runtime break):
   - `activate`/`deactivate` source `validate.sh`. Confirm whether activation actually invokes a `validate.sh` function or whether validation is a separate user-invoked verb; if unused, drop it from these arms.
   - `list` sources `status.sh` to reach `cmd_list`. If `cmd_list` is small and self-contained, consider moving it to its own `list.sh` so `list` stops dragging in the entire `status.sh` surface.
3. Circular dependency check: there are no `source` statements inside `lib/sciagent/*.sh` (all sourcing is centralized in `bin/sciagent`), so true source-time cycles are impossible by construction. The only "cycle" risk is *call-time* (function A in module X calls function B in module Y and vice versa); document the intended layering as a one-directional DAG:

   ```
   block.sh        (leaf: file I/O, no sciagent deps)
   roles.sh        (leaf: YAML reads)
   skill_deps.sh   (depends on roles.sh)
   stack.sh        (depends on roles.sh; render_block_body depends on block.sh via the orchestrator)
   collisions.sh   (leaf: toolkit enumeration; reads allowlist)
   symlinks.sh     (depends on roles.sh)
   validate.sh     (depends on skill_deps.sh, collisions.sh)
   claude_settings.sh (leaf-ish)
   inject.sh       (depends on roles.sh, stack.sh, block.sh, symlinks.sh, collisions.sh)
   activate.sh     (depends on ~everything above)
   eject.sh        (depends on inject.sh + stack.sh + block.sh)
   status.sh       (depends on stack.sh, collisions.sh, roles.sh)
   ```

4. Codify the graph: replace the per-verb hand-maintained `source` lists with a small declarative table at the top of `bin/sciagent` (verb → space-separated module list) plus a loop that sources `block.sh` first (no deps) and the verb's listed modules in dependency order. This makes the closure auditable in one place and is where the 5.2/5.4 edits land permanently.
5. Add a CI smoke test that runs every verb with `--help`/no-op args in a scratch project and asserts no "command not found" — the cheapest possible guard against an under-loaded arm.

### 5.5 Standardize error handling and propagation

**Files:** all of `lib/sciagent/*.sh`, `bin/sciagent`

**What:** Adopt one convention for how library functions signal failure and how the dispatcher surfaces it.

**Why:** The current mix of `return 1`, top-level `exit 1`, and silently-swallowed failures (`|| true`, unchecked `manifest_*` writes) makes failure modes unpredictable — a failed manifest write or block write can leave the project half-mounted with exit code 0.

**How:**

1. The rule: **library functions (`lib/sciagent/*.sh`) only ever `return <code>`; never `exit`.** Only `bin/sciagent` (and only at the top dispatch level) may `exit`. This keeps every library function sourceable and testable without it tearing down the caller's shell. (`block.sh` already documents per-function exit *codes* — formalize that every function in every module documents its non-zero return meanings in a header comment, as `block.sh` does.)
2. Propagation: a verb's `cmd_*` entrypoint returns the status of its first failing step; `bin/sciagent` does `cmd_<verb> "$@"; exit $?` so the process exit code mirrors the library return. Replace any mid-function `exit` in libraries with `return`.
3. No silent state-mutating failures. Side-effecting calls — `block_write`, `manifest_append_inject`, `manifest_update_stack`, `manifest_update_block_hash`, `ln -sfn`, `mkdir -p` — must be checked: `cmd || { echo "..." >&2; return 1; }`. Reserve `|| true` for genuinely best-effort, non-state operations (e.g. the best-effort `sed` substitution in `new.sh:_new_skill`) and add a `# best-effort: <reason>` comment on each surviving `|| true`.
4. Keep `set -u` in `bin/sciagent`. Do **not** turn on `set -e`/`set -o pipefail` globally — the codebase intentionally uses return-code dispatch (e.g. `_inject_detect_kind` returns 1 as a control-flow signal, idempotent guards return 0 after refusing) and `errexit` would break those patterns and the many `read`-from-process-substitution loops. Rely on explicit checks instead.
5. Standardize the user-facing error prefix to `sciagent <verb>: <message>` on stderr (already the dominant style in `inject.sh`/`new.sh`); fix the stragglers that use bare `error:` (e.g. `_inject_detect_kind`'s `"error: ambiguous …"`) to match, or deliberately keep `error:`/`note:` as severity tags and document the convention. Capture the chosen convention in `CONTRIBUTING.md`.
6. Add a tiny `tests/` assertion that greps `lib/sciagent/*.sh` for the token `exit ` and fails if any library module contains a top-level `exit` (allowing `exit` only inside `awk`/subshell strings) — a structural guard for rule (1).

## Open ADRs

### ADR-5.1: What is inject's behavior when it would create a cross-namespace collision?
**Options:**
- A — **Warn and proceed.** Emit `sciagent inject: '<name>' now collides across <csv> — disambiguation (/foo vs @foo vs foo) is load-bearing` on stderr, then mount. Exit 0.
- B — **Fail closed.** Refuse the inject (exit 1) unless the collision is on the allowlist; require an explicit flag to override.
- C — **Warn by default, `--force` to silence, optional `--strict` (or `SCIAGENT_STRICT_COLLISIONS=1`) to upgrade the warning to a hard fail.**

**Recommended:** C.
**Why:** The header of `collisions.sh` states that "most collisions in my practice are deliberate 'family overlaps'" — so failing closed (B) would be hostile to the common, intentional case and would fight the existing allowlist mechanism. Pure warn (A) is the right default but gives no escape hatch for scripted/CI runs that want collisions to be fatal. C keeps the friendly default, lets scripts opt into strictness, and lets an interactive user silence an already-understood overlap with `--force`. Allowlisted collisions are downgraded to a single informational line (or silence) under all three modes, consistent with how `validate.sh` already treats the allowlist.
**Blocking implementation:** no — 5.2 can land with the warn-and-proceed core; `--force`/`--strict` are additive flags on `_inject_handle_collision`.

### ADR-5.2: How is the associative array passed into `_sw_record()` — nameref vs. eval vs. global-by-convention?
**Options:**
- A — **`local -n` namerefs** (bash 4.3+).
- B — **`eval`-built array references.**
- C — **Convention: helper mutates well-known globals `_sw_*` directly**, taking only the kind tag as an argument.

**Recommended:** A.
**Why:** The codebase already requires bash 4.x (it uses `declare -A` throughout), so namerefs are available and are the idiomatic, readable way to pass arrays by reference. `eval` (B) reintroduces quoting hazards the refactor is trying to kill. Convention-globals (C) would make `_sw_record` non-reusable and re-couple it to `stack_walk`'s internal variable names — exactly the implicit coupling we are removing. The only caveat is nameref self-aliasing, mitigated by the local-naming hygiene in 5.1 step 4.
**Blocking implementation:** no — A is the default; if a future bash-3.2 target ever appears, fall back to C, but that is not a current constraint.

### ADR-5.3: Does `block_render_and_write` live in `stack.sh` or `block.sh`?
**Options:**
- A — In **`stack.sh`** (it needs `render_block_body`, which lives there).
- B — In **`block.sh`** (it is the "write the managed block" orchestrator).
- C — In a new thin **`render.sh`** orchestration module.

**Recommended:** A.
**Why:** `block.sh` must stay a pure, dependency-free leaf (file I/O + hashing only) so it remains the unambiguous single authority for block framing — pulling `render_block_body`/`manifest_*` knowledge into it would invert the dependency and pollute the leaf. The orchestrator inherently needs `render_block_body` (body content) and `manifest_*` (state), both of which already sit at the `stack.sh` layer, so A keeps the dependency arrow pointing the right way (`stack.sh` → `block.sh`, never the reverse). A new module (C) is over-engineering for a single 12-line function.
**Blocking implementation:** no — placement is mechanical; the function body is identical wherever it lands.

### ADR-5.4: Replace hand-maintained per-verb source lists with a declarative table?
**Options:**
- A — **Keep explicit per-verb `source` blocks** (current), just fix the gaps.
- B — **Declarative verb→modules table** + a sourcing loop (5.4 step 4).

**Recommended:** B.
**Why:** The hand-maintained lists already drifted (inject's missing `collisions.sh`). A single table makes each verb's closure auditable at a glance and makes the dependency DAG (5.4 step 3) enforceable by a test. The cost is one small indirection in `bin/sciagent`. Given Phase 5's whole premise is "internals should be sound and auditable," B aligns better.
**Blocking implementation:** no — 5.2's correctness fix works under either A or B; B is the cleanup we apply once the gaps are known.

## Dependencies
- **Depends on:** Phase 1 (project-structure / scaffolding moved into SciAgent-toolkit, establishing it as the single context-management authority — gives this internals work its mandate). Does not depend on Phases 2–4 at the code level; can proceed in parallel once Phase 1 lands. The `new.sh:_new_project` "too thin" fix is owned by the scaffolding phase, not here.
- **Enables:** Any later phase that adds verbs or extends inject/activate behavior, because it leaves a clean `block_render_and_write` orchestrator, a single `_sw_record` shadow algorithm, an inline-collision contract, and a declarative load graph to build on.

## Breaking Changes
- **New user-visible behavior on `inject`:** a collision warning (and, under ADR-5.1 option C, a possible `--strict` hard-fail and a `--force` flag). Scripts that parse inject's stderr may see new lines. Acceptable per the "single user, no downstream consumers" ground rule.
- **Internal API churn (not user-facing):** `_inject_rewrite_block` is removed/reduced to a call site; `block_render_and_write` is introduced in `stack.sh`; `_sw_record`, `collisions_for_name`, `collisions_is_allowlisted` are introduced. No external caller depends on these private functions.
- **No changes** to role/skill/agent/command content, to the AGENTS.md managed-block format/markers/hash, to the manifest format, or to the symlink layout. A re-`activate` after this phase produces a byte-identical managed block for the same stack.

## Estimated Scope
- `lib/sciagent/stack.sh`: +~15 (`_sw_record`), −~40 (collapsed arms) → net ~−25; plus +~12 for `block_render_and_write`. Net roughly flat, materially simpler.
- `lib/sciagent/collisions.sh`: +~20 (`collisions_for_name`, `collisions_is_allowlisted`).
- `lib/sciagent/inject.sh`: +~15 (collision check + `_inject_handle_collision`), −~10 (`_inject_rewrite_block` body collapses to a call). Net ~+5.
- `lib/sciagent/activate.sh` / `eject.sh`: −~10 each (open-coded block-rebuild recipe → call to `block_render_and_write`).
- `bin/sciagent`: ~+10/−~30 if converting to the declarative table (ADR-5.4 B); +1 line if keeping explicit blocks (option A).
- Error-handling sweep (5.5): touches every `lib/sciagent/*.sh`, mostly 1–3 lines each (checked side-effects, prefix normalization) → ~+40 spread thin.
- `tests/`: +4 new assertions (shadow-CSV multi-kind, marker-constant boundary grep, no-`exit`-in-libs grep, per-verb `--help` smoke).
- **Total:** ~8 source files + tests; net line delta near zero (roughly +60 / −90), with a large reduction in duplicated/branchy logic and a net gain in test coverage.
