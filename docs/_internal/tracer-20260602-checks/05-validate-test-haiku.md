# SciAgent-Toolkit Integrity Audit — 2026-06-02

**Audit Path:** Validate & Test (tracer-bullet mode)  
**Toolkit Location:** `/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit`  
**Test Command:** `bash tests/run-all.sh`  
**Audit Date:** 2026-06-02  
**Auditor:** Haiku 4.5

---

## Executive Summary

**Test Suite Result:** 68 PASS, 1 FAIL (scope-lint regression block)  
**Exit Code:** 1 (failed)  
**Status:** **BLOCKER-CLASS BUG** preventing green CI

- ✅ All 67 bash tests passing
- ❌ 1 skill violates legacy scope vocabulary (atomic → concept/implementation)
- ⚠️ 2 skills missing `last-reviewed` metadata (non-blocking)
- ⚠️ 3 orphaned references in README.md to deleted roles
- ⚠️ 1 silent template mismatch: code expects `templates/project/software-tool/` but directory is `software`
- ✅ No violations: exit() in libs, missing modules, dead modules

---

## Part A: Test Suite Results

### Overall Summary
```
passed: 68
failed: 1
failures:
  - test_skill_scope_lint.sh
```

### Test Result Table

| Test Name | Status | Category | Notes |
|-----------|--------|----------|-------|
| test_activate_aborts_on_unknown_tag.sh | ✅ PASS | activation | |
| test_activate_idempotent.sh | ✅ PASS | activation | |
| test_activate_max_stack.sh | ✅ PASS | activation | Stack depth cap (max 2) enforced |
| test_activate_solo.sh | ✅ PASS | activation | |
| test_activate_stack.sh | ✅ PASS | activation | |
| test_activate_warn_on_missing_complementary.sh | ✅ PASS | activation | |
| test_activate_warns_on_dropped_inject.sh | ✅ PASS | activation | |
| test_block_drift.sh | ✅ PASS | block markers | |
| test_block_marker_boundary.sh | ✅ PASS | block markers | Boundary condition handling |
| test_block_missing_marker.sh | ✅ PASS | block markers | |
| test_block_no_markers_append.sh | ✅ PASS | block markers | |
| test_block_roundtrip.sh | ✅ PASS | block markers | Inject→eject round-trip |
| test_claude_settings_lifecycle.sh | ✅ PASS | claude settings | .claude/ state management |
| test_collision_allowlist_blocks_unknown.sh | ✅ PASS | collisions | Rejects unlisted cross-namespace overlaps |
| test_collision_allowlist_passes_known.sh | ✅ PASS | collisions | Allows known overlaps (architect, architecture-treemap) |
| test_deactivate_full.sh | ✅ PASS | deactivation | |
| test_deactivate_partial.sh | ✅ PASS | deactivation | Tear down one role from stack |
| test_eject_agent.sh | ✅ PASS | injection/ejection | |
| test_eject_command.sh | ✅ PASS | injection/ejection | |
| test_eject_kind_forward_compat.sh | ✅ PASS | injection/ejection | Handles pre-PR3 manifest format |
| test_eject_named.sh | ✅ PASS | injection/ejection | |
| test_eject_pre_pr3_manifest.sh | ✅ PASS | injection/ejection | Backward compat with old manifest |
| test_eject_rejects_stack_mounted_command.sh | ✅ PASS | injection/ejection | Cannot eject stack-mounted items |
| test_eject_rejects_stack_mounted.sh | ✅ PASS | injection/ejection | |
| test_eject_rejects_stack_mounted_with_comments.sh | ✅ PASS | injection/ejection | |
| test_eject_tag.sh | ✅ PASS | injection/ejection | Eject all by tag |
| test_inject_agent_basic.sh | ✅ PASS | injection/ejection | |
| test_inject_ambiguous_hard_fail.sh | ✅ PASS | injection/ejection | Name collision disambiguation |
| test_inject_collision_warns.sh | ✅ PASS | injection/ejection | |
| test_inject_command_basic.sh | ✅ PASS | injection/ejection | |
| test_inject_command_with_companion_skill.sh | ✅ PASS | injection/ejection | |
| test_inject_creates_overlay.sh | ✅ PASS | injection/ejection | Creates overlay role if needed |
| test_inject_eject_roundtrip.sh | ✅ PASS | injection/ejection | Inject→eject round-trip |
| test_inject_explicit_flag_resolves.sh | ✅ PASS | injection/ejection | --skill/--agent/--command flags work |
| test_inject_extends_overlay.sh | ✅ PASS | injection/ejection | |
| test_inject_rejects_stack_mounted.sh | ✅ PASS | injection/ejection | Cannot inject stack-mounted |
| test_inject_tag_creates_overlay.sh | ✅ PASS | injection/ejection | Inject by tag creates overlay |
| test_inject_unknown_name.sh | ✅ PASS | injection/ejection | Rejects unknown names |
| test_list_deps_verbs.sh | ✅ PASS | listing | |
| test_manifest_ownership.sh | ✅ PASS | manifest | |
| test_new_project_types.sh | ✅ PASS | scaffolding | Creates analysis & software-tool projects |
| test_no_duplicate_basenames.sh | ✅ PASS | naming | Cross-namespace collision detection |
| test_no_exit_in_libs.sh | ✅ PASS | code quality | Libraries use `return`, not `exit` |
| test_results_gitignore_layout.sh | ✅ PASS | gitignore | |
| test_skill_frontmatter_valid.sh | ✅ PASS | skill metadata | YAML frontmatter validates |
| test_skill_requires_cycle.sh | ✅ PASS | dependency graph | Detects cycles |
| test_skill_requires_deep_chain.sh | ✅ PASS | dependency graph | Transitive closure |
| test_skill_requires_diamond.sh | ✅ PASS | dependency graph | Diamond dependencies resolved |
| test_skill_requires_missing.sh | ✅ PASS | dependency graph | Detects missing requires targets |
| test_skill_requires_resolution.sh | ✅ PASS | dependency graph | Resolution semantics |
| **test_skill_scope_lint.sh** | ❌ **FAIL** | **scope enforcement** | **Legacy scope 'atomic' in mllmcelltype-consensus-annotation** |
| test_stack_walk_multi_kind_shadow.sh | ✅ PASS | stack walks | |
| test_status_injected_agent_classified_correctly.sh | ✅ PASS | status | |
| test_status_no_notes_when_clean.sh | ✅ PASS | status | |
| test_status_shows_collisions.sh | ✅ PASS | status | |
| test_status_shows_inherited.sh | ✅ PASS | status | |
| test_tags_vocabulary.sh | ✅ PASS | tags | Valid tag vocabulary |
| test_validate_collision_quiet.sh | ✅ PASS | validation | |
| test_validate_collision_three_kinds.sh | ✅ PASS | validation | Detects 3-way collisions |
| test_validate_collision_warn.sh | ✅ PASS | validation | |
| test_validate_cycle.sh | ✅ PASS | validation | Detects cycles |
| test_validate_missing_requires.sh | ✅ PASS | validation | |
| test_validate_no_false_positives.sh | ✅ PASS | validation | |
| test_validate_output_style_drift.sh | ✅ PASS | validation | |
| test_validate_quiet_flag.sh | ✅ PASS | validation | |
| test_validate_skills_ref_absent.sh | ✅ PASS | validation | Checks for undefined skills |
| test_validate_unknown_tag.sh | ✅ PASS | validation | Rejects unknown tags |
| test_verb_smoke_load_graph.sh | ✅ PASS | smoke test | Loads CLI dependency DAG |
| [mllmct] smoke_check_versions | ✅ PASS | packaged skill | Version pinning + seams intact |

### Detailed Failure Analysis

#### ❌ test_skill_scope_lint.sh

**Exit Code:** 1  
**Output:**
```
FAIL [test_skill_scope_lint.sh] mllmcelltype-consensus-annotation: legacy scope 'atomic' — rename to concept or implementation (ADR-003)
FAIL [test_skill_scope_lint.sh] 1 skill(s) over cap (post-cutoff); 0 legacy warning(s); 69 ok
```

**Root Cause:** The skill `skills/mllmcelltype-consensus-annotation/SKILL.md` (lines 1–30) carries:
```yaml
metadata:
  scope: atomic       # ← LEGACY VOCABULARY (deprecated as of ADR-003, 2026-05-24)
  last-reviewed: 2026-05-29
```

**Why This Fails:**
- `test_skill_scope_lint.sh` enforces a regression block (lines 81–88) that rejects any skill carrying scope keywords from the old vocabulary: `atomic|orchestrator|foundation`.
- These were renamed to `concept` (≤500 body lines) or `implementation` (≤350 body lines) per ADR-001 and ADR-003 (kickoff.md §9, 2026-05-24).
- The `last-reviewed` date is **at or after** the cutoff (2026-05-29 ≥ 2026-05-24), so the legacy exemption does NOT apply.

**Is This a Real Bug?**  
✅ **YES — BLOCKER.** The skill uses forbidden vocabulary explicitly disallowed by the test suite to prevent accidental re-introduction of deprecated scope categories. This is a hard fail by design.

**Fix Required:**
- Change `scope: atomic` to `scope: concept` or `scope: implementation` in `skills/mllmcelltype-consensus-annotation/SKILL.md` (line 5).

---

## Part B: Validation Command

### sciagent validate (non-quiet)
```
validate: warning — name 'architect' appears as agent, command, and role (mounting both is supported; ensure the overlap is intentional)
validate: warning — name 'architecture-treemap' appears as both skill and command (mounting both is supported; ensure the overlap is intentional)
sciagent validate: all checks passed
EXIT CODE: 0
```

**Status:** ✅ PASS (warnings are expected and in collision-allowlist.txt)

### sciagent validate --quiet
```
EXIT CODE: 0
```

**Status:** ✅ PASS (no output, silent success)

---

## Part C: Static Integrity Sweeps

### C.1: Module & Code Quality

#### ✅ No exit() in library modules
**Status:** PASS  
All 15 modules in `lib/sciagent/` use `return`, never `exit`:
- activate.sh, block.sh, claude_settings.sh, collisions.sh, deactivate.sh, eject.sh, frontmatter.sh, inject.sh, new.sh, roles.sh, skill_deps.sh, stack.sh, status.sh, symlinks.sh, validate.sh

**Evidence:** Test `test_no_exit_in_libs.sh` passes.

#### ✅ All VERB_MODULES present and accounted for
**Status:** PASS

Declared in `bin/sciagent` (lines 82–91): 15 modules  
Found in `lib/sciagent/`: 15 modules  
**Missing modules:** 0  
**Unreferenced modules:** 0

---

### C.2: Orphaned References (5 findings)

| Reference | File:Line | Status | Category |
|-----------|-----------|--------|----------|
| `scrna-atlas` (deleted role) | README.md:80 | ⚠️ ORPHAN | README references non-existent role |
| `planning` (deleted role) | README.md:80 | ⚠️ ORPHAN | README references non-existent role |
| `templates/project/software-tool/` (mismatch) | lib/sciagent/new.sh:170 | ⚠️ ORPHAN | Code expects dir, actual name is `software/` |

**Details:**

#### README.md Line 80 (ORPHAN #1 & #2)
```markdown
Toggle whichever combination fits the session — `base` + `pathway-signature` 
for downstream interpretation, `scrna-atlas` + `planning` to layer a research 
stack on top of atlas work, `architect` solo for design sessions. 
`sciagent list roles` enumerates what's available.
```

**Problem:** References `scrna-atlas` and `planning`, both deleted per CHANGELOG.md v1.7.0 (2026-06-02):
```
- Deleted 6 redundant roles that bloated the UI and duplicated coverage: 
  `planning`, `annotator`, `scrna-atlas`, `multiome-analysis`, `multiome-grn`, 
  `dc-dictionary`. The clean framework is now 4 roles: `base` (general scRNA 
  foundation), `scatac-regulatory` (chromatin/ATAC overlay), `pathway-signature` 
  (pathway/functional overlay), `architect` (software architecture, standalone).
```

**Current Roles:** architect, base, pathway-signature, scatac-regulatory, software-tool

**Severity:** [NIT] — Documentation stale, but does not break functionality.

#### lib/sciagent/new.sh Line 170 (ORPHAN #3)
```bash
_render_tree "$proj_tpl/$type"   "$dir" "$force" || return 1
```

With `type="software-tool"`, this expands to:
```
/path/to/templates/project/software-tool
```

**But actual directory is:**
```
/path/to/templates/project/software
```

**Why It Works:** `_render_tree()` (line 57) has a guard:
```bash
[[ -d "$tpl_root" ]] || return 0
```

So when the directory doesn't exist, it silently succeeds, rendering zero files. This means `software-tool` projects are created with skeleton directories but **NO type-specific templates** (AGENTS.md.template, README.md.template, tool_config.yaml.template).

**Test Evidence:** `test_new_project_types.sh` passes, meaning the project scaffolds (directories created), but we didn't verify that type-specific templates were rendered.

**Manual Test Output:**
```
$ sciagent new project --type software-tool test
wrote: test/docs/_internal/README.md
wrote: test/CLAUDE.md
wrote: test/.gitignore
Project scaffolded at test (type: software-tool).
```

The rendered files come from `_common/` (CLAUDE.md, .gitignore, docs/_internal/README.md), NOT from `software/` (which should provide tool_config.yaml.template, AGENTS.md.template, etc.).

**Severity:** [BUG] — Silent failure: `software-tool` projects don't receive type-specific templates.

---

### C.3: Skill Metadata Completeness

#### Total Skills: 70 (excluding _TEMPLATE)
| Field | Present | Missing | Pct |
|-------|---------|---------|-----|
| scope | 70 | 0 | 100% |
| last-reviewed | 68 | 2 | 97% |

#### Missing last-reviewed (2 skills)
| Skill Name | Metadata Status |
|-----------|-----------------|
| architecture-first-dev | No `last-reviewed` entry |
| skill-creator | No `last-reviewed` entry |

**Impact:** These skills do not benefit from the legacy exemption in `test_skill_scope_lint.sh`. If they ever exceed the body-line cap **and** have a scope vocabulary at or after the 2026-05-24 cutoff, they will fail hard. Currently not at risk (modern scopes, reasonable body sizes).

**Severity:** [NIT] — Non-blocking metadata hygiene issue.

---

### C.4: Deprecated Content (4 files, 0 active wiring)

**Location:** `deprecated/`

| File | Purpose | Still Referenced |
|------|---------|-------------------|
| install_claude.sh | Legacy Claude CLI installer | ❌ No |
| install_codex.sh | Legacy Codex CLI installer | ✅ Yes (CHANGELOG only) |
| install_gemini.sh | Legacy Gemini CLI installer | ✅ Yes (CHANGELOG only) |
| README.md | Deprecated directory marker | ❌ No |

**Status:** ✅ CLEAN  
Deprecated content is appropriately isolated. All active references are in CHANGELOG.md (historical record), not in active code.

---

### C.5: Scope Violations (ADR-003 regression block)

**Scope Vocabulary (post-cutoff 2026-05-24):**
- `concept` — ≤500 body lines
- `implementation` — ≤350 body lines

**Legacy Vocabulary (rejected by test):**
- `atomic` — ❌ **FORBIDDEN** (exactly 1 violation found)
- `orchestrator` — ❌ **FORBIDDEN**
- `foundation` — ❌ **FORBIDDEN**

| Skill | Scope | Last-Reviewed | Status |
|-------|-------|---------------|--------|
| mllmcelltype-consensus-annotation | **atomic** | 2026-05-29 | ❌ **FAIL** (legacy vocab + post-cutoff date) |

---

## Part D: Findings & Recommendations

### [BLOCKER] #1: test_skill_scope_lint.sh Failure

**Issue:** Skill `mllmcelltype-consensus-annotation` uses deprecated scope vocabulary `atomic`.

**Files Affected:**
- `skills/mllmcelltype-consensus-annotation/SKILL.md:5`

**Remediation:**
```diff
  metadata:
-   scope: atomic
+   scope: concept
```

Or:
```diff
  metadata:
-   scope: atomic
+   scope: implementation
```

**Timeline:** **CRITICAL** — blocks CI/test suite. Must be fixed before next commit.

---

### [ORPHAN] #2: README.md references deleted roles

**Issue:** README.md line 80 references `scrna-atlas` and `planning`, both deleted in v1.7.0.

**Files Affected:**
- `README.md:80`

**Current (Stale):**
```
Toggle whichever combination fits the session — `base` + `pathway-signature` 
for downstream interpretation, `scrna-atlas` + `planning` to layer a research 
stack on top of atlas work, `architect` solo for design sessions.
```

**Suggested Fix:**
```
Toggle whichever combination fits the session — `base` + `pathway-signature` 
for downstream interpretation, `base` + `scatac-regulatory` for chromatin 
analysis, or `architect` solo for design sessions. Run `sciagent list roles` 
to see all available role combinations.
```

**Timeline:** LOW — documentation only, no functional impact.

---

### [BUG] #3: Template directory mismatch (software vs software-tool)

**Issue:** Code in `lib/sciagent/new.sh` expects `templates/project/software-tool/` but the actual directory is `templates/project/software/`.

**Files Affected:**
- `lib/sciagent/new.sh:170`  
- Template directory naming: `templates/project/software/` (actual) vs `software-tool/` (expected by code)

**Impact:** Projects created with `sciagent new project --type software-tool` silently skip type-specific templates because the directory doesn't exist. The guard clause `[[ -d "$tpl_root" ]] || return 0` in `_render_tree()` makes this silent (returns success without rendering).

**Evidence:**
```bash
$ sciagent new project --type software-tool test
# Renders files from _common/ (CLAUDE.md, .gitignore, docs/_internal/README.md)
# BUT SKIPS software/ (tool_config.yaml.template, AGENTS.md.template, README.md.template)
```

**Why Tests Pass:** `test_new_project_types.sh` validates that directories are created (✅), not that templates are rendered.

**Remediation:**
Either (A) rename the directory:
```bash
mv templates/project/software templates/project/software-tool
```

Or (B) update code to match actual directory name:
```diff
-    _render_tree "$proj_tpl/$type"   "$dir" "$force" || return 1
+    local type_dir="$type"
+    [[ "$type" == "software-tool" ]] && type_dir="software"
+    _render_tree "$proj_tpl/$type_dir"   "$dir" "$force" || return 1
```

**Recommendation:** Use (A) rename — semantic consistency is clearer. Update code comments if needed to explain the historical name.

**Timeline:** MEDIUM — affects new project scaffolding experience but doesn't break existing workflows.

---

### [NIT] #4: Missing last-reviewed metadata

**Issue:** 2 skills lack `metadata.last-reviewed` entries.

**Files Affected:**
- `skills/architecture-first-dev/SKILL.md`
- `skills/skill-creator/SKILL.md`

**Impact:** Non-blocking. These skills are not subject to the body-line cap if reviewed before 2026-05-24, but their absence means they'll be hard-failed if caps are ever exceeded after that date.

**Remediation:** Add entries:
```yaml
metadata:
  last-reviewed: 2026-06-02
```

**Timeline:** LOW — optional hygiene improvement.

---

## Part E: Validation Checklist

| Check | Status | Notes |
|-------|--------|-------|
| Test suite runs | ✅ PASS | 67/68 tests pass |
| Scope-lint regression block | ❌ FAIL | 1 skill with `scope: atomic` |
| sciagent validate (quiet & non-quiet) | ✅ PASS | All checks pass |
| exit() in libs | ✅ PASS | No violations |
| Missing modules | ✅ PASS | All declared modules present |
| Dead modules | ✅ PASS | No unreferenced modules |
| Deprecated content isolated | ✅ PASS | No active wiring |
| Skill metadata complete | ⚠️ PARTIAL | 2/70 missing last-reviewed |
| Cross-namespace collisions | ✅ PASS | 2 allowlisted (architect, architecture-treemap) |

---

## Audit Output Artifact

**Full test run output:**
```
[command output truncated for brevity — see run-all.sh output]
passed: 68
failed: 1
failures:
  - test_skill_scope_lint.sh
```

**sciagent validate output:**
```
validate: warning — name 'architect' appears as agent, command, and role (mounting both is supported; ensure the overlap is intentional)
validate: warning — name 'architecture-treemap' appears as both skill and command (mounting both is supported; ensure the overlap is intentional)
sciagent validate: all checks passed
```

---

## Summary by Severity

### 🔴 Critical / Blocker (Fix Before Release)
1. **[BLOCKER]** test_skill_scope_lint.sh fails: `mllmcelltype-consensus-annotation` uses deprecated scope `atomic`
   - File: `skills/mllmcelltype-consensus-annotation/SKILL.md:5`
   - Fix: Change `scope: atomic` → `scope: concept` (or `implementation`)

### 🟠 High (Fix Soon)
1. **[BUG]** Template directory mismatch: code expects `templates/project/software-tool/` but directory is `software/`
   - Files: `lib/sciagent/new.sh:170`, `templates/project/`
   - Impact: `software-tool` projects skip type-specific templates silently
   - Fix: Rename directory or update code reference

### 🟡 Medium (Document & Schedule)
1. **[ORPHAN]** README.md references deleted roles `scrna-atlas` and `planning`
   - File: `README.md:80`
   - Fix: Update example to use current roles (base, pathway-signature, scatac-regulatory, architect)

### 🔵 Low (Nice to Have)
1. **[NIT]** 2 skills missing `last-reviewed` metadata (architecture-first-dev, skill-creator)
   - Files: `skills/{architecture-first-dev,skill-creator}/SKILL.md`
   - Fix: Add `last-reviewed: 2026-06-02` to metadata

---

## Integrity Audit Conclusion

**Overall Health:** 🟡 CONDITIONAL PASS  
- Test suite **functionally robust** (67/68 passing)
- **One known blocker** in scope vocabulary (regression-block test, not a logic bug)
- **Silent template mismatch** masks a scaffolding bug
- **Stale documentation** references deleted roles

**Path to Green:**
1. Fix `mllmcelltype-consensus-annotation` scope (5 min)
2. Rename/update templates directory reference (5 min)
3. Update README.md roles (5 min)
4. Add last-reviewed metadata (optional, 2 min)

**Estimated remediation time:** 15 minutes

