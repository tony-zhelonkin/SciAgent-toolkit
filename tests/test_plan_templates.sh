#!/usr/bin/env bash
# tests/test_plan_templates.sh
# Verify that the plan templates exist, carry the required section headers,
# and render cleanly via _subst (no broken {{TOKEN}} after DATE substitution).

set -u
# shellcheck source=/dev/null
. "$(dirname "$0")/_lib.sh"

PLAN_TPL_DIR="$TOOLKIT_ROOT/templates/plan"
INDEX_TPL="$PLAN_TPL_DIR/00_INDEX.md.template"
PHASE_TPL="$PLAN_TPL_DIR/NN_slug.md.template"

# ---------------------------------------------------------------------------
# 1. Files exist
# ---------------------------------------------------------------------------
assert_file_exists "$INDEX_TPL" "00_INDEX.md.template missing"
assert_file_exists "$PHASE_TPL" "NN_slug.md.template missing"

# ---------------------------------------------------------------------------
# 2. INDEX template — required sections and structural elements
# ---------------------------------------------------------------------------
assert_grep "seq | slug | title"              "$INDEX_TPL" "INDEX missing phase table header"
assert_grep "## Phase table"                  "$INDEX_TPL" "INDEX missing ## Phase table"
assert_grep "tier"                            "$INDEX_TPL" "INDEX phase table missing tier column"
assert_grep "concern"                         "$INDEX_TPL" "INDEX phase table missing concern column"
assert_grep "depends_on"                      "$INDEX_TPL" "INDEX phase table missing depends_on column"
assert_grep "Opus REVIEW CHECKPOINT"          "$INDEX_TPL" "INDEX missing Opus review checkpoint row example"
assert_grep "## toolkit_additions"            "$INDEX_TPL" "INDEX missing ## toolkit_additions"
assert_grep "## config_additions"             "$INDEX_TPL" "INDEX missing ## config_additions"
assert_grep "## Global notes"                 "$INDEX_TPL" "INDEX missing ## Global notes"
assert_grep "load_or_compute"                 "$INDEX_TPL" "INDEX missing idempotency note (load_or_compute)"
assert_grep "compute.*viz\|viz.*compute"      "$INDEX_TPL" "INDEX missing compute→viz split note"
assert_grep "claim ladder\|Claim ladder"      "$INDEX_TPL" "INDEX missing claim ladder note"
assert_grep "{{DATE}}"                        "$INDEX_TPL" "INDEX missing {{DATE}} token"

# ---------------------------------------------------------------------------
# 3. Phase brief template — required fixed sections
# ---------------------------------------------------------------------------
assert_grep "## 1. Objective"                 "$PHASE_TPL" "phase brief missing ## 1. Objective"
assert_grep "## 2. Scope"                     "$PHASE_TPL" "phase brief missing ## 2. Scope"
assert_grep "## 3. Inputs"                    "$PHASE_TPL" "phase brief missing ## 3. Inputs"
assert_grep "## 4. Outputs"                   "$PHASE_TPL" "phase brief missing ## 4. Outputs"
assert_grep "## 5. Implementation"            "$PHASE_TPL" "phase brief missing ## 5. Implementation"
assert_grep "## 6. Captions"                  "$PHASE_TPL" "phase brief missing ## 6. Captions"
assert_grep "## 7. Acceptance checks"         "$PHASE_TPL" "phase brief missing ## 7. Acceptance checks"
assert_grep "## 8. Gotchas"                   "$PHASE_TPL" "phase brief missing ## 8. Gotchas"
assert_grep "Figure-style contract"           "$PHASE_TPL" "phase brief missing Figure-style contract declaration"
assert_grep "COMPUTE-ONLY\|VIZ-ONLY\|MIXED"  "$PHASE_TPL" "phase brief missing VIZ/COMPUTE declaration"
assert_grep "Out of scope"                    "$PHASE_TPL" "phase brief missing Out of scope with grep"
assert_grep "save_overview"                   "$PHASE_TPL" "phase brief missing save_overview reference"
assert_grep "{{DATE}}"                        "$PHASE_TPL" "phase brief missing {{DATE}} token"

# ---------------------------------------------------------------------------
# 4. _subst render — no unresolved {{TOKEN}} after DATE substitution
# ---------------------------------------------------------------------------
setup_tmpdir

rendered_index="$TMPDIR_TEST/00_INDEX.md"
rendered_phase="$TMPDIR_TEST/NN_slug.md"

# Simulate _subst: substitute {{DATE}} only (the only auto-fill token in plan templates)
sed -e "s|{{DATE}}|2099-01-01|g" "$INDEX_TPL" > "$rendered_index"
sed -e "s|{{DATE}}|2099-01-01|g" "$PHASE_TPL" > "$rendered_phase"

# Confirm substitution landed
assert_grep "2099-01-01"  "$rendered_index" "INDEX: {{DATE}} was not substituted"
assert_grep "2099-01-01"  "$rendered_phase" "phase: {{DATE}} was not substituted"

# No remaining {{TOKEN}} (only intentional <angle-bracket> placeholders should remain)
if grep -qE '\{\{[A-Z_]+\}\}' "$rendered_index"; then
    echo "FAIL [$_TEST_NAME] INDEX has unresolved {{TOKEN}} after _subst:" >&2
    grep -nE '\{\{[A-Z_]+\}\}' "$rendered_index" >&2
    exit 1
fi
if grep -qE '\{\{[A-Z_]+\}\}' "$rendered_phase"; then
    echo "FAIL [$_TEST_NAME] phase brief has unresolved {{TOKEN}} after _subst:" >&2
    grep -nE '\{\{[A-Z_]+\}\}' "$rendered_phase" >&2
    exit 1
fi

pass
