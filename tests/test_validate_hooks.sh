#!/usr/bin/env bash
# tests/test_validate_hooks.sh — opt-in `--check hooks` guardrail.
#
# Regression cover for a state that shipped: settings.json REGISTERS
# PreToolUse/Stop hooks by path, but the bodies were rendered only by
# `sciagent new project`, so any project retrofitted by `activate` got hook
# registration with no hook files. Meta-Aging/14616-DM ran that way with no
# .claude/hooks/ directory at all, silently disabling the (c) GUARDRAIL layer.
#
# Tests:
#   1. Registered + present (invoked via `bash <path>`, mode 0644) → clean.
#      Mode must NOT matter here: demanding +x would warn on every correctly
#      provisioned repo.
#   2. Registered + ABSENT → WARN on stderr, exit 0 by default; exit 1 --strict.
#   3. Registered for DIRECT execution + not executable → WARN.
#   4. Same, after chmod +x → clean.
#   5. No .claude/settings.json → clean no-op.
#   6. settings.json with no `hooks` key → clean no-op.
#   7. claude_settings_ensure_hooks materializes the bodies and is idempotent.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

fail() {
    echo "FAIL [$_TEST_NAME] $1" >&2
    shift
    [[ $# -gt 0 ]] && { echo "--- output ---" >&2; printf '%s\n' "$@" >&2; }
    exit 1
}

# write_settings <projdir> <command>
# The command routinely contains double quotes (`bash "$CLAUDE_PROJECT_DIR/..."`),
# so escape them for JSON — an unescaped quote makes the file unparseable and jq
# then yields no hooks at all, which would make this test vacuously pass.
write_settings() {
    local esc=${2//\"/\\\"}
    mkdir -p "$1/.claude"
    cat > "$1/.claude/settings.json" <<EOF
{
  "hooks": {
    "Stop": [
      { "matcher": "*", "hooks": [ { "type": "command", "command": "$esc" } ] }
    ]
  }
}
EOF
    # Guard against the failure mode above: the fixture must be valid JSON.
    if command -v jq >/dev/null 2>&1; then
        jq -e . "$1/.claude/settings.json" >/dev/null 2>&1 \
            || fail "test fixture is not valid JSON: $1/.claude/settings.json"
    fi
}

run_check() {   # run_check <projdir> [--strict] -> sets OUT / RC
    set +e
    OUT=$("$SCIAGENT" validate --check hooks --project-dir "$1" ${2:-} 2>&1)
    RC=$?
    set -e
}

# --- 1. registered + present, interpreter-invoked, mode 0644 → clean ---------
P1="$TMPDIR_TEST/p_present"
write_settings "$P1" 'bash "$CLAUDE_PROJECT_DIR/.claude/hooks/no_ephemeral.sh"'
mkdir -p "$P1/.claude/hooks"
echo '#!/usr/bin/env bash' > "$P1/.claude/hooks/no_ephemeral.sh"
chmod 644 "$P1/.claude/hooks/no_ephemeral.sh"
run_check "$P1"
[[ "$RC" -eq 0 ]] || fail "present-hook project exited $RC (expected 0)" "$OUT"
grep -q 'WARN hooks' <<<"$OUT" && fail "warned on a present, bash-invoked hook (mode 0644 is fine)" "$OUT"

# --- 2. registered + absent → WARN, soft by default, hard under --strict -----
P2="$TMPDIR_TEST/p_absent"
write_settings "$P2" 'bash "$CLAUDE_PROJECT_DIR/.claude/hooks/caption_sweep.sh"'
run_check "$P2"
[[ "$RC" -eq 0 ]] || fail "absent-hook project exited $RC by default (expected soft 0)" "$OUT"
grep -q 'caption_sweep.sh is registered' <<<"$OUT" \
    || fail "no WARN for a registered-but-absent hook" "$OUT"
run_check "$P2" --strict
[[ "$RC" -eq 1 ]] || fail "absent-hook project exited $RC under --strict (expected 1)" "$OUT"

# --- 3/4. direct execution requires +x --------------------------------------
P3="$TMPDIR_TEST/p_direct"
write_settings "$P3" '$CLAUDE_PROJECT_DIR/.claude/hooks/direct.sh'
mkdir -p "$P3/.claude/hooks"
echo '#!/usr/bin/env bash' > "$P3/.claude/hooks/direct.sh"
chmod 644 "$P3/.claude/hooks/direct.sh"
run_check "$P3"
grep -q 'registered for direct execution but is not executable' <<<"$OUT" \
    || fail "no WARN for a directly-invoked non-executable hook" "$OUT"
chmod +x "$P3/.claude/hooks/direct.sh"
run_check "$P3"
grep -q 'WARN hooks' <<<"$OUT" && fail "still warning after chmod +x" "$OUT"

# --- 5. no settings.json → clean no-op --------------------------------------
P5="$TMPDIR_TEST/p_nosettings"; mkdir -p "$P5"
run_check "$P5"
[[ "$RC" -eq 0 ]] || fail "project without settings.json exited $RC (expected 0)" "$OUT"
grep -q 'WARN hooks' <<<"$OUT" && fail "warned on a project with no settings.json" "$OUT"

# --- 6. settings.json without a hooks key → clean no-op ---------------------
P6="$TMPDIR_TEST/p_nohooks"; mkdir -p "$P6/.claude"
echo '{"editorMode":"vim"}' > "$P6/.claude/settings.json"
run_check "$P6"
[[ "$RC" -eq 0 ]] || fail "settings.json without hooks exited $RC (expected 0)" "$OUT"
grep -q 'WARN hooks' <<<"$OUT" && fail "warned on settings.json carrying no hooks" "$OUT"

# --- 7. ensure_hooks materializes the bodies, idempotently ------------------
# Uses the REAL toolkit templates (the fake toolkit has no templates/ tree), so
# this also asserts the shipped settings template and hook bodies stay in sync.
# $0 is relative and setup_tmpdir changed cwd; _lib.sh already resolved this.
REAL_TK="$TOOLKIT_ROOT"
if [[ -d "$REAL_TK/templates/project/_common/.claude/hooks" ]]; then
    P7="$TMPDIR_TEST/p_materialize"; mkdir -p "$P7"
    first=$(cd "$P7" && SCIAGENT_TOOLKIT="$REAL_TK" bash -c '
        . "$SCIAGENT_TOOLKIT/lib/sciagent/claude_settings.sh"
        claude_settings_ensure_project_defaults >/dev/null
        claude_settings_ensure_hooks' 2>&1)
    grep -q 'wrote: .claude/hooks/' <<<"$first" \
        || fail "ensure_hooks wrote no hook bodies" "$first"

    # Every hook the shipped settings template registers must now exist.
    SCIAGENT_TOOLKIT="$REAL_TK" run_check "$P7"
    grep -q 'WARN hooks' <<<"$OUT" \
        && fail "hooks still missing after ensure_hooks materialized them" "$OUT"

    second=$(cd "$P7" && SCIAGENT_TOOLKIT="$REAL_TK" bash -c '
        . "$SCIAGENT_TOOLKIT/lib/sciagent/claude_settings.sh"
        claude_settings_ensure_hooks' 2>&1)
    grep -q 'wrote:' <<<"$second" \
        && fail "ensure_hooks re-wrote existing hooks (not idempotent)" "$second"
fi

pass
