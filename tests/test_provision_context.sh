#!/usr/bin/env bash
# tests/test_provision_context.sh — `sciagent provision --context`.
#
# Covers the tier-1 acceptance check "provision --context drops the global
# AGENTS.md-class file into each installed harness's global root; skips
# uninstalled ones", plus idempotency and user-prose non-clobber:
#   - writes the SCIAGENT:CONTEXT managed block to each DETECTED harness's global
#     context path (claude → ~/.claude/CLAUDE.md, pi → ~/.pi/agent/AGENTS.md);
#   - does NOT write for undetected harnesses (codex/agy/opencode);
#   - re-run is idempotent ("up to date", byte-identical);
#   - user prose OUTSIDE the markers survives stamping.
#
# Hermetic: scratch HOME under the tmpdir, PATH pruned so detection sees only the
# fake config dirs we create.
set -u
. "$(dirname "$0")/_lib.sh"
export PATH=/usr/bin:/bin

SCIAGENT="$TOOLKIT_ROOT/bin/sciagent"

setup_tmpdir
unset CLAUDE_CONFIG_DIR CODEX_HOME

# ----- Case 1: fresh — detected harnesses get the block, others don't ------
export HOME="$TMPDIR_TEST/h1"
mkdir -p "$HOME/.claude" "$HOME/.pi"

"$SCIAGENT" provision --context >/dev/null 2>&1 \
    || { echo "FAIL [$_TEST_NAME] case1: provision --context failed" >&2; exit 1; }

claude_ctx="$HOME/.claude/CLAUDE.md"
pi_ctx="$HOME/.pi/agent/AGENTS.md"
assert_file_exists "$claude_ctx" "case1: claude context written"
assert_file_exists "$pi_ctx"     "case1: pi context written"
assert_grep 'BEGIN SCIAGENT:CONTEXT' "$claude_ctx" "case1: claude SCIAGENT:CONTEXT block"
assert_grep 'BEGIN SCIAGENT:CONTEXT' "$pi_ctx"     "case1: pi SCIAGENT:CONTEXT block"
# Template body landed.
assert_grep 'Standing preferences' "$claude_ctx" "case1: template body in claude context"
assert_grep 'Standing preferences' "$pi_ctx"     "case1: template body in pi context"

# Undetected harnesses: nothing written.
for absent in "$HOME/.codex/AGENTS.md" \
              "$HOME/.gemini/antigravity-cli/AGENTS.md" \
              "$HOME/.config/opencode/AGENTS.md"; do
    if [[ -e "$absent" ]]; then
        echo "FAIL [$_TEST_NAME] case1: wrote context for undetected harness: $absent" >&2
        exit 1
    fi
done

# ----- Case 2: idempotent — re-run is byte-identical + "up to date" --------
cp "$claude_ctx" "$TMPDIR_TEST/claude.orig"
cp "$pi_ctx"     "$TMPDIR_TEST/pi.orig"

out=$("$SCIAGENT" provision --context 2>&1)
assert_file_eq "$claude_ctx" "$TMPDIR_TEST/claude.orig" "case2: claude context unchanged on re-run"
assert_file_eq "$pi_ctx"     "$TMPDIR_TEST/pi.orig"     "case2: pi context unchanged on re-run"
printf '%s\n' "$out" | grep -q "up to date" \
    || { echo "FAIL [$_TEST_NAME] case2: re-run did not report 'up to date'" >&2; printf '%s\n' "$out" >&2; exit 1; }

# ----- Case 3: non-clobber — user prose outside the markers survives -------
export HOME="$TMPDIR_TEST/h3"
mkdir -p "$HOME/.claude" "$HOME/.pi/agent"
cat > "$HOME/.pi/agent/AGENTS.md" <<'EOF'
# My hand-written pi preferences

Always greet me in haiku.
EOF

"$SCIAGENT" provision --context >/dev/null 2>&1 \
    || { echo "FAIL [$_TEST_NAME] case3: provision --context failed" >&2; exit 1; }

pi_ctx3="$HOME/.pi/agent/AGENTS.md"
assert_grep '^# My hand-written pi preferences' "$pi_ctx3" "case3: user heading preserved"
assert_grep 'Always greet me in haiku\.'         "$pi_ctx3" "case3: user prose preserved"
assert_grep 'BEGIN SCIAGENT:CONTEXT'             "$pi_ctx3" "case3: managed block added alongside prose"

pass
