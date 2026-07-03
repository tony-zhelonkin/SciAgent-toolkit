#!/usr/bin/env bash
# tests/test_user_settings_seed.sh — Group-A user-scope seeders in
# claude_settings.sh: claude_settings_ensure_user_defaults + _statusline.
#
# Covers the tier-0 acceptance checks:
#   - absent  → creates ${CLAUDE_CONFIG_DIR:-$HOME/.claude}/settings.json with the
#               gold defaults + an executable statusline.sh.
#   - non-clobber → a user-set value (effortLevel=low) and a custom key both
#               survive a re-run, while still backfilling missing template keys.
#   - idempotent → a second run on an already-canonical file is byte-identical.
#
# Hermetic: scratch HOME/CLAUDE_CONFIG_DIR under the tmpdir; nothing touches the
# real user env.
set -u
. "$(dirname "$0")/_lib.sh"
. "$TOOLKIT_ROOT/lib/sciagent/claude_settings.sh"

command -v jq >/dev/null 2>&1 || { echo "SKIP [$_TEST_NAME] jq unavailable"; exit 0; }

setup_tmpdir
export HOME="$TMPDIR_TEST/home"
mkdir -p "$HOME"

# ----- Case 1: absent → create defaults + executable statusline ------------
export CLAUDE_CONFIG_DIR="$TMPDIR_TEST/cfg1"
claude_settings_ensure_user_defaults   >/dev/null
claude_settings_ensure_user_statusline >/dev/null

settings="$CLAUDE_CONFIG_DIR/settings.json"
statusline="$CLAUDE_CONFIG_DIR/statusline.sh"
assert_file_exists "$settings"   "case1: settings.json created"
assert_file_exists "$statusline" "case1: statusline.sh created"
[[ -x "$statusline" ]] || { echo "FAIL [$_TEST_NAME] case1: statusline.sh not executable" >&2; exit 1; }

# Gold defaults present.
assert_eq "$(jq -r '.editorMode' "$settings")"            "vim"   "case1: editorMode"
assert_eq "$(jq -r '.effortLevel' "$settings")"           "xhigh" "case1: effortLevel"
assert_eq "$(jq -r '.alwaysThinkingEnabled' "$settings")" "true"  "case1: alwaysThinkingEnabled"
assert_eq "$(jq -r '.autoMemoryEnabled' "$settings")"     "false" "case1: autoMemoryEnabled"
# User-level statusLine path must be user-scoped (~/.claude/...), not $CLAUDE_PROJECT_DIR.
assert_eq "$(jq -r '.statusLine.command' "$settings")" "~/.claude/statusline.sh" "case1: statusLine.command user-scoped"

# ----- Case 2: non-clobber — user value + custom key survive, gaps filled ---
export CLAUDE_CONFIG_DIR="$TMPDIR_TEST/cfg2"
mkdir -p "$CLAUDE_CONFIG_DIR"
settings2="$CLAUDE_CONFIG_DIR/settings.json"
# A partial file: one overridden template key + one wholly-custom key, and
# deliberately missing editorMode / statusLine so backfill has work to do.
cat > "$settings2" <<'EOF'
{
  "effortLevel": "low",
  "myCustomKey": "keepme"
}
EOF

claude_settings_ensure_user_defaults >/dev/null

# User-set + custom keys preserved (existing-wins reverse merge).
assert_eq "$(jq -r '.effortLevel' "$settings2")" "low"    "case2: user effortLevel preserved"
assert_eq "$(jq -r '.myCustomKey' "$settings2")" "keepme" "case2: custom key preserved"
# Missing template keys backfilled.
assert_eq "$(jq -r '.editorMode' "$settings2")"          "vim"                     "case2: editorMode backfilled"
assert_eq "$(jq -r '.statusLine.command' "$settings2")"  "~/.claude/statusline.sh" "case2: statusLine backfilled"

# ----- Case 3: idempotent — steady-state re-run is byte-identical ----------
export CLAUDE_CONFIG_DIR="$TMPDIR_TEST/cfg3"
# Reach jq-canonical steady state (first run copies the raw template verbatim,
# the second run canonicalises it via the merge).
claude_settings_ensure_user_defaults >/dev/null
claude_settings_ensure_user_defaults >/dev/null
settings3="$CLAUDE_CONFIG_DIR/settings.json"
cp "$settings3" "$TMPDIR_TEST/cfg3.orig"

out=$(claude_settings_ensure_user_defaults)
assert_file_eq "$settings3" "$TMPDIR_TEST/cfg3.orig" "case3: steady-state re-run byte-identical"
if printf '%s\n' "$out" | grep -qE 'wrote:|updated:'; then
    echo "FAIL [$_TEST_NAME] case3: steady-state re-run reported a change" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

pass
