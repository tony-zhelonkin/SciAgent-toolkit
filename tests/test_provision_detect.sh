#!/usr/bin/env bash
# tests/test_provision_detect.sh — harness.sh detection API.
#
# Covers the tier-1 acceptance check "harness_detect() correctly reports
# installed vs absent for all five":
#   - scratch HOME with fake ~/.pi/agent + ~/.claude, harness bins off PATH →
#     harness_detect lists exactly `claude` then `pi`, canonical order, one each.
#   - empty CLAUDE_CONFIG_DIR= / CODEX_HOME= fall back to the $HOME defaults.
#   - harness_is_present returns 0/1/2 for present/absent/unknown.
#
# Hermetic: scratch HOME under the tmpdir, PATH pruned so none of the five real
# harness binaries resolve.
set -u
. "$(dirname "$0")/_lib.sh"

# Prune PATH to standard system dirs so claude/pi/agy etc. (installed under the
# real user's home) are invisible; keep coreutils/jq available for _lib.sh.
export PATH=/usr/bin:/bin
. "$TOOLKIT_ROOT/lib/sciagent/harness.sh"

setup_tmpdir
export HOME="$TMPDIR_TEST/home"
unset CLAUDE_CONFIG_DIR CODEX_HOME
mkdir -p "$HOME/.pi/agent" "$HOME/.claude"

# ----- harness_detect: exactly claude then pi, canonical order -------------
detected=$(harness_detect)
assert_eq "$detected" $'claude\npi' "detect lists exactly claude,pi in canonical order"

# codex / agy / opencode must be absent from the listing.
for absent in codex agy opencode; do
    if printf '%s\n' "$detected" | grep -qx "$absent"; then
        echo "FAIL [$_TEST_NAME] $absent should not be detected" >&2
        printf '%s\n' "$detected" >&2
        exit 1
    fi
done

# One line each (no duplicates).
assert_eq "$(printf '%s\n' "$detected" | grep -cx claude)" "1" "claude listed once"
assert_eq "$(printf '%s\n' "$detected" | grep -cx pi)"     "1" "pi listed once"

# ----- empty CLAUDE_CONFIG_DIR / CODEX_HOME fall back to $HOME defaults -----
export CLAUDE_CONFIG_DIR="" CODEX_HOME=""
assert_eq "$(harness_config_dir claude)" "$HOME/.claude"       "empty CLAUDE_CONFIG_DIR falls back"
assert_eq "$(harness_config_dir codex)"  "$HOME/.codex"        "empty CODEX_HOME falls back"
assert_eq "$(harness_global_context_path claude)" "$HOME/.claude/CLAUDE.md" "empty CLAUDE_CONFIG_DIR context path falls back"
# Detection still resolves without error under the empty-var fallback.
assert_eq "$(harness_detect)" $'claude\npi' "detect unchanged under empty env fallback"
unset CLAUDE_CONFIG_DIR CODEX_HOME

# ----- harness_is_present: present / absent / unknown exit codes -----------
assert_exit 0 harness_is_present claude     # config dir exists
assert_exit 0 harness_is_present pi         # config dir exists
assert_exit 1 harness_is_present codex      # bin off PATH, no config dir
assert_exit 1 harness_is_present opencode   # bin off PATH, no config dir
assert_exit 2 harness_is_present bogus      # unknown harness name

pass
