#!/usr/bin/env bash
# tests/test_enforcement_hooks.sh — P14 enforcement-hooks test suite.
#
# Tests:
#   1. Scaffolding: sciagent new project creates .claude/{settings.json,hooks/*.sh}
#   2. no_ephemeral — python3 -c + 03_results write: default exit 0 + advisory;
#                     SCIAGENT_STRICT=1: exit 2.
#   3. no_ephemeral — _scratch/ target: exit 0, no advisory (sanctioned).
#   4. no_ephemeral — committed-script run (Rscript 02_analysis/scripts/…): exit 0.
#   5. caption_sweep — figure present, no caption: exit 0, reminder on stderr.
#   6. caption_sweep — conformant caption present: exit 0, no reminder.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir

# ---------------------------------------------------------------------------
# Test 1: scaffolding.
# ---------------------------------------------------------------------------
PROJ="$TMPDIR_TEST/P"
"$SCIAGENT_TOOLKIT/bin/sciagent" new project "$PROJ" --type analysis >/dev/null 2>&1 \
    || { echo "FAIL [$_TEST_NAME] test1: scaffold failed" >&2; exit 1; }

assert_file_exists "$PROJ/.claude/hooks/no_ephemeral.sh"   "no_ephemeral.sh not scaffolded"
assert_file_exists "$PROJ/.claude/hooks/caption_sweep.sh"  "caption_sweep.sh not scaffolded"
assert_file_exists "$PROJ/.claude/settings.json"           "settings.json not scaffolded"

# settings.json must register PreToolUse and Stop hooks.
assert_grep 'PreToolUse'      "$PROJ/.claude/settings.json" "settings.json missing PreToolUse"
assert_grep 'Stop'            "$PROJ/.claude/settings.json" "settings.json missing Stop"
assert_grep 'no_ephemeral'    "$PROJ/.claude/settings.json" "settings.json missing no_ephemeral"
assert_grep 'caption_sweep'   "$PROJ/.claude/settings.json" "settings.json missing caption_sweep"

# ---------------------------------------------------------------------------
# Helper: build a minimal PreToolUse JSON payload.
# ---------------------------------------------------------------------------
_make_pretool_json() {
    local cmd="$1"
    # Use python3 for reliable JSON quoting.
    python3 -c "
import json, sys
payload = {
    'tool_name': 'Bash',
    'tool_input': {'command': sys.argv[1]},
    'cwd': '$PROJ'
}
print(json.dumps(payload))
" "$cmd"
}

# ---------------------------------------------------------------------------
# Helper: build a minimal Stop JSON payload.
# ---------------------------------------------------------------------------
_make_stop_json() {
    printf '{"stop_reason":"end_turn","session_id":"test"}\n'
}

HOOK_EPHEMERAL="$PROJ/.claude/hooks/no_ephemeral.sh"
HOOK_CAPTION="$PROJ/.claude/hooks/caption_sweep.sh"

# ---------------------------------------------------------------------------
# Test 2a: python3 -c with 03_results write → default: exit 0 with advisory.
# ---------------------------------------------------------------------------
CMD_EPHEMERAL="python3 -c \"import matplotlib; fig=matplotlib.pyplot.figure(); fig.savefig('03_results/01_qc/figures/probe.png')\""
PAYLOAD=$(_make_pretool_json "$CMD_EPHEMERAL")

set +e
out2a=$(printf '%s\n' "$PAYLOAD" | bash "$HOOK_EPHEMERAL" 2>&1)
rc2a=$?
set -e

if [[ "$rc2a" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] test2a: default mode expected exit 0, got $rc2a" >&2
    printf '%s\n' "$out2a" >&2
    exit 1
fi

if ! printf '%s\n' "$out2a" | grep -qi 'reproducible\|02_analysis/scripts\|ephemeral\|no-ephemeral'; then
    echo "FAIL [$_TEST_NAME] test2a: expected advisory message in stderr output" >&2
    printf '%s\n' "$out2a" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 2b: same command + SCIAGENT_STRICT=1 → exit 2 (BLOCK).
# ---------------------------------------------------------------------------
set +e
out2b=$(printf '%s\n' "$PAYLOAD" | SCIAGENT_STRICT=1 bash "$HOOK_EPHEMERAL" 2>&1)
rc2b=$?
set -e

if [[ "$rc2b" -ne 2 ]]; then
    echo "FAIL [$_TEST_NAME] test2b: SCIAGENT_STRICT=1 expected exit 2, got $rc2b" >&2
    printf '%s\n' "$out2b" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 3: write under _scratch/ → exit 0, NO advisory (sanctioned).
# ---------------------------------------------------------------------------
CMD_SCRATCH="python3 -c \"open('03_results/_scratch/probe.csv','w').write('x')\""
PAYLOAD_SCRATCH=$(_make_pretool_json "$CMD_SCRATCH")

set +e
out3=$(printf '%s\n' "$PAYLOAD_SCRATCH" | bash "$HOOK_EPHEMERAL" 2>&1)
rc3=$?
set -e

if [[ "$rc3" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] test3: _scratch/ sanctioned, expected exit 0, got $rc3" >&2
    printf '%s\n' "$out3" >&2
    exit 1
fi

if printf '%s\n' "$out3" | grep -qi 'reproducible\|02_analysis/scripts\|ephemeral\|no-ephemeral'; then
    echo "FAIL [$_TEST_NAME] test3: _scratch/ target must NOT produce advisory" >&2
    printf '%s\n' "$out3" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 4: normal committed-script run → exit 0, no block even under strict.
# ---------------------------------------------------------------------------
CMD_COMMITTED="Rscript 02_analysis/scripts/10_qc_viz.R"
PAYLOAD_COMMITTED=$(_make_pretool_json "$CMD_COMMITTED")

set +e
out4=$(printf '%s\n' "$PAYLOAD_COMMITTED" | SCIAGENT_STRICT=1 bash "$HOOK_EPHEMERAL" 2>&1)
rc4=$?
set -e

if [[ "$rc4" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] test4: committed script must exit 0 even under STRICT, got $rc4" >&2
    printf '%s\n' "$out4" >&2
    exit 1
fi

if printf '%s\n' "$out4" | grep -qi 'no-ephemeral.*convention'; then
    echo "FAIL [$_TEST_NAME] test4: committed script must not produce ephemeral advisory" >&2
    printf '%s\n' "$out4" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 5: caption_sweep — figure present, no caption → exit 0, reminder.
# ---------------------------------------------------------------------------
FIG_PROJ="$TMPDIR_TEST/figproj"
mkdir -p "$FIG_PROJ/03_results/01_qc/figures/_overview"
touch "$FIG_PROJ/03_results/01_qc/figures/_overview/x.screen.png"
# No README.md → no caption.

set +e
out5=$(printf '%s\n' "$(_make_stop_json)" \
    | CLAUDE_PROJECT_DIR="$FIG_PROJ" bash "$HOOK_CAPTION" 2>&1)
rc5=$?
set -e

if [[ "$rc5" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] test5: uncaptioned default expected exit 0, got $rc5" >&2
    printf '%s\n' "$out5" >&2
    exit 1
fi

if ! printf '%s\n' "$out5" | grep -qi 'caption\|uncaptioned\|caption-sweep'; then
    echo "FAIL [$_TEST_NAME] test5: expected caption reminder in output" >&2
    printf '%s\n' "$out5" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 6: caption_sweep — conformant caption present → exit 0, no reminder.
# ---------------------------------------------------------------------------
CAP_PROJ="$TMPDIR_TEST/capproj"
mkdir -p "$CAP_PROJ/03_results/01_qc/figures/_overview"
touch "$CAP_PROJ/03_results/01_qc/figures/_overview/x.screen.png"
# Conformant README.md with path-qualified heading.
cat > "$CAP_PROJ/03_results/01_qc/README.md" <<'MD'
# 01_qc

## figures/_overview/x.screen.png

A figure.
MD

set +e
out6=$(printf '%s\n' "$(_make_stop_json)" \
    | CLAUDE_PROJECT_DIR="$CAP_PROJ" bash "$HOOK_CAPTION" 2>&1)
rc6=$?
set -e

if [[ "$rc6" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] test6: conformant caption expected exit 0, got $rc6" >&2
    printf '%s\n' "$out6" >&2
    exit 1
fi

if printf '%s\n' "$out6" | grep -qi 'uncaptioned\|no caption\|caption-sweep'; then
    echo "FAIL [$_TEST_NAME] test6: conformant project must NOT produce caption reminder" >&2
    printf '%s\n' "$out6" >&2
    exit 1
fi

pass
