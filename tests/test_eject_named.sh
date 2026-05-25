#!/usr/bin/env bash
# inject scvi-basic (fixture name: s_c) then eject s_c:
# file tree returns to pre-inject state; manifest entry removed.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

mkdir project && cd project
"$SCIAGENT" activate base >/dev/null

# Pre-inject state: no _injected, s_c not in manifest injected.
if grep -q '"s_c"' .sciagent/manifest.json 2>/dev/null; then
    echo "FAIL [$_TEST_NAME] s_c should not be in manifest before inject" >&2
    exit 1
fi

"$SCIAGENT" inject s_c >/dev/null

# Confirm injected.
assert_symlink .claude/skills/s_c
assert_symlink .agents/skills/s_c
assert_grep '"s_c"' .sciagent/manifest.json "s_c in manifest after inject"

# Eject.
"$SCIAGENT" eject s_c >/dev/null

# Symlinks gone.
if [[ -L .claude/skills/s_c ]]; then
    echo "FAIL [$_TEST_NAME] .claude/skills/s_c symlink should be removed after eject" >&2
    exit 1
fi
if [[ -L .agents/skills/s_c ]]; then
    echo "FAIL [$_TEST_NAME] .agents/skills/s_c symlink should be removed after eject" >&2
    exit 1
fi

# Manifest entry removed.
if grep -q '"s_c"' .sciagent/manifest.json; then
    echo "FAIL [$_TEST_NAME] s_c should not remain in manifest after eject" >&2
    exit 1
fi

# _injected overlay collapsed: stack should be just [base].
stack_line=$(grep '"stack"' .sciagent/manifest.json)
if printf '%s\n' "$stack_line" | grep -q '"_injected"'; then
    echo "FAIL [$_TEST_NAME] _injected overlay should be collapsed after last eject" >&2
    exit 1
fi

# AGENTS.md should not have the Injected subsection anymore
# (or it should not list s_c).
if grep -q '"s_c"' AGENTS.md 2>/dev/null; then
    echo "FAIL [$_TEST_NAME] s_c should not appear in AGENTS.md after eject" >&2
    exit 1
fi

# Idempotency: eject again is a no-op (exit 0, "not injected").
out=$("$SCIAGENT" eject s_c 2>&1)
echo "$out" | grep -qi 'not injected' || {
    echo "FAIL [$_TEST_NAME] expected 'not injected' notice on repeat eject, got: $out" >&2
    exit 1
}

pass
