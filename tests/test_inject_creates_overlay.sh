#!/usr/bin/env bash
# Inject a skill on top of solo base: synthesizes _injected overlay.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

mkdir project && cd project
"$SCIAGENT" activate base >/dev/null
"$SCIAGENT" inject s_c >/dev/null

# Stack rewritten to include synthetic overlay.
assert_grep '^STACK base _injected$' .sciagent/manifest.json "synthetic overlay recorded"
assert_grep '^INJECTED _injected s_c$' .sciagent/manifest.json "INJECTED record present"

# Dual symlinks for injected skill.
assert_symlink .claude/skills/s_c
assert_symlink .agents/skills/s_c

# Managed block has the Injected subsection.
assert_grep '## Injected (overlay)' AGENTS.md "Injected subsection present"
assert_grep 's_c' AGENTS.md "injected skill listed"

# Idempotency: second inject is a no-op.
out=$("$SCIAGENT" inject s_c 2>&1)
echo "$out" | grep -q 'already injected' || {
    echo "FAIL [$_TEST_NAME] expected idempotency message, got: $out" >&2
    exit 1
}

pass
