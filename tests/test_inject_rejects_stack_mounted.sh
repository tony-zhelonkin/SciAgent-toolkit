#!/usr/bin/env bash
# Regression: `sciagent inject <skill>` against a skill already mounted via
# the active role stack must not create a duplicate injected manifest entry.
# Prior bug: the idempotency check keyed on (target_overlay, skill); since the
# synthetic `_injected` overlay differs from `base`, a duplicate slipped
# through and a subsequent `eject` would tear down the base symlink.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

mkdir project && cd project
"$SCIAGENT" activate base >/dev/null

# s_a is in roles/base.yaml. Inject must refuse-but-exit-0 (idempotent
# end-state) and must not extend the injected[] array.
out=$("$SCIAGENT" inject s_a 2>&1)
rc=$?

if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] expected exit 0 (idempotent) injecting stack-mounted s_a, got $rc" >&2
    echo "output: $out" >&2
    exit 1
fi

printf '%s\n' "$out" | grep -qi 'already mounted' || {
    echo "FAIL [$_TEST_NAME] expected message naming stack-mount; got: $out" >&2
    exit 1
}

# Manifest must have no injected entries for s_a.
inj_count=$(grep -c '"skill": *"s_a"' .sciagent/manifest.json 2>/dev/null || true)
if [[ "${inj_count:-0}" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] injected[] gained an s_a entry (count=$inj_count)" >&2
    cat .sciagent/manifest.json >&2
    exit 1
fi

# Symlinks must not be duplicated either: exactly the two from base activation.
sym_count=$(grep -Eoc '"\.(claude|agents)/skills/s_a"' .sciagent/manifest.json || true)
if [[ "$sym_count" -ne 2 ]]; then
    echo "FAIL [$_TEST_NAME] expected 2 s_a symlinks in manifest, got $sym_count" >&2
    cat .sciagent/manifest.json >&2
    exit 1
fi

# A follow-up eject must still refuse (because s_a is stack-mounted), and the
# base symlinks must remain intact afterwards.
eject_out=$("$SCIAGENT" eject s_a 2>&1)
eject_rc=$?
if [[ "$eject_rc" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] eject of stack-mounted s_a should fail; got 0" >&2
    echo "output: $eject_out" >&2
    exit 1
fi

assert_symlink .claude/skills/s_a "base symlink damaged after refused eject"
assert_symlink .agents/skills/s_a "base agents symlink damaged after refused eject"

pass
