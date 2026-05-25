#!/usr/bin/env bash
# tests/test_validate_skills_ref_absent.sh — when `skills-ref` is not in PATH,
# `sciagent validate` must exit 0 and produce no error output.
#
# Per kickoff.md §9 ADR-007 resolution (2026-05-24): the skills-ref check is
# silently skipped when the binary is absent. This is an explicit no-op path,
# not implicit fallthrough.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

# Ensure skills-ref is not in the PATH used for this test by prepending a
# private bin/ that does not contain it. If skills-ref is already absent from
# the real PATH this is a no-op; if it exists elsewhere the shadow hides it.
mkdir -p "$TMPDIR_TEST/empty-bin"
PATH="$TMPDIR_TEST/empty-bin:$PATH"

if command -v skills-ref >/dev/null 2>&1; then
    echo "SKIP [$_TEST_NAME] skills-ref is still reachable in sanitized PATH — cannot isolate absence" >&2
    exit 0
fi

set +e
out=$("$SCIAGENT" validate 2>&1)
rc=$?
set -e

if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] validate exited $rc (expected 0) when skills-ref is absent" >&2
    echo "--- output ---" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# No error or warning output attributable to skills-ref absence.
if printf '%s\n' "$out" | grep -qi "skills.ref"; then
    echo "FAIL [$_TEST_NAME] unexpected skills-ref mention in output when binary absent" >&2
    echo "--- output ---" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

pass
