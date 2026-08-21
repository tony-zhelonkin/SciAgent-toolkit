#!/usr/bin/env bash
# tests/test_lint_harness_links.sh — the `harness-links` project check, and the
# relative targets `link` now writes (#37).
#
# Tests:
#   1. A project with no mounts at all → CLEAN no-op, even under --strict.
#   2. `scio link` writes RELATIVE targets for all six category mounts when the
#      toolkit is vendored inside the project — the fix for #37. Verified by
#      readlink, not by resolution: an absolute link resolves correctly on the
#      side that wrote it, which is exactly why the bug was invisible.
#   3. A freshly linked project passes `--check harness-links --strict`.
#   4. An absolute in-project target → finding, and `--strict` promotes it.
#   5. Re-running `link` REPAIRS that absolute target back to relative. Without
#      this, the fix would never reach the 25 already-bound copies.
#   6. A dangling mount → finding naming the unresolved target.
#   7. An absolute target OUTSIDE the project → silent. That is the
#      global-install channel, where relative would break on a move.
#   8. `link` is idempotent: a second run reports no replacement.
#   9. harness-links is a member of `all`.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir

# A project with the toolkit vendored inside it, as the canonical channel does.
# bin/scio resolves its own symlinks to find the toolkit root, so the fake copy
# is selected through SCIO_TOOLKIT rather than by invoking its shimmed bin/.
PROJ="$TMPDIR_TEST/proj"
mkdir -p "$PROJ/01_modules"
FAKE="$PROJ/01_modules/SciAgent-toolkit"
build_fake_toolkit "$FAKE"
ln -sfn "$TOOLKIT_ROOT/templates" "$FAKE/templates"
SCIO="$TOOLKIT_ROOT/bin/scio"
export SCIO_TOOLKIT="$FAKE"

# ---------------------------------------------------------------------------
# Test 1: nothing mounted → silent.
# ---------------------------------------------------------------------------
BARE="$TMPDIR_TEST/bare"
mkdir -p "$BARE"
set +e
out1=$("$SCIO" lint --check harness-links --strict --project-dir "$BARE" 2>&1)
rc1=$?
set -e
if [[ "$rc1" -ne 0 || -n "$out1" ]]; then
    echo "FAIL [$_TEST_NAME] test1: an unmounted project must be silent (rc=$rc1)" >&2
    printf '%s\n' "$out1" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 2: link writes relative targets.
# ---------------------------------------------------------------------------
set +e
outlink=$(cd "$PROJ" && "$SCIO" link 2>&1)
rclink=$?
set -e
if [[ "$rclink" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] test2: scio link failed (rc=$rclink)" >&2
    printf '%s\n' "$outlink" >&2
    exit 1
fi

for m in .claude/skills .claude/agents .claude/commands \
         .agents/skills .agents/agents .agents/commands; do
    if [[ ! -L "$PROJ/$m" ]]; then
        echo "FAIL [$_TEST_NAME] test2: $m is not a symlink" >&2
        exit 1
    fi
    t=$(readlink "$PROJ/$m")
    case "$t" in
        /*)
            echo "FAIL [$_TEST_NAME] test2: $m target is absolute: $t" >&2
            exit 1
            ;;
    esac
    if [[ ! -e "$PROJ/$m" ]]; then
        echo "FAIL [$_TEST_NAME] test2: $m is relative but does not resolve: $t" >&2
        exit 1
    fi
done

# ---------------------------------------------------------------------------
# Test 3: a freshly linked project is clean under --strict.
# ---------------------------------------------------------------------------
set +e
out3=$("$SCIO" lint --check harness-links --strict --project-dir "$PROJ" 2>&1)
rc3=$?
set -e
if [[ "$rc3" -ne 0 || -n "$out3" ]]; then
    echo "FAIL [$_TEST_NAME] test3: a freshly linked project must be clean (rc=$rc3)" >&2
    printf '%s\n' "$out3" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 4: an absolute in-project target is a finding; --strict promotes.
# ---------------------------------------------------------------------------
ln -sfn "$FAKE/skills" "$PROJ/.claude/skills"

set +e
out4=$("$SCIO" lint --check harness-links --project-dir "$PROJ" 2>&1)
rc4=$?
set -e
if [[ "$rc4" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] test4: a soft warn must exit 0, got $rc4" >&2
    printf '%s\n' "$out4" >&2
    exit 1
fi
if ! printf '%s\n' "$out4" | grep -q 'WARN harness-links: .claude/skills names an absolute path inside this project'; then
    echo "FAIL [$_TEST_NAME] test4: expected the absolute-target WARN" >&2
    printf '%s\n' "$out4" >&2
    exit 1
fi

set +e
out4s=$("$SCIO" lint --check harness-links --strict --project-dir "$PROJ" 2>&1)
rc4s=$?
set -e
if [[ "$rc4s" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] test4: --strict must exit 1" >&2
    exit 1
fi
if ! printf '%s\n' "$out4s" | grep -q 'ERROR harness-links: .claude/skills names an absolute path'; then
    echo "FAIL [$_TEST_NAME] test4: --strict must emit ERROR" >&2
    printf '%s\n' "$out4s" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 5: re-running link repairs the absolute target.
# ---------------------------------------------------------------------------
set +e
out5=$(cd "$PROJ" && "$SCIO" link 2>&1)
rc5=$?
set -e
if [[ "$rc5" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] test5: repair link run failed (rc=$rc5)" >&2
    printf '%s\n' "$out5" >&2
    exit 1
fi
t5=$(readlink "$PROJ/.claude/skills")
case "$t5" in
    /*)
        echo "FAIL [$_TEST_NAME] test5: link did not repair the absolute target: $t5" >&2
        exit 1
        ;;
esac
if ! printf '%s\n' "$out5" | grep -q 'replaced link: .claude/skills'; then
    echo "FAIL [$_TEST_NAME] test5: the repair must be reported, not silent" >&2
    printf '%s\n' "$out5" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 6: a dangling mount.
# ---------------------------------------------------------------------------
DANGLE="$TMPDIR_TEST/dangle"
mkdir -p "$DANGLE/.claude"
ln -s ../nowhere/skills "$DANGLE/.claude/skills"

set +e
out6=$("$SCIO" lint --check harness-links --project-dir "$DANGLE" 2>&1)
set -e
if ! printf '%s\n' "$out6" | grep -q 'WARN harness-links: .claude/skills does not resolve'; then
    echo "FAIL [$_TEST_NAME] test6: expected the dangling-mount WARN" >&2
    printf '%s\n' "$out6" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 7: an absolute target outside the project is the global-install channel.
# ---------------------------------------------------------------------------
GLOBAL="$TMPDIR_TEST/global"
mkdir -p "$GLOBAL/.claude"
ln -s "$FAKE/skills" "$GLOBAL/.claude/skills"

set +e
out7=$("$SCIO" lint --check harness-links --strict --project-dir "$GLOBAL" 2>&1)
rc7=$?
set -e
if [[ "$rc7" -ne 0 || -n "$out7" ]]; then
    echo "FAIL [$_TEST_NAME] test7: an out-of-project absolute target must be silent (rc=$rc7)" >&2
    printf '%s\n' "$out7" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 8: link stays idempotent.
# ---------------------------------------------------------------------------
set +e
out8=$(cd "$PROJ" && "$SCIO" link 2>&1)
set -e
if printf '%s\n' "$out8" | grep -qE 'replaced link|^linked:'; then
    echo "FAIL [$_TEST_NAME] test8: a second link run must not rewrite the mounts" >&2
    printf '%s\n' "$out8" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 9: harness-links runs under `all`.
# ---------------------------------------------------------------------------
ln -sfn "$FAKE/skills" "$PROJ/.claude/skills"
set +e
out9=$("$SCIO" lint --check all --project-dir "$PROJ" 2>&1)
set -e
if ! printf '%s\n' "$out9" | grep -q 'harness-links'; then
    echo "FAIL [$_TEST_NAME] test9: harness-links must be a member of all" >&2
    printf '%s\n' "$out9" >&2
    exit 1
fi

pass
