#!/usr/bin/env bash
# tests/test_verb_help.sh — every verb answers -h/--help, and help never mutates.
#
# `activate` and `status` were the two holdouts: `-h` fell through activate's
# positional-arg loop and came back as "role not found: -h", and hit status's
# unknown-flag arm with exit 1. The dangerous half of that gap is not the
# missing text — it is a help flag that still MOUNTS. Hence the emptiness
# assertions below: a help invocation must leave the project directory exactly
# as it found it.
#
# `list` is deliberately absent from VERBS: it has never had a -h branch
# (`sciagent list -h` reports an unknown category), and giving it one is a
# separate change from this one.
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

VERBS=(activate deactivate validate lint status new craft gitignore update provision)

# ---------------------------------------------------------------------------
# 1. Every verb: -h and --help exit 0 and print a usage line naming the verb.
# ---------------------------------------------------------------------------
WORK="$TMPDIR_TEST/work"; mkdir -p "$WORK"
for v in "${VERBS[@]}"; do
    for flag in -h --help; do
        set +e
        out=$(cd "$WORK" && "$SCIAGENT" "$v" "$flag" 2>&1); rc=$?
        set -e
        [[ "$rc" -eq 0 ]] || fail "sciagent $v $flag exited $rc (expected 0)" "$out"
        [[ -n "$out" ]]   || fail "sciagent $v $flag printed nothing"
        printf '%s\n' "$out" | grep -q "sciagent $v" \
            || fail "sciagent $v $flag does not name the verb in its usage" "$out"
    done
done

# ---------------------------------------------------------------------------
# 2. Help never mutates. Run every verb's help inside a pristine project dir
#    and assert the directory is still pristine afterwards. activate is the
#    one that matters, but the invariant is worth holding for all of them.
# ---------------------------------------------------------------------------
[[ -z "$(ls -A "$WORK")" ]] \
    || fail "a --help invocation created files in the project dir: $(ls -A "$WORK" | tr '\n' ' ')"
for artifact in .claude .agents .sciagent AGENTS.md CLAUDE.md .gitignore; do
    [[ -e "$WORK/$artifact" ]] && fail "sciagent <verb> --help created $artifact"
done

# ---------------------------------------------------------------------------
# 3. activate --help specifically: documents the two positionals, omits the
#    retired output-style flag, and mounts nothing even though a role exists.
# ---------------------------------------------------------------------------
set +e
out=$(cd "$WORK" && "$SCIAGENT" activate --help 2>&1); rc=$?
set -e
[[ "$rc" -eq 0 ]] || fail "activate --help exited $rc" "$out"
if printf '%s\n' "$out" | grep -q -- '--output-style'; then
    fail "activate --help still advertises retired --output-style" "$out"
fi
printf '%s\n' "$out" | grep -qi 'overlay' \
    || fail "activate --help omits the overlay positional" "$out"
printf '%s\n' "$out" | grep -qi 'exit code' \
    || fail "activate --help omits the exit-code line the other verbs carry" "$out"
[[ -e "$WORK/.claude" || -e "$WORK/AGENTS.md" ]] \
    && fail "activate --help mounted something"

set +e
out=$(cd "$WORK" && "$SCIAGENT" activate base --output-style ghost 2>&1); rc=$?
set -e
[[ "$rc" -ne 0 ]] || fail "retired --output-style flag is still accepted" "$out"
[[ -e "$WORK/.claude" || -e "$WORK/AGENTS.md" ]] \
    && fail "rejected --output-style invocation mounted something"

# `activate -h` must not be read as a role name (the old failure mode).
set +e
out=$(cd "$WORK" && "$SCIAGENT" activate -h 2>&1); rc=$?
set -e
printf '%s\n' "$out" | grep -q 'role not found' \
    && fail "activate -h is still parsed as a role name" "$out"

# ---------------------------------------------------------------------------
# 4. status --help: documents every mode, exits 0, and does not require an
#    active stack (nor create one).
# ---------------------------------------------------------------------------
set +e
out=$(cd "$WORK" && "$SCIAGENT" status --help 2>&1); rc=$?
set -e
[[ "$rc" -eq 0 ]] || fail "status --help exited $rc" "$out"
for mode in -- --json --effective --source; do
    [[ "$mode" == "--" ]] && continue
    printf '%s\n' "$out" | grep -q -- "$mode" \
        || fail "status --help omits $mode" "$out"
done
printf '%s\n' "$out" | grep -q 'No role active' \
    && fail "status --help fell through to the actual status report" "$out"

# ---------------------------------------------------------------------------
# 5. Per-verb help composes with the dispatcher's top-level help rather than
#    shadowing it: `sciagent --help` still prints the verb table, and the two
#    outputs are different documents.
# ---------------------------------------------------------------------------
set +e
top=$("$SCIAGENT" --help 2>&1); trc=$?
set -e
[[ "$trc" -eq 0 ]] || fail "sciagent --help exited $trc" "$top"
printf '%s\n' "$top" | grep -q 'per-project AI-harness context manager' \
    || fail "top-level --help no longer prints the dispatcher usage" "$top"
vhelp=$("$SCIAGENT" status --help 2>&1)
[[ "$top" != "$vhelp" ]] || fail "status --help returned the top-level usage instead of its own"

# `sciagent help` and a bare `sciagent` still behave as before (0 / 1).
set +e
"$SCIAGENT" help >/dev/null 2>&1; hrc=$?
"$SCIAGENT"      >/dev/null 2>&1; brc=$?
set -e
[[ "$hrc" -eq 0 ]] || fail "sciagent help exited $hrc (expected 0)"
[[ "$brc" -eq 1 ]] || fail "bare sciagent exited $brc (expected 1)"

pass
