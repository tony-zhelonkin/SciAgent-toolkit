#!/usr/bin/env bash
# tests/test_lint_one_door.sh — `scio lint --check harness-links` names one door
# to the catalog.
#
# A vendored project holds every skill twice: .claude/skills/ and .agents/skills/
# are the mounts the harness routes through, and 01_modules/SciAgent-toolkit/skills/
# is a real directory that ls, glob and rg reach first. The bytes are identical, so
# nothing breaks — what breaks is progressive disclosure. Routed through a mount a
# skill offers its summary and the reader pulls detail on demand; named by vendor
# path it is an ordinary file, and the whole SKILL.md lands in context with assets/
# unseen. Observed in JR-MC 2026-08-21: 6,610 bytes in two reads. The path is also
# wrong under a global install, where no vendored tree exists.
#
# Tests:
#   1. a tracked file citing the vendor skills path      -> finding, names the file
#   2. the same citation under docs/_internal/           -> silent (memory records)
#   3. the same citation under tests/                    -> silent (fixtures build it)
#   4. citing 01_modules/SciAgent-toolkit/ without skills/ -> silent (that is where
#      the toolkit lives; the defect is reaching a SKILL through it)
#   5. the capitalised 01_Modules/ spelling              -> finding
#   6. an untracked file                                 -> silent (git grep scope)
#   7. not a git repo at all                             -> silent, no error
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIO_TOOLKIT="$FAKE"
SCIO="$FAKE/bin/scio"

# _proj <dir> — a git-tracked project skeleton.
_proj() {
    mkdir -p "$1"
    git -C "$1" init -q
    git -C "$1" config user.email t@t
    git -C "$1" config user.name t
}

# _run <projdir>
_run() {
    set +e
    out=$("$SCIO" lint --check harness-links --project-dir "$1" 2>&1)
    rc=$?
    set -e
}

VENDOR='01_modules/SciAgent-toolkit/skills/louper-seurat-conversion/SKILL.md'

# --- 1. a tracked instruction reaching a skill by vendor path ----------------
P=$TMPDIR_TEST/p1; _proj "$P"
printf 'Read %s before converting.\n' "$VENDOR" > "$P/AGENTS.md"
git -C "$P" add AGENTS.md && git -C "$P" commit -qm init
_run "$P"
if ! printf '%s\n' "$out" | grep -q 'AGENTS.md reaches a skill by vendor path'; then
    echo "FAIL [$_TEST_NAME] case1: expected the vendor-path finding" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# --- 2. memory records what happened, including a wrong path -----------------
P=$TMPDIR_TEST/p2; _proj "$P"
mkdir -p "$P/docs/_internal/_project"
printf 'An agent read %s unrouted.\n' "$VENDOR" > "$P/docs/_internal/_project/session.md"
git -C "$P" add -f docs && git -C "$P" commit -qm init
_run "$P"
if printf '%s\n' "$out" | grep -q 'vendor path'; then
    echo "FAIL [$_TEST_NAME] case2: memory must be exempt" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# --- 3. a test's job is to build the state under audit ----------------------
P=$TMPDIR_TEST/p3; _proj "$P"
mkdir -p "$P/tests"
printf 'fixture="%s"\n' "$VENDOR" > "$P/tests/test_thing.sh"
git -C "$P" add tests && git -C "$P" commit -qm init
_run "$P"
if printf '%s\n' "$out" | grep -q 'vendor path'; then
    echo "FAIL [$_TEST_NAME] case3: tests must be exempt" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# --- 4. naming where the toolkit lives is not the defect --------------------
P=$TMPDIR_TEST/p4; _proj "$P"
printf 'The toolkit is vendored at 01_modules/SciAgent-toolkit/.\n' > "$P/AGENTS.md"
git -C "$P" add AGENTS.md && git -C "$P" commit -qm init
_run "$P"
if printf '%s\n' "$out" | grep -q 'vendor path'; then
    echo "FAIL [$_TEST_NAME] case4: the mount location itself is not a finding" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# --- 5. the capitalised spelling is the same defect -------------------------
P=$TMPDIR_TEST/p5; _proj "$P"
printf 'See 01_Modules/SciAgent-toolkit/skills/figure-style/SKILL.md\n' > "$P/README.md"
git -C "$P" add README.md && git -C "$P" commit -qm init
_run "$P"
if ! printf '%s\n' "$out" | grep -q 'README.md reaches a skill by vendor path'; then
    echo "FAIL [$_TEST_NAME] case5: 01_Modules/ must be caught too" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# --- 6. untracked files are outside the audit ------------------------------
P=$TMPDIR_TEST/p6; _proj "$P"
printf 'x\n' > "$P/seed.md"
git -C "$P" add seed.md && git -C "$P" commit -qm init
printf 'Read %s\n' "$VENDOR" > "$P/scratch.md"
_run "$P"
if printf '%s\n' "$out" | grep -q 'vendor path'; then
    echo "FAIL [$_TEST_NAME] case6: an untracked file must not be audited" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# --- 7. a project outside git is silence, not an error ----------------------
P=$TMPDIR_TEST/p7; mkdir -p "$P"
printf 'Read %s\n' "$VENDOR" > "$P/AGENTS.md"
_run "$P"
if [[ "$rc" -ne 0 ]] || [[ -n "$out" ]]; then
    echo "FAIL [$_TEST_NAME] case7: a non-repo project produced output" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

pass
