#!/usr/bin/env bash
# tests/test_lint_internal_memory.sh — the `internal-memory` project check.
#
# Tests:
#   1. No docs/_internal/ at all → CLEAN no-op: exit 0, NO output, even under
#      --strict. Absence is the majority state across the fleet and is
#      legitimate; the check audits the SHAPE of a tree that exists.
#   2. _project/ with content plus a stage dir holding a non-empty session.md
#      → clean under --strict. A sibling stage dir holding only a flat
#      <topic>.md is equally valid. The fixture nests a git repo because that
#      is the topology ADR-D10 recommends, and test 12 covers its absence.
#   3. A stage dir holding no non-empty Markdown → the empty-scaffold finding
#      this check exists for.
#   4. A .gitkeep anywhere beneath the tree → finding.
#   5. handoffs/ as an immediate child → retired flat namespace finding.
#   6. A directory matching no stage stem → finding.
#   7. A .venv/ directory and a .parquet file → non-memory payload findings,
#      and the venv is reported once rather than per contained file.
#   8. A nested docs/_internal/.git/ → NOT a finding (ADR-D10 topology).
#   9. --strict promotes WARN to ERROR and exit 1; without it, exit 0.
#  10. The check is opt-in: `--check all --strict` never mentions it.
#  11. An unknown --check name lists internal-memory among the valid names.
#  12. A populated tree that is not its own repository → finding. The memory is
#      always its own repo: that is what makes an in-place session.md update
#      safe, and what lets its visibility be chosen apart from the code. A tree
#      of topic notes with no session.md is reported too — the requirement is
#      the topology, not the filename. Observed, not mandated: no verb creates
#      the repo.
#  13. A plan directory with no 00_INDEX.md → finding; with one → clean.
#  14. The nested repo may be a .git FILE (worktree or submodule) — not an
#      absence. A parent that tracks the tree as plain files ALONGSIDE that
#      repository is its own finding: two histories of one tree. The submodule
#      form (a gitlink) is the publication route and stays silent.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIO_TOOLKIT="$FAKE"
SCIO="$FAKE/bin/scio"

# Build a project with one stage file per requested stem.
make_proj() {
    local dir="$1"; shift
    mkdir -p "$dir/02_analysis/stages"
    local stem
    for stem in "$@"; do
        printf 'x <- 1\n' > "$dir/02_analysis/stages/$stem.R"
    done
}

# ---------------------------------------------------------------------------
# Test 1: no docs/_internal/ → silent.
# ---------------------------------------------------------------------------
P1="$TMPDIR_TEST/p1"
make_proj "$P1" 30_grn
mkdir -p "$P1/docs"

set +e
out1=$("$SCIO" lint --check internal-memory --strict --project-dir "$P1" 2>&1)
rc1=$?
set -e
if [[ "$rc1" -ne 0 || -n "$out1" ]]; then
    echo "FAIL [$_TEST_NAME] test1: absent docs/_internal/ must be silent (rc=$rc1)" >&2
    printf '%s\n' "$out1" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 2: a well-formed tree passes under --strict.
# ---------------------------------------------------------------------------
P2="$TMPDIR_TEST/p2"
make_proj "$P2" 30_grn 40_peaks
mkdir -p "$P2/docs/_internal/_project" \
         "$P2/docs/_internal/30_grn" \
         "$P2/docs/_internal/40_peaks"
printf 'Cohort scope decided 2026-08-01.\n' > "$P2/docs/_internal/_project/session.md"
printf 'Stopped after the GRN pass.\n' > "$P2/docs/_internal/30_grn/session.md"
printf 'Why the peak floor is 0.05.\n' > "$P2/docs/_internal/40_peaks/floor.md"
git -C "$P2/docs/_internal" init -q

set +e
out2=$("$SCIO" lint --check internal-memory --strict --project-dir "$P2" 2>&1)
rc2=$?
set -e
if [[ "$rc2" -ne 0 || -n "$out2" ]]; then
    echo "FAIL [$_TEST_NAME] test2: a well-formed memory tree must be clean (rc=$rc2)" >&2
    printf '%s\n' "$out2" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 3: a stage dir with no content, and test 9's severity pair.
# ---------------------------------------------------------------------------
P3="$TMPDIR_TEST/p3"
make_proj "$P3" 30_grn
mkdir -p "$P3/docs/_internal/30_grn"

set +e
out3=$("$SCIO" lint --check internal-memory --project-dir "$P3" 2>&1)
rc3=$?
set -e
if [[ "$rc3" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] test3: soft warn must exit 0, got $rc3" >&2
    printf '%s\n' "$out3" >&2
    exit 1
fi
if ! printf '%s\n' "$out3" | grep -q 'WARN internal-memory: docs/_internal/30_grn/ holds no non-empty Markdown'; then
    echo "FAIL [$_TEST_NAME] test3: expected the empty-scaffold WARN" >&2
    printf '%s\n' "$out3" >&2
    exit 1
fi

set +e
out3s=$("$SCIO" lint --check internal-memory --strict --project-dir "$P3" 2>&1)
rc3s=$?
set -e
if [[ "$rc3s" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] test9: --strict must promote to exit 1" >&2
    printf '%s\n' "$out3s" >&2
    exit 1
fi
if ! printf '%s\n' "$out3s" | grep -q 'ERROR internal-memory: docs/_internal/30_grn/ holds no non-empty Markdown'; then
    echo "FAIL [$_TEST_NAME] test9: --strict must emit ERROR" >&2
    printf '%s\n' "$out3s" >&2
    exit 1
fi

# An empty file does not count as content.
printf '' > "$P3/docs/_internal/30_grn/session.md"
set +e
out3e=$("$SCIO" lint --check internal-memory --project-dir "$P3" 2>&1)
set -e
if ! printf '%s\n' "$out3e" | grep -q 'WARN internal-memory: docs/_internal/30_grn/ holds no non-empty Markdown'; then
    echo "FAIL [$_TEST_NAME] test3: an empty session.md must not satisfy the content rule" >&2
    printf '%s\n' "$out3e" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 4: a .gitkeep anywhere beneath the tree.
# ---------------------------------------------------------------------------
P4="$TMPDIR_TEST/p4"
make_proj "$P4" 30_grn
mkdir -p "$P4/docs/_internal/30_grn"
printf 'Stopped mid-pass.\n' > "$P4/docs/_internal/30_grn/session.md"
touch "$P4/docs/_internal/30_grn/.gitkeep"

set +e
out4=$("$SCIO" lint --check internal-memory --project-dir "$P4" 2>&1)
set -e
if ! printf '%s\n' "$out4" | grep -q 'WARN internal-memory: docs/_internal/30_grn/.gitkeep claims a directory'; then
    echo "FAIL [$_TEST_NAME] test4: expected a .gitkeep finding" >&2
    printf '%s\n' "$out4" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 5: handoffs/ as an immediate child.
# ---------------------------------------------------------------------------
P5="$TMPDIR_TEST/p5"
make_proj "$P5" 30_grn
mkdir -p "$P5/docs/_internal/handoffs"
printf 'Where the work stopped.\n' > "$P5/docs/_internal/handoffs/2026-08-01.md"

set +e
out5=$("$SCIO" lint --check internal-memory --project-dir "$P5" 2>&1)
set -e
if ! printf '%s\n' "$out5" | grep -q 'WARN internal-memory: docs/_internal/handoffs/ is a retired flat namespace'; then
    echo "FAIL [$_TEST_NAME] test5: expected the retired-namespace finding" >&2
    printf '%s\n' "$out5" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 6: a directory matching no stage stem.
# ---------------------------------------------------------------------------
P6="$TMPDIR_TEST/p6"
make_proj "$P6" 30_grn
mkdir -p "$P6/docs/_internal/99_nowhere"
printf 'Notes with no stage.\n' > "$P6/docs/_internal/99_nowhere/session.md"

set +e
out6=$("$SCIO" lint --check internal-memory --project-dir "$P6" 2>&1)
set -e
if ! printf '%s\n' "$out6" | grep -q 'WARN internal-memory: docs/_internal/99_nowhere/ matches no stage stem'; then
    echo "FAIL [$_TEST_NAME] test6: expected the unmatched-stem finding" >&2
    printf '%s\n' "$out6" >&2
    exit 1
fi

# The 02_analysis/scripts/ spelling is accepted for the pre-rename window.
P6B="$TMPDIR_TEST/p6b"
mkdir -p "$P6B/02_analysis/scripts" "$P6B/docs/_internal/30_grn"
printf 'x <- 1\n' > "$P6B/02_analysis/scripts/30_grn.R"
printf 'Stopped after the GRN pass.\n' > "$P6B/docs/_internal/30_grn/session.md"
git -C "$P6B/docs/_internal" init -q

set +e
out6b=$("$SCIO" lint --check internal-memory --strict --project-dir "$P6B" 2>&1)
rc6b=$?
set -e
if [[ "$rc6b" -ne 0 || -n "$out6b" ]]; then
    echo "FAIL [$_TEST_NAME] test6b: the 02_analysis/scripts/ spelling must be accepted (rc=$rc6b)" >&2
    printf '%s\n' "$out6b" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 7: non-memory payloads, and the venv reported once.
# ---------------------------------------------------------------------------
P7="$TMPDIR_TEST/p7"
make_proj "$P7" 30_grn
mkdir -p "$P7/docs/_internal/30_grn" "$P7/docs/_internal/_project/.venv/lib"
printf 'Stopped mid-pass.\n' > "$P7/docs/_internal/30_grn/session.md"
printf 'x\n' > "$P7/docs/_internal/_project/.venv/lib/site.pyc"
printf 'x\n' > "$P7/docs/_internal/_project/.venv/pyvenv.cfg"
printf 'x\n' > "$P7/docs/_internal/30_grn/cells.parquet"

set +e
out7=$("$SCIO" lint --check internal-memory --project-dir "$P7" 2>&1)
set -e
if ! printf '%s\n' "$out7" | grep -q 'docs/_internal/_project/.venv is not memory'; then
    echo "FAIL [$_TEST_NAME] test7: expected the virtualenv finding" >&2
    printf '%s\n' "$out7" >&2
    exit 1
fi
if ! printf '%s\n' "$out7" | grep -q 'docs/_internal/30_grn/cells.parquet is not memory'; then
    echo "FAIL [$_TEST_NAME] test7: expected the parquet finding" >&2
    printf '%s\n' "$out7" >&2
    exit 1
fi
if printf '%s\n' "$out7" | grep -q 'site.pyc'; then
    echo "FAIL [$_TEST_NAME] test7: a matched directory must be reported once, not descended into" >&2
    printf '%s\n' "$out7" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 8: a nested docs/_internal/.git/ is the recommended topology.
# ---------------------------------------------------------------------------
P8="$TMPDIR_TEST/p8"
make_proj "$P8" 30_grn
mkdir -p "$P8/docs/_internal/30_grn"
printf 'Stopped mid-pass.\n' > "$P8/docs/_internal/30_grn/session.md"
git -C "$P8/docs/_internal" init -q

set +e
out8=$("$SCIO" lint --check internal-memory --strict --project-dir "$P8" 2>&1)
rc8=$?
set -e
if [[ "$rc8" -ne 0 || -n "$out8" ]]; then
    echo "FAIL [$_TEST_NAME] test8: a nested memory repo must not be reported (rc=$rc8)" >&2
    printf '%s\n' "$out8" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 8b: a continuity filename outside the live-plus-archive shape.
# ---------------------------------------------------------------------------
P8B="$TMPDIR_TEST/p8b"
make_proj "$P8B" 30_grn
mkdir -p "$P8B/docs/_internal/30_grn"
printf 'Current.\n' > "$P8B/docs/_internal/30_grn/session.md"
printf 'Older.\n' > "$P8B/docs/_internal/30_grn/session_20260801.md"

set +e
out8b=$("$SCIO" lint --check internal-memory --project-dir "$P8B" 2>&1)
set -e
if ! printf '%s\n' "$out8b" | grep -q 'session_20260801.md: retired continuity filename'; then
    echo "FAIL [$_TEST_NAME] test8b: expected the retired-filename finding" >&2
    printf '%s\n' "$out8b" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 8c: session.md beside its dated archives is the CORRECT state. A handoff
# rewrites the live file and keeps the copy it supersedes, stamped with the last
# day that copy was actual, so several files here must stay silent. The -2 form
# covers a second handoff on the same day.
# ---------------------------------------------------------------------------
P8C="$TMPDIR_TEST/p8c"
make_proj "$P8C" 30_grn
mkdir -p "$P8C/docs/_internal/30_grn"
git -C "$P8C/docs/_internal" init -q 2>/dev/null
printf 'Current.\n'          > "$P8C/docs/_internal/30_grn/session.md"
printf 'Actual to Aug 1.\n'  > "$P8C/docs/_internal/30_grn/session-2026-08-01.md"
printf 'Later that day.\n'   > "$P8C/docs/_internal/30_grn/session-2026-08-01-2.md"
printf 'Actual to Aug 9.\n'  > "$P8C/docs/_internal/30_grn/session-2026-08-09.md"

set +e
out8c=$("$SCIO" lint --check internal-memory --project-dir "$P8C" 2>&1)
rc8c=$?
set -e
if [[ "$rc8c" -ne 0 ]] || [[ -n "$out8c" ]]; then
    echo "FAIL [$_TEST_NAME] test8c: archives beside a live session.md must be silent" >&2
    printf '%s\n' "$out8c" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 8d: archives with no live session.md — the rewrite lost its subject, and
# a reader opening the scope finds only records that stopped being true.
# ---------------------------------------------------------------------------
P8D="$TMPDIR_TEST/p8d"
make_proj "$P8D" 30_grn
mkdir -p "$P8D/docs/_internal/30_grn"
git -C "$P8D/docs/_internal" init -q 2>/dev/null
printf 'Actual to Aug 1.\n' > "$P8D/docs/_internal/30_grn/session-2026-08-01.md"

set +e
out8d=$("$SCIO" lint --check internal-memory --project-dir "$P8D" 2>&1)
set -e
if ! printf '%s\n' "$out8d" | grep -q 'archives with no live session.md'; then
    echo "FAIL [$_TEST_NAME] test8d: expected the orphaned-archive finding" >&2
    printf '%s\n' "$out8d" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 10: opt-in — `--check all --strict` never runs it.
# ---------------------------------------------------------------------------
set +e
out10=$("$SCIO" lint --check all --strict --project-dir "$P3" 2>&1)
set -e
if printf '%s\n' "$out10" | grep -q 'internal-memory'; then
    echo "FAIL [$_TEST_NAME] test10: internal-memory must not be a member of all" >&2
    printf '%s\n' "$out10" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 11: the name is selectable and appears in the valid-name list.
# ---------------------------------------------------------------------------
set +e
out11=$("$SCIO" lint --check bogus --project-dir "$P1" 2>&1)
rc11=$?
set -e
if [[ "$rc11" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] test11: an unknown --check name must exit 1" >&2
    exit 1
fi
if ! printf '%s\n' "$out11" | grep -q 'internal-memory'; then
    echo "FAIL [$_TEST_NAME] test11: the valid-name list must include internal-memory" >&2
    printf '%s\n' "$out11" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 12: a populated tree that is not its own repository.
# ---------------------------------------------------------------------------
P12="$TMPDIR_TEST/p12"
make_proj "$P12" 30_grn
mkdir -p "$P12/docs/_internal/30_grn"
printf 'Stopped mid-pass.\n' > "$P12/docs/_internal/30_grn/session.md"

set +e
out12=$("$SCIO" lint --check internal-memory --project-dir "$P12" 2>&1)
set -e
if ! printf '%s\n' "$out12" | grep -q 'WARN internal-memory: docs/_internal/ is not its own repository'; then
    echo "FAIL [$_TEST_NAME] test12: expected the not-a-repository finding" >&2
    printf '%s\n' "$out12" >&2
    exit 1
fi

# Topic notes with no session.md are reported too: the requirement is that the
# memory is its own repository, which does not depend on which file is present.
P12B="$TMPDIR_TEST/p12b"
make_proj "$P12B" 30_grn
mkdir -p "$P12B/docs/_internal/30_grn"
printf 'Why the floor is 0.05.\n' > "$P12B/docs/_internal/30_grn/floor.md"

set +e
out12b=$("$SCIO" lint --check internal-memory --project-dir "$P12B" 2>&1)
set -e
if ! printf '%s\n' "$out12b" | grep -q 'is not its own repository'; then
    echo "FAIL [$_TEST_NAME] test12b: a topic-note tree with no repo must be reported" >&2
    printf '%s\n' "$out12b" >&2
    exit 1
fi

# Giving it a repository clears the finding.
git -C "$P12B/docs/_internal" init -q
set +e
out12c=$("$SCIO" lint --check internal-memory --strict --project-dir "$P12B" 2>&1)
rc12c=$?
set -e
if [[ "$rc12c" -ne 0 || -n "$out12c" ]]; then
    echo "FAIL [$_TEST_NAME] test12c: a nested repo must clear the finding (rc=$rc12c)" >&2
    printf '%s\n' "$out12c" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 13: a plan directory earns itself with a phase map.
# ---------------------------------------------------------------------------
P13="$TMPDIR_TEST/p13"
make_proj "$P13" 30_grn
mkdir -p "$P13/docs/_internal/_project/plans/2026-08-20-slug"
printf 'Phase one.\n' > "$P13/docs/_internal/_project/plans/2026-08-20-slug/01_first.md"
git -C "$P13/docs/_internal" init -q

set +e
out13=$("$SCIO" lint --check internal-memory --project-dir "$P13" 2>&1)
set -e
if ! printf '%s\n' "$out13" | grep -q 'docs/_internal/_project/plans/2026-08-20-slug needs a non-empty 00_INDEX.md'; then
    echo "FAIL [$_TEST_NAME] test13: expected the missing-index finding" >&2
    printf '%s\n' "$out13" >&2
    exit 1
fi

printf 'The phase map.\n' > "$P13/docs/_internal/_project/plans/2026-08-20-slug/00_INDEX.md"
set +e
out13b=$("$SCIO" lint --check internal-memory --strict --project-dir "$P13" 2>&1)
rc13b=$?
set -e
if [[ "$rc13b" -ne 0 || -n "$out13b" ]]; then
    echo "FAIL [$_TEST_NAME] test13b: a plan with an index must be clean (rc=$rc13b)" >&2
    printf '%s\n' "$out13b" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 14: the two other ways history can already exist.
# ---------------------------------------------------------------------------
# (a) .git as a FILE, which is what a worktree or a submodule leaves behind.
P14="$TMPDIR_TEST/p14"
make_proj "$P14" 30_grn
mkdir -p "$P14/docs/_internal/30_grn"
printf 'Stopped mid-pass.\n' > "$P14/docs/_internal/30_grn/session.md"
printf 'gitdir: /elsewhere/memory.git\n' > "$P14/docs/_internal/.git"

set +e
out14=$("$SCIO" lint --check internal-memory --project-dir "$P14" 2>&1)
set -e
if printf '%s\n' "$out14" | grep -q 'not its own repository'; then
    echo "FAIL [$_TEST_NAME] test14a: a .git file is a nested repo, not an absence" >&2
    printf '%s\n' "$out14" >&2
    exit 1
fi

# (b) the parent force-adds the tree as plain files while the nested repository
# also versions it: two histories of one tree, which is its own finding.
P14B="$TMPDIR_TEST/p14b"
make_proj "$P14B" 30_grn
mkdir -p "$P14B/docs/_internal/30_grn"
printf 'Stopped mid-pass.\n' > "$P14B/docs/_internal/30_grn/session.md"
# Order matters: git refuses to add plain paths inside an existing embedded
# repository, so the parent tracks them first and the nested repo arrives after.
git -C "$P14B" init -q
git -C "$P14B" config user.email "test@example.com"
git -C "$P14B" config user.name "Test"
printf 'docs/_internal/\n' > "$P14B/.gitignore"
git -C "$P14B" add -f docs/_internal/30_grn/session.md
git -C "$P14B" commit -q -m "memory" >/dev/null
git -C "$P14B/docs/_internal" init -q

set +e
out14b=$("$SCIO" lint --check internal-memory --project-dir "$P14B" 2>&1)
set -e
if ! printf '%s\n' "$out14b" | grep -q 'two histories of one tree'; then
    echo "FAIL [$_TEST_NAME] test14b: plain-tracked files beside the nested repo must be reported" >&2
    printf '%s\n' "$out14b" >&2
    exit 1
fi

pass
