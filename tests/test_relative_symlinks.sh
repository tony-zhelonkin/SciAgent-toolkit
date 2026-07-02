#!/usr/bin/env bash
# tests/test_relative_symlinks.sh
# Reproducibility fix 1: activated skill/agent/command symlinks are RELATIVE
# (and resolve to the in-repo toolkit) when the toolkit lives inside the
# project tree, and ABSOLUTE when the toolkit is external.
#
# Covers:
#   (a) in-repo toolkit (project/01_modules/SciAgent-toolkit): activate creates
#       relative links in BOTH .claude/ and .agents/ for skills/agents/commands,
#       links resolve to the in-repo toolkit, and deactivate removes them.
#   (b) external toolkit (sibling dir, override to bypass the locality guard):
#       links fall back to absolute targets.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir

# ---------------------------------------------------------------------------
# (a) in-repo toolkit → relative links
# ---------------------------------------------------------------------------
PROJ="$TMPDIR_TEST/proj"
IN_REPO="$PROJ/01_modules/SciAgent-toolkit"
mkdir -p "$IN_REPO"
build_fake_toolkit "$IN_REPO"

cd "$PROJ"
export SCIAGENT_TOOLKIT="$IN_REPO"
SCIAGENT="$IN_REPO/bin/sciagent"

"$SCIAGENT" activate base >/dev/null

# All three namespaces, both mirrors, must be relative (not starting with /).
for link in \
    .claude/skills/s_a .agents/skills/s_a \
    .claude/agents/ag_a.md .agents/agents/ag_a.md \
    .claude/commands/c_a.md .agents/commands/c_a.md; do
    assert_symlink "$link" "expected symlink at $link"
    target="$(readlink "$link")"
    case "$target" in
        /*) echo "FAIL [$_TEST_NAME] $link is ABSOLUTE ($target), expected relative" >&2; exit 1 ;;
    esac
    # The relative link must resolve to the in-repo toolkit copy.
    resolved="$(realpath "$link")"
    case "$resolved" in
        "$IN_REPO"/*) : ;;
        *) echo "FAIL [$_TEST_NAME] $link resolves to $resolved, expected under $IN_REPO" >&2; exit 1 ;;
    esac
done

# Files readable through a relative skill link.
assert_file_exists ".claude/skills/s_a/SKILL.md" "s_a not readable through relative link"

# Deactivate removes them (round-trips through manifest).
"$SCIAGENT" deactivate >/dev/null
for link in .claude/skills/s_a .agents/skills/s_a .claude/agents/ag_a.md; do
    if [[ -L "$link" || -e "$link" ]]; then
        echo "FAIL [$_TEST_NAME] $link survived deactivate" >&2
        exit 1
    fi
done

# ---------------------------------------------------------------------------
# (b) external toolkit → absolute links (guard overridden)
# ---------------------------------------------------------------------------
EXT_PROJ="$TMPDIR_TEST/ext_proj"
EXT_TK="$TMPDIR_TEST/ext-toolkit"
mkdir -p "$EXT_PROJ"
build_fake_toolkit "$EXT_TK"

cd "$EXT_PROJ"
export SCIAGENT_TOOLKIT="$EXT_TK"
EXT_SCIAGENT="$EXT_TK/bin/sciagent"

# No ./01_modules/SciAgent-toolkit here, so the guard is a no-op; links must be
# absolute because the toolkit is outside the project tree.
"$EXT_SCIAGENT" activate base >/dev/null

target="$(readlink .claude/skills/s_a)"
case "$target" in
    /*) : ;;
    *) echo "FAIL [$_TEST_NAME] external-toolkit link is relative ($target), expected absolute" >&2; exit 1 ;;
esac

pass
