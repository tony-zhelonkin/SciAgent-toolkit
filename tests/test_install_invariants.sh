#!/usr/bin/env bash
# tests/test_install_invariants.sh
#
# ADR-D3 asks for a CHECKABLE rule, not an aspirational one: "grep install.sh
# for any harness name and the invariant either holds or it doesn't." This test
# is that grep, plus the two behavioural halves of the same rule.
#
#   1. Static — the executable body of install.sh (comments and the help text
#      stripped) contains no fetcher (curl/wget/npm/npx/git clone|fetch|pull),
#      no registry or telemetry word, and no harness name. Comments are stripped
#      because a comment cannot open a socket; the header deliberately DISCUSSES
#      curl and npm in order to say it does not use them, and a test that could
#      not tell those apart would force the honesty out of the file.
#   2. Behavioural — installing writes strictly inside --prefix, links only the
#      executable, and creates no harness or project artifact anywhere.
#   3. Behavioural — a project directory next to the install is byte-identical
#      before and after. `install.sh` puts a program on the machine; `scio`
#      decides what a project gets (00_INDEX.md §2).
set -u
. "$(dirname "$0")/_lib.sh"
. "$(dirname "$0")/_release_lib.sh"

INSTALL="$TOOLKIT_ROOT/install.sh"

setup_tmpdir
release_git_env

# --- 1. static invariants --------------------------------------------------
# Strip full-line comments and the two heredoc help texts' prose is left in —
# it is code-adjacent output, and it must not name a harness either.
code=$(grep -v '^[[:space:]]*#' "$INSTALL")

check_absent() {
    local pattern="$1" why="$2"
    if printf '%s\n' "$code" | grep -Eqi "$pattern"; then
        echo "FAIL [$_TEST_NAME] install.sh code matches /$pattern/ — $why" >&2
        printf '%s\n' "$code" | grep -Eni "$pattern" >&2
        exit 1
    fi
}

check_absent '\b(curl|wget|npm|npx|pip install|nc)\b'      "no fetcher may appear in the code path"
check_absent 'git[[:space:]]+(clone|fetch|pull|remote)'    "transport is not this script's layer (ADR-D2)"
check_absent '\b(telemetry|registry|marketplace)\b'        "no registry/telemetry (ADR-D2)"
check_absent '\b(claude|codex|opencode|gemini|cursor)\b'   "packaging must not select or configure harnesses (ADR-D3)"
check_absent 'settings\.json|CLAUDE\.md|AGENTS\.md|\.agents/|\.claude/' \
                                                           "install.sh must not know about project artifacts (ADR-D3)"

# Positive control: the grep would in fact fire if such a line existed. Without
# this, all five checks above could be passing because the pattern is broken.
if ! printf 'x=$(curl -s http://example.invalid)\n' | grep -Eqi '\b(curl|wget|npm|npx|pip install|nc)\b'; then
    echo "FAIL [$_TEST_NAME] the fetcher pattern does not match an obvious fetcher — the check is vacuous" >&2
    exit 1
fi

# --- 2/3. behavioural ------------------------------------------------------
REPO="$TMPDIR_TEST/repo"
mkdir -p "$REPO"
fixture_repo "$REPO"
release_build "$REPO" "$TMPDIR_TEST/out" HEAD >/dev/null 2>&1 \
    || { echo "FAIL [$_TEST_NAME] fixture build failed" >&2; exit 1; }
ART=$(artifact_in "$TMPDIR_TEST/out")
SHA=$(git -C "$REPO" rev-parse HEAD)

# A project that must survive untouched, and a scratch HOME that must too.
PROJ="$TMPDIR_TEST/project"
mkdir -p "$PROJ/.claude/skills" "$PROJ/.agents"
echo "# project agents" > "$PROJ/AGENTS.md"
echo '{"hooks":{}}'     > "$PROJ/.claude/settings.json"
export HOME="$TMPDIR_TEST/home"
mkdir -p "$HOME"
proj_before=$(tree_snapshot "$PROJ")
home_before=$(tree_snapshot "$HOME")

P="$TMPDIR_TEST/prefix"
( cd "$PROJ" && "$INSTALL" --archive "$ART" --checksum "$ART.sha256" --prefix "$P" >/dev/null 2>&1 ) \
    || { echo "FAIL [$_TEST_NAME] install failed" >&2; exit 1; }

assert_eq "$(tree_snapshot "$PROJ")" "$proj_before" "install.sh must never mutate a project"
assert_eq "$(tree_snapshot "$HOME")" "$home_before" "install.sh must not write into \$HOME (only --prefix)"

# Everything it wrote is under the prefix, and only three kinds of thing.
paths=$(cd "$P" && find . \( -type f -o -type l \) | LC_ALL=C sort)
while IFS= read -r p; do
    case "$p" in
        ./bin/scio) ;;
        ./share/scio/versions/"$SHA"/*) ;;
        ./share/scio/receipts/"$SHA".json) ;;
        *) echo "FAIL [$_TEST_NAME] unexpected path written by install.sh: $p" >&2; exit 1 ;;
    esac
done <<< "$paths"

# The executable is a symlink INTO the version dir — not a copy, not a wrapper
# script that could pin \$SCIO_TOOLKIT globally and subvert fleet precedence.
assert_symlink "$P/bin/scio"
assert_eq "$(readlink "$P/bin/scio")" "../share/scio/versions/$SHA/bin/scio" \
    "the linked executable must resolve into the content-addressed version dir"

pass
