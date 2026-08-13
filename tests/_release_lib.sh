# tests/_release_lib.sh — fixtures shared by the release/install tests.
#
# Sourced AFTER tests/_lib.sh. Not named test_*.sh, so run-all.sh never executes
# it as a test.
#
# WHY FIXTURE REPOS AND NOT THIS REPO — the recursion hazard, resolved.
# `scripts/build-release.sh` runs `tests/run-all.sh` as a release gate. A test
# that built THIS repo would therefore re-enter the suite (and, inside that
# nested suite, re-enter it again). Rather than add a `--skip-checks` escape
# hatch to the release gate — a flag that would eventually be used for a real
# release — every release test builds a throwaway git repo whose `bin/sciagent`
# and `tests/run-all.sh` are one-line stubs with a committed exit code. The gate
# path is then exercised in BOTH directions (passing and failing) in
# milliseconds, and recursion is impossible by construction.
#
# Hermetic: every fixture lives under the test's own mktemp dir, git is run with
# system/global config neutralised and with fixed identity and dates (so a
# fixture commit SHA is reproducible), and nothing reads or writes $HOME.

# shellcheck shell=bash

# release_git_env — neutralise ambient git configuration and pin identity/time.
# GIT_CONFIG_GLOBAL/SYSTEM need git >= 2.32; on older git they are ignored and
# the -c flags below still keep commits from failing on a signing default.
release_git_env() {
    export GIT_CONFIG_GLOBAL=/dev/null
    export GIT_CONFIG_SYSTEM=/dev/null
    export GIT_AUTHOR_NAME="fixture"
    export GIT_AUTHOR_EMAIL="fixture@example.invalid"
    export GIT_COMMITTER_NAME="fixture"
    export GIT_COMMITTER_EMAIL="fixture@example.invalid"
    export GIT_AUTHOR_DATE="2026-01-01T00:00:00+0000"
    export GIT_COMMITTER_DATE="2026-01-01T00:00:00+0000"
}

_fixture_git() {
    local dir="$1"; shift
    git -C "$dir" -c commit.gpgsign=false -c init.defaultBranch=main "$@"
}

# fixture_commit <repo> <message> [date]
fixture_commit() {
    local repo="$1" msg="$2" date="${3:-}"
    if [[ -n "$date" ]]; then
        GIT_AUTHOR_DATE="$date" GIT_COMMITTER_DATE="$date" \
            _fixture_git "$repo" add -A
        GIT_AUTHOR_DATE="$date" GIT_COMMITTER_DATE="$date" \
            _fixture_git "$repo" commit -q -m "$msg"
    else
        _fixture_git "$repo" add -A
        _fixture_git "$repo" commit -q -m "$msg"
    fi
}

# fixture_repo <dir> [catalog_rc] [tests_rc]
#
# Build a minimal releasable "toolkit": an executable bin/sciagent, a
# tests/run-all.sh, a copy of the real scripts/build-release.sh, a little
# tracked content, and — deliberately — a fat GITIGNORED .venv that must never
# reach an artifact (mirroring the 171 MB skill virtualenv in the real repo,
# which is the whole reason the contract mandates `git archive`).
#
# catalog_rc / tests_rc are baked into the stubs as literal exit codes, so the
# release gate's pass and fail paths are both reachable without any environment
# dependence.
fixture_repo() {
    local dir="$1" catalog_rc="${2:-0}" tests_rc="${3:-0}"

    mkdir -p "$dir"/{bin,tests,scripts,lib,docs,.venv/lib}

    cat > "$dir/bin/sciagent" <<EOF
#!/usr/bin/env bash
# fixture CLI stub — enough surface for the release gate and for install tests.
case "\${1:-}" in
    lint) echo "fixture: toolkit lint"; exit $catalog_rc ;;
    --help|help) echo "fixture sciagent"; exit 0 ;;
    *) echo "fixture sciagent: \$*"; exit 0 ;;
esac
EOF
    chmod +x "$dir/bin/sciagent"

    cat > "$dir/tests/run-all.sh" <<EOF
#!/usr/bin/env bash
echo "fixture: test suite"
exit $tests_rc
EOF
    chmod +x "$dir/tests/run-all.sh"

    cp "$TOOLKIT_ROOT/scripts/build-release.sh" "$dir/scripts/build-release.sh"
    chmod +x "$dir/scripts/build-release.sh"

    echo "fixture toolkit" > "$dir/README.md"
    echo "lib content"     > "$dir/lib/thing.sh"
    echo "doc content"     > "$dir/docs/guide.md"

    printf '%s\n' ".venv/" "dist/" > "$dir/.gitignore"

    # 2 MiB of ignored development residue. If this ever lands in an artifact,
    # the size assertions in test_build_release_contents.sh fail loudly.
    dd if=/dev/zero of="$dir/.venv/lib/big.bin" bs=1024 count=2048 status=none

    _fixture_git "$dir" init -q
    fixture_commit "$dir" "fixture release commit"
}

# release_build <repo> <outdir> [extra args...] — run the fixture's copy of
# build-release.sh. Echoes its stdout; returns its exit code.
release_build() {
    local repo="$1" out="$2"; shift 2
    ( cd "$repo" && ./scripts/build-release.sh "$@" --out "$out" --quiet )
}

# artifact_in <dir> — print the single .tar.gz in <dir>, or nothing.
artifact_in() {
    local d="$1" f
    for f in "$d"/*.tar.gz; do
        [[ -e "$f" ]] && { printf '%s\n' "$f"; return 0; }
    done
    return 1
}

sha256_hex() {
    if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | awk '{print $1}'
    else shasum -a 256 "$1" | awk '{print $1}'; fi
}

# tree_snapshot <dir> — path list + per-file content hashes, for "wrote nothing"
# assertions. Prints nothing when <dir> does not exist.
tree_snapshot() {
    local d="$1"
    [[ -d "$d" ]] || { echo "(absent)"; return 0; }
    ( cd "$d" && find . | LC_ALL=C sort && find . -type f -exec sha256sum {} \; | LC_ALL=C sort )
}
