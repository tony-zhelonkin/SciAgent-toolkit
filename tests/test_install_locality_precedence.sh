#!/usr/bin/env bash
# tests/test_install_locality_precedence.sh
#
# "Fleet precedence is unchanged and must stay that way: a project-local
# 01_modules/SciAgent-toolkit overrides any global installation"
# (10_packaging_contracts.md §3). Exact-commit reproducibility for the 22 fleet
# projects rests on that, so an installed layout that could subvert it would be
# a silent correctness bug, not a packaging nit.
#
# This is the end-to-end check, not a code reading: a REAL `bin/scio` plus a
# real `lib/scio/` are packaged into a release, installed into a scratch
# prefix, and then invoked from a project.
#
#   - In a project that ships its own toolkit, the globally installed
#     `scio link` REFUSES, because the
#     symlinked executable resolves $SCIO_TOOLKIT to the version directory,
#     which is not the project's in-repo toolkit.
#   - In a project that ships none, the same global binary works. Without this
#     control the refusal above could be caused by a broken install rather than
#     by precedence.
set -u
. "$(dirname "$0")/_lib.sh"
. "$(dirname "$0")/_release_lib.sh"

INSTALL="$TOOLKIT_ROOT/install.sh"

setup_tmpdir
release_git_env

# --- a fixture release carrying the REAL dispatcher and libs ---------------
SRC="$TMPDIR_TEST/src"
build_fake_toolkit "$SRC"                       # fixture skills/agents/commands
rm -f "$SRC/lib/scio" "$SRC/bin/scio"   # build_fake_toolkit symlinks these
cp -R "$TOOLKIT_ROOT/lib/scio" "$SRC/lib/scio"
cp    "$TOOLKIT_ROOT/bin/scio" "$SRC/bin/scio"
chmod +x "$SRC/bin/scio"

mkdir -p "$SRC/tests" "$SRC/scripts"
printf '#!/usr/bin/env bash\necho "fixture: test suite"\nexit 0\n' > "$SRC/tests/run-all.sh"
chmod +x "$SRC/tests/run-all.sh"
cp "$TOOLKIT_ROOT/scripts/build-release.sh" "$SRC/scripts/build-release.sh"
chmod +x "$SRC/scripts/build-release.sh"
printf '%s\n' "dist/" > "$SRC/.gitignore"

git -C "$SRC" -c init.defaultBranch=main init -q
git -C "$SRC" -c commit.gpgsign=false add -A
git -C "$SRC" -c commit.gpgsign=false commit -q -m "packaged toolkit fixture"
SHA=$(git -C "$SRC" rev-parse HEAD)

release_build "$SRC" "$TMPDIR_TEST/out" HEAD >/dev/null 2>&1 \
    || { echo "FAIL [$_TEST_NAME] packaging the real dispatcher failed" >&2; exit 1; }
ART=$(artifact_in "$TMPDIR_TEST/out")

P="$TMPDIR_TEST/prefix"
"$INSTALL" --archive "$ART" --checksum "$ART.sha256" --prefix "$P" >/dev/null 2>&1 \
    || { echo "FAIL [$_TEST_NAME] install failed" >&2; exit 1; }
GLOBAL="$P/bin/scio"
assert_symlink "$GLOBAL"

# --- control: no in-repo toolkit → the global install is usable ------------
PLAIN="$TMPDIR_TEST/plain-project"
mkdir -p "$PLAIN"
out=$( cd "$PLAIN" && env -u SCIO_TOOLKIT "$GLOBAL" link 2>&1 )
rc=$?
if [[ $rc -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] the globally installed scio could not link a plain project" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi
assert_symlink "$PLAIN/.claude/skills" "control link produced no skill binding"

# --- the real claim: a project shipping its own toolkit wins ---------------
FLEET="$TMPDIR_TEST/fleet-project"
mkdir -p "$FLEET/01_modules"
cp -R "$SRC" "$FLEET/01_modules/SciAgent-toolkit"
rm -rf "$FLEET/01_modules/SciAgent-toolkit/.git"

out=$( cd "$FLEET" && env -u SCIO_TOOLKIT "$GLOBAL" link 2>&1 )
rc=$?
if [[ $rc -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] the global install mutated a project that ships its own toolkit" >&2
    echo "  Fleet precedence is what keeps a published analysis pinned to its commit." >&2
    printf '%s\n' "$out" >&2
    exit 1
fi
printf '%s\n' "$out" | grep -q "refusing to link against an external toolkit" \
    || { echo "FAIL [$_TEST_NAME] refusal did not come from the locality guard" >&2; printf '%s\n' "$out" >&2; exit 1; }
[[ -e "$FLEET/.claude/skills" ]] \
    && { echo "FAIL [$_TEST_NAME] the refused run still created bindings" >&2; exit 1; }

# The project's OWN toolkit can of course do it.
out=$( cd "$FLEET" && env -u SCIO_TOOLKIT ./01_modules/SciAgent-toolkit/bin/scio link 2>&1 )
rc=$?
if [[ $rc -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] the project's in-repo toolkit could not link it" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi
assert_symlink "$FLEET/.claude/skills" "in-repo link produced no skill binding"

pass
