#!/usr/bin/env bash
# tests/test_teardown_target_ownership.sh
# Phase 5c follow-up ("Option D"): symlink_teardown_all derives ownership from
# each mount's resolved LINK TARGET, not from the manifest's bookkeeping.
# Covers the four scenarios the manifest-driven teardown could get wrong:
#   (a) a mount the manifest has lost track of is still cleaned up (recovery);
#   (b) a BROKEN toolkit-owned link (canonical file deleted out from under a
#       live mount) is still recognised and removed, not stranded;
#   (c) a user symlink pointing OUTSIDE the toolkit survives, even though it
#       sits right next to toolkit mounts in the same directory;
#   (d) a real (non-symlink) user file in each of the four `.claude/` trees
#       survives untouched.
# Each assertion below is paired with a fixture actually capable of failing
# it — e.g. (a) is proven by hand-truncating the manifest's "symlinks" array
# before deactivating, so a still-manifest-driven teardown would leave that
# entry behind.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

mkdir project && cd project

"$SCIAGENT" activate base >/dev/null
assert_symlink .claude/skills/s_a "sanity: s_a mounted"
assert_symlink .claude/skills/s_b "sanity: s_b mounted"
assert_symlink .claude/agents/ag_a.md "sanity: ag_a mounted"

# --- (a) manifest loses track of a mount; teardown must recover it anyway ---
# Drop .claude/skills/s_b's entry from the manifest's "symlinks" array (both
# the string and its dual .agents/ mirror) — simulating a crashed activate,
# hand-edited manifest, or manual repair that never got recorded.
python3 - <<'PYEOF'
import json
m = json.load(open('.sciagent/manifest.json'))
m['symlinks'] = [s for s in m['symlinks'] if 's_b' not in s]
json.dump(m, open('.sciagent/manifest.json', 'w'), indent=2)
PYEOF
if grep -q 's_b' .sciagent/manifest.json; then
    echo "FAIL [$_TEST_NAME] fixture setup: s_b still referenced in manifest" >&2
    exit 1
fi
assert_symlink .claude/skills/s_b "fixture: s_b still mounted, but now absent from manifest"

# --- (b) a BROKEN toolkit-owned link must still be recognised and removed ---
# Delete the canonical skill directory s_c points at, so .claude/skills/s_c
# becomes dangling while remaining a toolkit-owned link (target lexically
# under $SCIAGENT_TOOLKIT even though the final component no longer exists).
rm -rf "$FAKE/skills/s_c"
if [[ -e .claude/skills/s_c ]]; then
    echo "FAIL [$_TEST_NAME] fixture setup: s_c should be dangling now" >&2
    exit 1
fi
assert_symlink .claude/skills/s_c "fixture: s_c is a dangling-but-present symlink"

# --- (c) a user symlink pointing OUTSIDE the toolkit must survive ---
mkdir -p ../outside
echo "user content" > ../outside/note.md
ln -s "$PWD/../outside/note.md" .claude/skills/user-owned-link

# --- (d) a real (non-symlink) user file in each of the four trees ---
mkdir -p .claude/output-styles
echo "user skill note"    > .claude/skills/user-real.md
echo "user agent note"    > .claude/agents/user-real.md
echo "user command note"  > .claude/commands/user-real.md
echo "user style note"    > .claude/output-styles/user-real.md

"$SCIAGENT" deactivate >/dev/null

# (a) recovered despite manifest amnesia.
if [[ -e .claude/skills/s_b || -L .claude/skills/s_b ]]; then
    echo "FAIL [$_TEST_NAME] s_b (dropped from manifest) survived teardown" >&2
    exit 1
fi

# (b) broken toolkit-owned link cleaned up.
if [[ -e .claude/skills/s_c || -L .claude/skills/s_c ]]; then
    echo "FAIL [$_TEST_NAME] broken toolkit-owned link s_c survived teardown" >&2
    exit 1
fi

# (c) user symlink outside the toolkit survives with its target unchanged.
assert_symlink .claude/skills/user-owned-link "user symlink outside toolkit must survive deactivate"
assert_eq "$(readlink .claude/skills/user-owned-link)" "$PWD/../outside/note.md" "target unchanged"

# (d) real user files in every tree survive untouched.
assert_file_exists .claude/skills/user-real.md         "real user file in skills/ survives"
assert_file_exists .claude/agents/user-real.md          "real user file in agents/ survives"
assert_file_exists .claude/commands/user-real.md        "real user file in commands/ survives"
assert_file_exists .claude/output-styles/user-real.md   "real user file in output-styles/ survives"
assert_eq "$(cat .claude/skills/user-real.md)"   "user skill note"   "skills/ content untouched"
assert_eq "$(cat .claude/agents/user-real.md)"   "user agent note"   "agents/ content untouched"
assert_eq "$(cat .claude/commands/user-real.md)" "user command note" "commands/ content untouched"

# No manifest, no other toolkit-owned mounts left.
[[ -e .sciagent/manifest.json ]] && { echo "FAIL [$_TEST_NAME] manifest still exists" >&2; exit 1; }
assert_symlink .claude/skills/user-owned-link "user link still the only thing left in .claude/skills"

# --- (e) the manifest FILE is gone entirely, not just an entry ---
# Distinct from (a): (a) drops an entry from a manifest that still exists, so a
# manifest-gated teardown still runs and merely misses one mount. Here the file
# itself is deleted — `rm -rf .sciagent`, a botched clean, a partially restored
# checkout. Teardown used to `return 0` on a missing manifest and strand EVERY
# mount permanently, since nothing else removes them. Target-based ownership
# needs no manifest, so this must recover fully.
cd "$TMPDIR_TEST"
mkdir project2 && cd project2
"$SCIAGENT" activate base >/dev/null
assert_symlink .claude/skills/s_a "sanity: mounted in project2"
mounted_before=$(find .claude .agents -mindepth 2 -maxdepth 2 -type l | wc -l)
[[ "$mounted_before" -gt 0 ]] || { echo "FAIL [$_TEST_NAME] fixture: nothing mounted in project2" >&2; exit 1; }

# A user file that must survive teardown even on this path.
echo "keep me" > .claude/skills/user-real.md

rm -rf .sciagent
[[ -e .sciagent/manifest.json ]] && { echo "FAIL [$_TEST_NAME] fixture: manifest not removed" >&2; exit 1; }

out=$("$SCIAGENT" deactivate 2>&1)
left=$(find .claude .agents -mindepth 2 -maxdepth 2 -type l 2>/dev/null | wc -l)
[[ "$left" -eq 0 ]] || {
    echo "FAIL [$_TEST_NAME] manifest-less teardown stranded $left of $mounted_before mounts" >&2
    printf '%s\n' "$out" >&2
    exit 1
}
assert_file_exists .claude/skills/user-real.md "user file survives manifest-less teardown"
assert_eq "$(cat .claude/skills/user-real.md)" "keep me" "user file content untouched"
# Reporting must reflect what happened: a removal count, not a bare
# "deactivated" (which would imply a manifest existed) and not "no active
# stack" (which would be a lie with N mounts on disk).
printf '%s\n' "$out" | grep -q 'orphaned mounts removed' || {
    echo "FAIL [$_TEST_NAME] expected an orphaned-mount count in output, got: $out" >&2
    exit 1
}

# --- (f) genuinely nothing mounted: must NOT claim to have deactivated ---
cd "$TMPDIR_TEST"
mkdir project3 && cd project3
out=$("$SCIAGENT" deactivate 2>&1)
assert_eq "$out" "no active stack" "empty project reports no active stack, not a phantom teardown"

pass
