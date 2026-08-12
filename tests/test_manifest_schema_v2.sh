#!/usr/bin/env bash
# tests/test_manifest_schema_v2.sh — .sciagent/manifest.json schema v2.
#
# v1 carried a `block_hash` field that `activate` wrote and NOTHING read: the
# drift guard reads the hash from the AGENTS.md BEGIN marker (block_hash_check),
# never from the manifest. It was removed. The interesting risk is not the
# removal but the fleet: 22 consumers hold v1 manifests written by an older
# toolkit, and those must keep working until their next activate rewrites them.
#
# Tests:
#   1. a freshly written manifest is v2, has no block_hash, and is valid JSON
#   2. manifest_finalize takes no argument (a stray one is not written anywhere)
#   3. a v1 manifest on disk still reads: status reports its stack
#   4. a v1 manifest still tears down: deactivate removes the links it lists
#   5. re-activating over a v1 manifest upgrades it to v2 in place
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
    [[ $# -gt 0 ]] && { echo "--- context ---" >&2; printf '%s\n' "$@" >&2; }
    exit 1
}

mkdir project && cd project
"$SCIAGENT" activate base >/dev/null

M=.sciagent/manifest.json

# ---------------------------------------------------------------------------
# 1. v2, no block_hash, valid JSON.
# ---------------------------------------------------------------------------
grep -q '"version": 2' "$M" || fail "manifest is not schema v2" "$(cat "$M")"
grep -q '"block_hash"'  "$M" && fail "manifest still carries block_hash" "$(cat "$M")"
if command -v jq >/dev/null 2>&1; then
    jq -e . "$M" >/dev/null 2>&1 || fail "manifest is not valid JSON" "$(cat "$M")"
    [[ "$(jq -r '.stack | join(" ")' "$M")" == "base" ]] \
        || fail "jq cannot read .stack back" "$(cat "$M")"
    [[ "$(jq -r '.block_hash // "absent"' "$M")" == "absent" ]] \
        || fail "block_hash resolvable via jq" "$(cat "$M")"
fi

# ---------------------------------------------------------------------------
# 2. manifest_finalize takes no argument. Call it directly with a stray one and
#    assert nothing of it reaches the file — this is what keeps a future caller
#    from quietly reintroducing the field through the old signature.
# ---------------------------------------------------------------------------
(
    # shellcheck source=/dev/null
    . "$FAKE/lib/sciagent/roles.sh"
    # shellcheck source=/dev/null
    . "$FAKE/lib/sciagent/symlinks.sh"
    manifest_begin "base"
    manifest_finalize "deadbeefdeadbeefdeadbeef" || exit 1
) || fail "manifest_finalize with a stray argument failed"
grep -q 'deadbeef' "$M" && fail "a stray manifest_finalize argument was written to the manifest" "$(cat "$M")"
grep -q '"version": 2' "$M" || fail "direct manifest_finalize did not write v2" "$(cat "$M")"

# ---------------------------------------------------------------------------
# 3-4. Backward tolerance: hand-plant a v1 manifest (the shape 22 consumers
#      hold) over a live activation and confirm both readers still work.
# ---------------------------------------------------------------------------
"$SCIAGENT" activate base >/dev/null      # restore the real symlink list
python3 - "$M" <<'PY' 2>/dev/null || {
import json, sys
p = sys.argv[1]
d = json.load(open(p))
d["version"] = 1
d["block_hash"] = "0123456789abcdef0123456789abcdef01234567"
json.dump(d, open(p, "w"), indent=2)
PY
    fail "could not synthesize a v1 manifest (python3 missing?)"
}
grep -q '"block_hash"' "$M" || fail "v1 fixture was not written"

out=$("$SCIAGENT" status 2>&1) || fail "status failed against a v1 manifest" "$out"
printf '%s\n' "$out" | grep -q 'base' \
    || fail "status could not read the stack out of a v1 manifest" "$out"

# The v1 manifest's symlinks list must still drive teardown.
[[ -L .claude/skills/s_a ]] || fail "fixture precondition: s_a not mounted"
out=$("$SCIAGENT" deactivate 2>&1) || fail "deactivate failed against a v1 manifest" "$out"
[[ -e .claude/skills/s_a ]] && fail "deactivate left a mount behind when driven by a v1 manifest"
[[ -e "$M" ]] && fail "deactivate left the v1 manifest in place"

# ---------------------------------------------------------------------------
# 5. Re-activating over a v1 manifest rewrites it as v2.
# ---------------------------------------------------------------------------
"$SCIAGENT" activate base >/dev/null
python3 - "$M" <<'PY'
import json, sys
p = sys.argv[1]
d = json.load(open(p))
d["version"] = 1
d["block_hash"] = "0123456789abcdef0123456789abcdef01234567"
json.dump(d, open(p, "w"), indent=2)
PY
"$SCIAGENT" activate base >/dev/null
grep -q '"version": 2' "$M" || fail "re-activate did not upgrade a v1 manifest" "$(cat "$M")"
grep -q '"block_hash"'  "$M" && fail "re-activate preserved the removed field" "$(cat "$M")"

pass
