#!/usr/bin/env bash
# tests/test_craft_verb.sh
# `scio craft` renders the SCIO:CRAFT block with no role and no mounts:
#   - works in a repo that has never been linked
#   - mounts nothing: no .claude/, .agents/, .scio/, no ROLES block
#   - re-render is byte-identical (no-op)
#   - --project-dir targets another directory
#   - a hand-edit inside the markers is refused, and --force overrides
#   - a valid SCIAGENT:CRAFT block is rewritten in place under the same guard
#   - user-owned bytes outside the markers are preserved
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"

cat > "$FAKE/craft.yaml" <<'EOF'
version: 1

floors:
  figure_base_size: 16

body: |
  # Craft standards
  - Figures: base >= {{figure_base_size}}pt.
EOF

export SCIO_TOOLKIT="$FAKE"
SCIO="$FAKE/bin/scio"

# ── 1. Unbound repo gains the block. ───────────────────────────────────────
mkdir umbrella && cd umbrella
cat > AGENTS.md <<'EOF'
# Umbrella AGENTS.md

Coordination repo. No analysis skills belong here.
EOF

out=$("$SCIO" craft) || { echo "FAIL [$_TEST_NAME] craft exited non-zero: $out" >&2; exit 1; }
assert_grep 'BEGIN SCIO:CRAFT' AGENTS.md "CRAFT block written"
assert_grep 'base >= 16pt'         AGENTS.md "floor token substituted"
assert_grep 'Coordination repo'    AGENTS.md "user content preserved"

# ── 2. Mounted nothing. ─────────────────────────────────────────────────────
for d in .claude .agents .scio; do
    if [[ -e "$d" ]]; then
        echo "FAIL [$_TEST_NAME] craft created $d — it must mount nothing" >&2
        exit 1
    fi
done
if grep -qF 'SCIO:ROLES' AGENTS.md; then
    echo "FAIL [$_TEST_NAME] craft wrote a ROLES block" >&2
    exit 1
fi

# ── 3. Re-render is a byte-identical no-op. ─────────────────────────────────
cp AGENTS.md AGENTS.md.first
out=$("$SCIO" craft)
assert_file_eq AGENTS.md AGENTS.md.first "re-render leaves AGENTS.md byte-identical"
printf '%s\n' "$out" | grep -q 'already current' \
    || { echo "FAIL [$_TEST_NAME] expected 'already current' on no-op re-render, got: $out" >&2; exit 1; }

# ── 4. Drift guard: hand-edit inside the markers is refused. ────────────────
sed -i 's/Craft standards/Craft STANDARDS/' AGENTS.md
cp AGENTS.md AGENTS.md.drifted
set +e
out=$("$SCIO" craft 2>&1); rc=$?
set -e
[[ "$rc" -eq 1 ]] || { echo "FAIL [$_TEST_NAME] drifted block: expected exit 1 got $rc" >&2; exit 1; }
printf '%s\n' "$out" | grep -q 'drifted' \
    || { echo "FAIL [$_TEST_NAME] expected a drift message, got: $out" >&2; exit 1; }
assert_file_eq AGENTS.md AGENTS.md.drifted "refusal left the file untouched"

# --force overwrites the drift and restores the canonical render.
"$SCIO" craft --force >/dev/null
assert_file_eq AGENTS.md AGENTS.md.first "--force restores the canonical block"

# ── 5. A legacy block migrates in place; legacy drift still refuses. ────────
cd "$TMPDIR_TEST"
mkdir legacy legacy-drift
sed 's/SCIO:CRAFT/SCIAGENT:CRAFT/g' umbrella/AGENTS.md.first > legacy/AGENTS.md
printf '\n# Content after the managed block.\n' >> legacy/AGENTS.md
chmod 664 legacy/AGENTS.md
cp legacy/AGENTS.md legacy.expected
sed -i 's/SCIAGENT:CRAFT/SCIO:CRAFT/g' legacy.expected

legacy_out=$("$SCIO" craft --project-dir "$TMPDIR_TEST/legacy" 2>&1) || {
    echo "FAIL [$_TEST_NAME] legacy CRAFT migration failed: $legacy_out" >&2
    exit 1
}
assert_file_eq legacy/AGENTS.md legacy.expected \
    "legacy rewrite changes only the marker prefix in place"
assert_eq "$(stat -c '%a' legacy/AGENTS.md)" "664" "legacy rewrite preserves mode"
assert_eq "$(grep -c 'BEGIN SCIO:CRAFT' legacy/AGENTS.md)" "1" \
    "legacy rewrite leaves exactly one CRAFT block"
if grep -q 'SCIAGENT:CRAFT' legacy/AGENTS.md; then
    echo "FAIL [$_TEST_NAME] legacy marker survived migration" >&2
    exit 1
fi

cp legacy.expected legacy-drift/AGENTS.md
sed -i 's/SCIO:CRAFT/SCIAGENT:CRAFT/g; s/Craft standards/Craft STANDARDS/' \
    legacy-drift/AGENTS.md
cp legacy-drift/AGENTS.md legacy-drift/AGENTS.md.before
set +e
legacy_drift_out=$("$SCIO" craft --project-dir "$TMPDIR_TEST/legacy-drift" 2>&1); legacy_drift_rc=$?
set -e
[[ "$legacy_drift_rc" -eq 1 ]] || {
    echo "FAIL [$_TEST_NAME] drifted legacy block: expected exit 1 got $legacy_drift_rc" >&2
    exit 1
}
printf '%s\n' "$legacy_drift_out" | grep -q 'drifted' \
    || { echo "FAIL [$_TEST_NAME] legacy drift produced no refusal: $legacy_drift_out" >&2; exit 1; }
assert_file_eq legacy-drift/AGENTS.md legacy-drift/AGENTS.md.before \
    "legacy drift refusal leaves the file untouched"
assert_grep 'BEGIN SCIAGENT:CRAFT' legacy-drift/AGENTS.md \
    "legacy drift refusal does not migrate the marker"

# ── 6. --project-dir targets elsewhere; AGENTS.md is created if absent. ─────
cd "$TMPDIR_TEST"
mkdir elsewhere
"$SCIO" craft --project-dir "$TMPDIR_TEST/elsewhere" >/dev/null
assert_file_exists "$TMPDIR_TEST/elsewhere/AGENTS.md" "craft created AGENTS.md in --project-dir"
assert_grep 'BEGIN SCIO:CRAFT' "$TMPDIR_TEST/elsewhere/AGENTS.md" "block written to --project-dir"

# A nonexistent --project-dir is an error, not a silent mkdir.
set +e
"$SCIO" craft --project-dir "$TMPDIR_TEST/nope" >/dev/null 2>&1; rc=$?
set -e
[[ "$rc" -eq 1 ]] || { echo "FAIL [$_TEST_NAME] missing --project-dir: expected exit 1 got $rc" >&2; exit 1; }

# ── 7. A toolkit with no craft.yaml fails loudly (the verb was asked for). ──
rm "$FAKE/craft.yaml"
set +e
out=$("$SCIO" craft --project-dir "$TMPDIR_TEST/elsewhere" 2>&1); rc=$?
set -e
[[ "$rc" -eq 1 ]] || { echo "FAIL [$_TEST_NAME] absent craft.yaml: expected exit 1 got $rc" >&2; exit 1; }
printf '%s\n' "$out" | grep -q 'no craft.yaml' \
    || { echo "FAIL [$_TEST_NAME] expected a no-craft.yaml message, got: $out" >&2; exit 1; }

pass
