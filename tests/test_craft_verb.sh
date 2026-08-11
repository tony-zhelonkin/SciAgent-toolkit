#!/usr/bin/env bash
# tests/test_craft_verb.sh
# `sciagent craft` renders the SCIAGENT:CRAFT block with no role and no mounts:
#   - works in a repo that has never been activated (the umbrella case)
#   - mounts nothing: no .claude/, .agents/, .sciagent/, no ROLES block
#   - re-render is byte-identical (no-op)
#   - --project-dir targets another directory
#   - a hand-edit inside the markers is refused, and --force overrides
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

export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

# ── 1. Never-activated repo gains the block. ────────────────────────────────
mkdir umbrella && cd umbrella
cat > AGENTS.md <<'EOF'
# Umbrella AGENTS.md

Coordination repo. No analysis skills belong here.
EOF

out=$("$SCIAGENT" craft) || { echo "FAIL [$_TEST_NAME] craft exited non-zero: $out" >&2; exit 1; }
assert_grep 'BEGIN SCIAGENT:CRAFT' AGENTS.md "CRAFT block written"
assert_grep 'base >= 16pt'         AGENTS.md "floor token substituted"
assert_grep 'Coordination repo'    AGENTS.md "user content preserved"

# ── 2. Mounted nothing. ─────────────────────────────────────────────────────
for d in .claude .agents .sciagent; do
    if [[ -e "$d" ]]; then
        echo "FAIL [$_TEST_NAME] craft created $d — it must mount nothing" >&2
        exit 1
    fi
done
if grep -qF 'SCIAGENT:ROLES' AGENTS.md; then
    echo "FAIL [$_TEST_NAME] craft wrote a ROLES block" >&2
    exit 1
fi

# ── 3. Re-render is a byte-identical no-op. ─────────────────────────────────
cp AGENTS.md AGENTS.md.first
out=$("$SCIAGENT" craft)
assert_file_eq AGENTS.md AGENTS.md.first "re-render leaves AGENTS.md byte-identical"
printf '%s\n' "$out" | grep -q 'already current' \
    || { echo "FAIL [$_TEST_NAME] expected 'already current' on no-op re-render, got: $out" >&2; exit 1; }

# ── 4. Drift guard: hand-edit inside the markers is refused. ────────────────
sed -i 's/Craft standards/Craft STANDARDS/' AGENTS.md
cp AGENTS.md AGENTS.md.drifted
set +e
out=$("$SCIAGENT" craft 2>&1); rc=$?
set -e
[[ "$rc" -eq 1 ]] || { echo "FAIL [$_TEST_NAME] drifted block: expected exit 1 got $rc" >&2; exit 1; }
printf '%s\n' "$out" | grep -q 'drifted' \
    || { echo "FAIL [$_TEST_NAME] expected a drift message, got: $out" >&2; exit 1; }
assert_file_eq AGENTS.md AGENTS.md.drifted "refusal left the file untouched"

# --force overwrites the drift and restores the canonical render.
"$SCIAGENT" craft --force >/dev/null
assert_file_eq AGENTS.md AGENTS.md.first "--force restores the canonical block"

# ── 5. --project-dir targets elsewhere; AGENTS.md is created if absent. ─────
cd "$TMPDIR_TEST"
mkdir elsewhere
"$SCIAGENT" craft --project-dir "$TMPDIR_TEST/elsewhere" >/dev/null
assert_file_exists "$TMPDIR_TEST/elsewhere/AGENTS.md" "craft created AGENTS.md in --project-dir"
assert_grep 'BEGIN SCIAGENT:CRAFT' "$TMPDIR_TEST/elsewhere/AGENTS.md" "block written to --project-dir"

# A nonexistent --project-dir is an error, not a silent mkdir.
set +e
"$SCIAGENT" craft --project-dir "$TMPDIR_TEST/nope" >/dev/null 2>&1; rc=$?
set -e
[[ "$rc" -eq 1 ]] || { echo "FAIL [$_TEST_NAME] missing --project-dir: expected exit 1 got $rc" >&2; exit 1; }

# ── 6. A toolkit with no craft.yaml fails loudly (the verb was asked for). ──
rm "$FAKE/craft.yaml"
set +e
out=$("$SCIAGENT" craft --project-dir "$TMPDIR_TEST/elsewhere" 2>&1); rc=$?
set -e
[[ "$rc" -eq 1 ]] || { echo "FAIL [$_TEST_NAME] absent craft.yaml: expected exit 1 got $rc" >&2; exit 1; }
printf '%s\n' "$out" | grep -q 'no craft.yaml' \
    || { echo "FAIL [$_TEST_NAME] expected a no-craft.yaml message, got: $out" >&2; exit 1; }

pass
