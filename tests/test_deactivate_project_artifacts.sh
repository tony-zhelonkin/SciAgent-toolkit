#!/usr/bin/env bash
# tests/test_deactivate_project_artifacts.sh
#
# `sciagent deactivate` must be a true inverse of `sciagent activate` for the
# project-level artifacts activate ensures unconditionally (independent of
# which role is active): .claude/statusline.sh, .claude/hooks/*.sh, the
# settings.json key-backfill, and the CLAUDE.md @AGENTS.md import shim.
# Before this fix, deactivate reversed none of these — a repo left with
# PreToolUse/Stop hooks registered AND materialized, permanently, with no
# verb to undo it.
#
# Uses the REAL toolkit (like test_claude_settings_lifecycle.sh /
# test_enforcement_hooks.sh), not the fake fixture toolkit — the fake
# toolkit ships no templates/project/_common/, so ensure_statusline/
# ensure_hooks/ensure_project_defaults silently no-op against it and the
# defect this test guards against would never be exercised.
set -u
. "$(dirname "$0")/_lib.sh"

SCIAGENT="$TOOLKIT_ROOT/bin/sciagent"
export SCIAGENT_TOOLKIT="$TOOLKIT_ROOT"

# _seed_project — AGENTS.md content activate's ROLES block gets embedded
# into; matches test_deactivate_full.sh's convention. AGENTS.md's own
# create-vs-preexist ownership is a separate, pre-existing, out-of-scope
# question (block_write creates the file if absent with no ownership record
# of the FILE itself, only of the markers) — every case here pre-seeds
# AGENTS.md so it never enters play.
_seed_agents_md() {
    cat > AGENTS.md <<'EOF'
# Project AGENTS.md

User-owned content.
EOF
}

setup_tmpdir

# ===========================================================================
# Case 1: virgin project — activate/deactivate leaves NO .claude/, NO
# .sciagent/, and NO CLAUDE.md residue (CLAUDE.md did not exist beforehand).
# ===========================================================================
mkdir case1 && cd case1
_seed_agents_md
cp AGENTS.md AGENTS.md.orig

"$SCIAGENT" activate base >/tmp/.c1-activate.out 2>&1
assert_file_exists .claude/statusline.sh   "case1: statusline materialized by activate"
assert_file_exists .claude/settings.json   "case1: settings.json materialized by activate"
assert_file_exists .claude/hooks/no_ephemeral.sh "case1: hook body materialized by activate"
assert_file_exists CLAUDE.md               "case1: CLAUDE.md created by activate"

"$SCIAGENT" deactivate >/tmp/.c1-deactivate.out 2>&1

[[ -e .claude ]]    && { echo "FAIL [$_TEST_NAME] case1: .claude/ residue: $(find .claude)"; exit 1; }
[[ -e .sciagent ]]  && { echo "FAIL [$_TEST_NAME] case1: .sciagent/ residue: $(find .sciagent)"; exit 1; }
[[ -e CLAUDE.md ]]  && { echo "FAIL [$_TEST_NAME] case1: CLAUDE.md residue (we created it, should be gone)"; exit 1; }
assert_file_eq AGENTS.md AGENTS.md.orig "case1: AGENTS.md restored byte-for-byte"
cd ..

# ===========================================================================
# Case 2: pre-existing CLAUDE.md with user prose (no @AGENTS.md import) —
# activate prepends the import; deactivate must strip EXACTLY that header,
# leaving the user's prose byte-identical.
# ===========================================================================
mkdir case2 && cd case2
_seed_agents_md
cat > CLAUDE.md <<'EOF'
# My hand-written project notes

Some prose the user wrote. Multiple lines.
Do not touch this.
EOF
cp CLAUDE.md CLAUDE.md.orig

"$SCIAGENT" activate base >/dev/null 2>&1
assert_eq "$(head -1 CLAUDE.md)" "@AGENTS.md" "case2: import prepended by activate"

"$SCIAGENT" deactivate >/dev/null 2>&1

assert_file_exists CLAUDE.md "case2: pre-existing CLAUDE.md must survive deactivate"
assert_file_eq CLAUDE.md CLAUDE.md.orig "case2: user prose byte-identical after activate+deactivate"
cd ..

# ===========================================================================
# Case 3: pre-existing settings.json with the user's own keys — activate
# backfills the gold defaults (including hook registration); deactivate must
# retain EXACTLY the user's original keys and remove the hook registration
# (not just the hook files on disk — the registration itself).
# ===========================================================================
mkdir case3 && cd case3
_seed_agents_md
mkdir -p .claude
cat > .claude/settings.json <<'EOF'
{
  "myCustomKey": "keepme",
  "permissions": {
    "allow": ["Bash(git:*)"]
  }
}
EOF

"$SCIAGENT" activate base >/dev/null 2>&1
assert_grep 'PreToolUse' .claude/settings.json "case3: activate registers PreToolUse hook"

"$SCIAGENT" deactivate >/dev/null 2>&1

assert_file_exists .claude/settings.json "case3: pre-existing settings.json must survive"
if command -v jq >/dev/null 2>&1; then
    keys=$(jq -r 'keys | sort | join(",")' .claude/settings.json)
    assert_eq "$keys" "myCustomKey,permissions" "case3: exactly the user's own keys remain"
    val=$(jq -r '.myCustomKey' .claude/settings.json)
    assert_eq "$val" "keepme" "case3: user's custom key value untouched"
fi
if grep -q 'PreToolUse\|no_ephemeral\|caption_sweep' .claude/settings.json; then
    echo "FAIL [$_TEST_NAME] case3: hook registration still present after deactivate" >&2
    cat .claude/settings.json >&2
    exit 1
fi
cd ..

# ===========================================================================
# Case 4: a user-authored hook in .claude/hooks/ survives deactivate; the
# hooks WE materialized (unedited) are removed; .claude/hooks/ itself
# survives because real content (the user's hook) remains.
# ===========================================================================
mkdir case4 && cd case4
_seed_agents_md
mkdir -p .claude/hooks
printf '#!/bin/bash\necho "user hook"\n' > .claude/hooks/my_custom_hook.sh
chmod +x .claude/hooks/my_custom_hook.sh
cp .claude/hooks/my_custom_hook.sh /tmp/.c4-hook.orig

"$SCIAGENT" activate base >/dev/null 2>&1
assert_file_exists .claude/hooks/no_ephemeral.sh "case4: activate materializes our hooks alongside the user's"

"$SCIAGENT" deactivate >/dev/null 2>&1

assert_file_exists .claude/hooks/my_custom_hook.sh "case4: user's own hook survives deactivate"
assert_file_eq .claude/hooks/my_custom_hook.sh /tmp/.c4-hook.orig "case4: user's hook byte-identical"
[[ -e .claude/hooks/no_ephemeral.sh ]] && { echo "FAIL [$_TEST_NAME] case4: our hook not removed"; exit 1; }
[[ -e .claude/hooks/caption_sweep.sh ]] && { echo "FAIL [$_TEST_NAME] case4: our hook not removed"; exit 1; }
rm -f /tmp/.c4-hook.orig
cd ..

# ===========================================================================
# Case 5: the user edits a value WE backfilled into settings.json before
# deactivating — that key must survive with its user-set value, and a
# warning must be emitted; the OTHER backfilled keys we didn't touch are
# still removed normally (per-key granularity, not all-or-nothing).
# ===========================================================================
mkdir case5 && cd case5
_seed_agents_md
mkdir -p .claude
cat > .claude/settings.json <<'EOF'
{
  "permissions": {"allow": ["Bash(git:*)"]}
}
EOF

command -v jq >/dev/null 2>&1 || { echo "SKIP [$_TEST_NAME] case5: jq unavailable"; cd ..; }
if command -v jq >/dev/null 2>&1; then
    "$SCIAGENT" activate base >/dev/null 2>&1
    jq '.effortLevel = "low"' .claude/settings.json > /tmp/.c5.json && mv /tmp/.c5.json .claude/settings.json

    out=$("$SCIAGENT" deactivate 2>&1)
    printf '%s\n' "$out" | grep -q "effortLevel.*modified" || {
        echo "FAIL [$_TEST_NAME] case5: expected a warning about the edited effortLevel key" >&2
        printf '%s\n' "$out" >&2
        exit 1
    }
    val=$(jq -r '.effortLevel' .claude/settings.json)
    assert_eq "$val" "low" "case5: user-edited value not clobbered"
    # editorMode was backfilled and never touched by the user — must be gone.
    if jq -e 'has("editorMode")' .claude/settings.json >/dev/null 2>&1; then
        echo "FAIL [$_TEST_NAME] case5: untouched backfilled key 'editorMode' not removed" >&2
        cat .claude/settings.json >&2
        exit 1
    fi
    cd ..
fi

# ===========================================================================
# Case 6: idempotent — deactivating an already-deactivated (or never
# activated) project is a clean no-op: rc 0, no stderr, no filesystem change.
# ===========================================================================
mkdir case6 && cd case6
_seed_agents_md
"$SCIAGENT" activate base >/dev/null 2>&1
"$SCIAGENT" deactivate >/dev/null 2>&1

find . -type f | sort > /tmp/.c6-before.txt
out2=$("$SCIAGENT" deactivate 2>&1)
rc2=$?
find . -type f | sort > /tmp/.c6-after.txt

assert_eq "$rc2" "0" "case6: second deactivate exits 0"
assert_eq "$out2" "no active stack" "case6: second deactivate reports no active stack"
diff -u /tmp/.c6-before.txt /tmp/.c6-after.txt >/tmp/.c6-diff.txt || {
    echo "FAIL [$_TEST_NAME] case6: second deactivate changed the filesystem" >&2
    cat /tmp/.c6-diff.txt >&2
    exit 1
}
rm -f /tmp/.c6-before.txt /tmp/.c6-after.txt /tmp/.c6-diff.txt

# Never-activated project: deactivate must not create anything.
mkdir ../case6b && cd ../case6b
out3=$("$SCIAGENT" deactivate 2>&1)
rc3=$?
assert_eq "$rc3" "0" "case6b: deactivate on never-activated project exits 0"
[[ -e .claude || -e .sciagent ]] && { echo "FAIL [$_TEST_NAME] case6b: deactivate created state on a virgin dir"; exit 1; }
cd ..

# ===========================================================================
# Case 7: drift — the user edits a file we CREATED (not backfilled-into)
# before deactivating: statusline.sh and a hook body. Both must survive with
# a warning, exactly like the settings.json per-key case above.
# ===========================================================================
mkdir case7 && cd case7
_seed_agents_md
"$SCIAGENT" activate base >/dev/null 2>&1
echo "# user addition" >> .claude/statusline.sh
echo "# user addition" >> .claude/hooks/no_ephemeral.sh

out4=$("$SCIAGENT" deactivate 2>&1)
printf '%s\n' "$out4" | grep -q "statusline.sh was modified" || {
    echo "FAIL [$_TEST_NAME] case7: expected statusline drift warning" >&2
    printf '%s\n' "$out4" >&2
    exit 1
}
printf '%s\n' "$out4" | grep -q "no_ephemeral.sh was modified" || {
    echo "FAIL [$_TEST_NAME] case7: expected hook drift warning" >&2
    printf '%s\n' "$out4" >&2
    exit 1
}
assert_file_exists .claude/statusline.sh "case7: edited statusline survives"
assert_file_exists .claude/hooks/no_ephemeral.sh "case7: edited hook survives"
assert_grep '# user addition' .claude/statusline.sh "case7: statusline edit preserved"
assert_grep '# user addition' .claude/hooks/no_ephemeral.sh "case7: hook edit preserved"
# The untouched hook (caption_sweep.sh) must still be removed normally.
[[ -e .claude/hooks/caption_sweep.sh ]] && { echo "FAIL [$_TEST_NAME] case7: untouched hook not removed"; exit 1; }
cd ..

pass
