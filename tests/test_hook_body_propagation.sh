#!/usr/bin/env bash
# tests/test_hook_body_propagation.sh
#
# A hook body fixed in the toolkit must actually REACH a project that was
# already provisioned. Before 2026-08-11 it never did: ensure_hooks wrote a body
# only `if [[ ! -f "$dst" ]]`, so the very first
# provisioning froze that file forever. Measured across the real fleet the same
# day: all 8 consumers carrying no_ephemeral.sh still held the af20086 version
# (12dafe8a3e39) and none had received de5719e's fix. A guardrail that cannot
# be updated in the field is only as good as the day it was installed.
#
# The fix cannot key off the .sciagent ownership records, which is what makes
# this subtle: NOT ONE consumer has a .sciagent/hook_state directory. Those
# records were introduced long after those projects were provisioned, so
# "refresh iff the record matches" is inert in exactly the population that
# needs repairing. Provenance is therefore established from CONTENT — the
# bodies are materialized with a plain `cp`, so a project's copy is
# byte-identical to whichever template version produced it, and
# templates/PROVENANCE.sha1 lists every version ever shipped.
#
# Cases (1)-(6) drive ownership_ensure_body directly against a SYNTHETIC
# toolkit, so every branch is exercised deterministically without depending on
# this repo's git history. Case (7) is the end-to-end wiring check against the
# REAL toolkit and the REAL manifest, reproducing the actual fleet defect.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir

# ---------------------------------------------------------------------------
# A synthetic toolkit: just enough for the ownership logic — one template and
# a provenance manifest we control completely.
# ---------------------------------------------------------------------------
FAKE_TK="$TMPDIR_TEST/fake-toolkit"
REL="project/_common/.claude/hooks/probe.sh.template"
mkdir -p "$FAKE_TK/templates/$(dirname "$REL")"

_write_tpl() { printf '#!/bin/bash\n# probe hook v%s\necho v%s\n' "$1" "$1" > "$FAKE_TK/templates/$REL"; }
_tpl_hash()  { sha1sum "$FAKE_TK/templates/$REL" | cut -d' ' -f1; }

# v1 is the "already shipped once" version; v2 is the current one.
_write_tpl 1; V1_BODY="$TMPDIR_TEST/v1.body"; cp "$FAKE_TK/templates/$REL" "$V1_BODY"; V1=$(_tpl_hash)
_write_tpl 2; V2=$(_tpl_hash)

cat > "$FAKE_TK/templates/PROVENANCE.sha1" <<EOF
# synthetic manifest
$V2  $REL
$V1  $REL
EOF

export SCIAGENT_TOOLKIT="$FAKE_TK"
# shellcheck source=/dev/null
. "$TOOLKIT_ROOT/lib/sciagent/block.sh"
# ownership.sh is where the discipline under test lives (it moved out of
# claude_settings.sh once the 02_analysis/helpers shims became its second
# consumer — an R helper is not a Claude artifact). claude_settings.sh is still
# sourced because the hook paths and their state-file layout are
# its half of the contract.
# shellcheck source=/dev/null
. "$TOOLKIT_ROOT/lib/sciagent/ownership.sh"
# shellcheck source=/dev/null
. "$TOOLKIT_ROOT/lib/sciagent/claude_settings.sh"

WORK="$TMPDIR_TEST/work"; mkdir -p "$WORK"; cd "$WORK"
SRC="$FAKE_TK/templates/$REL"
DST=".claude/hooks/probe.sh"
STATE=".sciagent/hook_state/probe.sh.sha1"
CEDED=".sciagent/hook_state/probe.sh.ceded"

_reset() { rm -rf .claude .sciagent; }
_hash()  { sha1sum "$1" | cut -d' ' -f1; }

# ---------------------------------------------------------------------------
# (1) Absent → created, recorded, executable. The original behaviour, retained.
# ---------------------------------------------------------------------------
_reset
out=$(ownership_ensure_body "$SRC" "$DST" "$REL" "$STATE" "$CEDED" 2>&1)
assert_file_exists "$DST" "(1) body created when absent"
assert_eq "$(_hash "$DST")" "$V2" "(1) created body is the current template"
assert_eq "$(cat "$STATE")" "$V2" "(1) ownership record written"
[[ -x "$DST" ]] || { echo "FAIL [$_TEST_NAME] (1) created body not executable" >&2; exit 1; }
case "$out" in *"wrote: $DST"*) : ;; *) echo "FAIL [$_TEST_NAME] (1) no 'wrote:' line: $out" >&2; exit 1 ;; esac

# ---------------------------------------------------------------------------
# (2) THE DEFECT. A stale body that is byte-identical to a version we shipped
# earlier, with NO ownership record — precisely the fleet's state — must be
# refreshed to the current template.
# ---------------------------------------------------------------------------
_reset
mkdir -p "$(dirname "$DST")"
cp "$V1_BODY" "$DST"
out=$(ownership_ensure_body "$SRC" "$DST" "$REL" "$STATE" "$CEDED" 2>&1)
assert_eq "$(_hash "$DST")" "$V2" "(2) stale-but-ours body REFRESHED to current template"
assert_eq "$(cat "$STATE")" "$V2" "(2) record updated to the refreshed content"
case "$out" in *"refreshed: $DST"*) : ;; *) echo "FAIL [$_TEST_NAME] (2) no 'refreshed:' line: $out" >&2; exit 1 ;; esac

# Mutation check: the refresh must depend on the manifest, not fire blindly.
# With the manifest hidden, the SAME stale body must NOT be overwritten.
_reset
mkdir -p "$(dirname "$DST")"
cp "$V1_BODY" "$DST"
mv "$FAKE_TK/templates/PROVENANCE.sha1" "$FAKE_TK/templates/PROVENANCE.sha1.hidden"
ownership_ensure_body "$SRC" "$DST" "$REL" "$STATE" "$CEDED" >/dev/null 2>&1
mv "$FAKE_TK/templates/PROVENANCE.sha1.hidden" "$FAKE_TK/templates/PROVENANCE.sha1"
assert_eq "$(_hash "$DST")" "$V1" "(2b) without the manifest the body is NOT overwritten (fix is non-vacuous)"

# ---------------------------------------------------------------------------
# (3) A USER-EDITED body — content matching no version we ever shipped — must
# survive untouched, warn once, and leave a ceded marker.
# ---------------------------------------------------------------------------
_reset
mkdir -p "$(dirname "$DST")"
printf '#!/bin/bash\n# I customized this myself\necho mine\n' > "$DST"
USER_HASH=$(_hash "$DST")
err=$(ownership_ensure_body "$SRC" "$DST" "$REL" "$STATE" "$CEDED" 2>&1 >/dev/null)
assert_eq "$(_hash "$DST")" "$USER_HASH" "(3) user-edited body NOT overwritten"
case "$err" in *"will not manage it"*) : ;; *) echo "FAIL [$_TEST_NAME] (3) no cede warning: $err" >&2; exit 1 ;; esac
assert_file_exists "$CEDED" "(3) ceded marker written"
assert_eq "$(cat "$CEDED")" "$USER_HASH" "(3) marker records the ceded content"
[[ -f "$STATE" ]] && { echo "FAIL [$_TEST_NAME] (3) ownership record should have been dropped" >&2; exit 1; }

# ---------------------------------------------------------------------------
# (4) Ceding is PERMANENT and quiet: a second run over the same content must
# emit nothing. Re-warning on every link is the noise the marker prevents.
# ---------------------------------------------------------------------------
err2=$(ownership_ensure_body "$SRC" "$DST" "$REL" "$STATE" "$CEDED" 2>&1 >/dev/null)
assert_eq "$err2" "" "(4) second run over a ceded body is silent"
assert_eq "$(_hash "$DST")" "$USER_HASH" "(4) still untouched"

# ---------------------------------------------------------------------------
# (5) ...but if the user REPLACES it with something else, that is a new
# situation and must be reported again.
# ---------------------------------------------------------------------------
printf '#!/bin/bash\n# different customization\necho other\n' > "$DST"
OTHER_HASH=$(_hash "$DST")
err3=$(ownership_ensure_body "$SRC" "$DST" "$REL" "$STATE" "$CEDED" 2>&1 >/dev/null)
case "$err3" in *"will not manage it"*) : ;; *) echo "FAIL [$_TEST_NAME] (5) re-cede not reported: $err3" >&2; exit 1 ;; esac
assert_eq "$(cat "$CEDED")" "$OTHER_HASH" "(5) marker follows the new content"
assert_eq "$(_hash "$DST")" "$OTHER_HASH" "(5) new customization untouched"

# ---------------------------------------------------------------------------
# (6) Already current, no record (a body predating the ownership records) →
# adopted silently for future refreshes. No write, no message.
# ---------------------------------------------------------------------------
_reset
mkdir -p "$(dirname "$DST")"
cp "$SRC" "$DST"
out6=$(ownership_ensure_body "$SRC" "$DST" "$REL" "$STATE" "$CEDED" 2>&1)
assert_eq "$out6" "" "(6) up-to-date body produces no output"
assert_eq "$(cat "$STATE")" "$V2" "(6) adopted: record written for a pre-existing current body"

# ===========================================================================
# (7) END-TO-END against the REAL toolkit and the REAL manifest: reproduce the
# fleet defect exactly. Seed a project with the OLDEST shipped version of
# no_ephemeral.sh and no ownership record, then run a real `link`.
# ===========================================================================
cd "$TMPDIR_TEST"
REAL_REL="project/_common/.claude/hooks/no_ephemeral.sh.template"
REAL_TPL="$TOOLKIT_ROOT/templates/$REAL_REL"
OLD_SHA=$(git -C "$TOOLKIT_ROOT" log --all --format=%H -- "templates/$REAL_REL" 2>/dev/null | tail -1)

if [[ -z "$OLD_SHA" ]] || ! git -C "$TOOLKIT_ROOT" show "$OLD_SHA:templates/$REAL_REL" >/dev/null 2>&1; then
    # Loudly, not silently: a shallow clone would otherwise make this vacuous.
    echo "SKIP [$_TEST_NAME] (7) no git history for $REAL_REL — end-to-end case not run" >&2
else
    git -C "$TOOLKIT_ROOT" show "$OLD_SHA:templates/$REAL_REL" > "$TMPDIR_TEST/old_hook.body"
    OLD_HASH=$(_hash "$TMPDIR_TEST/old_hook.body")
    CUR_HASH=$(_hash "$REAL_TPL")

    if [[ "$OLD_HASH" == "$CUR_HASH" ]]; then
        echo "SKIP [$_TEST_NAME] (7) template has only ever had one version" >&2
    else
        mkdir -p e2e/.claude/hooks && cd e2e
        printf '# Project AGENTS.md\n\nUser-owned content.\n' > AGENTS.md
        cp "$TMPDIR_TEST/old_hook.body" .claude/hooks/no_ephemeral.sh
        # No .sciagent/hook_state — exactly the fleet's state.
        [[ -e .sciagent/hook_state ]] && { echo "FAIL [$_TEST_NAME] (7) fixture seeded a record" >&2; exit 1; }

        SCIAGENT_TOOLKIT="$TOOLKIT_ROOT" "$TOOLKIT_ROOT/bin/sciagent" link >/dev/null 2>&1

        assert_eq "$(_hash .claude/hooks/no_ephemeral.sh)" "$CUR_HASH" \
            "(7) real link refreshed the stale hook to the current template"
        [[ -x .claude/hooks/no_ephemeral.sh ]] || {
            echo "FAIL [$_TEST_NAME] (7) refreshed hook not executable" >&2; exit 1; }

        # And the refreshed body must still accept the pre-rename layout: this
        # propagates into 8 repos that have NOT done scripts/ -> stages/ yet.
        assert_grep 'stages|scripts' .claude/hooks/no_ephemeral.sh \
            "(7) refreshed hook still accepts BOTH stage dir names during the migration"

        cd ..
    fi
fi

pass
