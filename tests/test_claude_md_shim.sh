#!/usr/bin/env bash
# tests/test_claude_md_shim.sh — activate.sh's ensure_claude_md_shim.
#
# Guarantees the project CLAUDE.md carries the `@AGENTS.md` import so Claude Code
# (which doesn't read AGENTS.md natively) picks up the canonical context:
#   - absent CLAUDE.md → created containing exactly `@AGENTS.md`;
#   - present WITH the import → exact no-op (no duplicate line);
#   - present WITHOUT the import → import prepended, existing content preserved.
#
# Each case also asserts the OWNERSHIP RECORD (.sciagent/claude_md.state) that
# `deactivate` later reads to reverse exactly what we did — without those
# assertions this test was vacuous: it passed while `sciagent_sha1_file` was
# undefined (`command not found` twice) because the hash lands in a state file
# nothing here looked at, so an empty hash was indistinguishable from a correct
# one. A wrong/empty hash means teardown silently declines to reverse and warns
# instead, which is the failure mode the assertions below now catch.
#
# Sourced directly (smallest surface); hermetic scratch cwd via setup_tmpdir.
# block.sh comes first: it is the single home for the hashing primitive (the 5.3
# invariant) and activate.sh calls sciagent_sha1_file/sciagent_sha1_stream from
# it, per the load graph in bin/sciagent.
set -u
. "$(dirname "$0")/_lib.sh"
. "$TOOLKIT_ROOT/lib/sciagent/block.sh"
. "$TOOLKIT_ROOT/lib/sciagent/activate.sh"

setup_tmpdir

STATE=".sciagent/claude_md.state"

# assert_state <expected-tag> <expected-hash> <case>
# Both lines matter: teardown dispatches on the tag and refuses to act unless
# the hash still matches, so an empty second line fails safe into "cede and
# warn" — a silent regression at the unit level.
assert_state() {
    local tag="$1" hash="$2" label="$3"
    # Guard the guard: the expected hash must itself be a 40-hex sha1, never the
    # empty string an undefined hashing helper leaves behind — otherwise
    # "" == "" would pass and this test would be vacuous again.
    case "$hash" in
        [0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f]) ;;
        *) echo "FAIL [$_TEST_NAME] $label: expected hash is not a sha1: '$hash'" >&2; exit 1 ;;
    esac
    assert_file_exists "$STATE" "$label: ownership record written"
    assert_eq "$(sed -n '1p' "$STATE")" "$tag"  "$label: state tag"
    assert_eq "$(sed -n '2p' "$STATE")" "$hash" "$label: state hash"
}

# ----- Case 1: absent → created with @AGENTS.md ----------------------------
mkdir case1 && cd case1
ensure_claude_md_shim >/dev/null
assert_file_exists CLAUDE.md "case1: CLAUDE.md created"
assert_grep '^@AGENTS\.md$' CLAUDE.md "case1: import present"
assert_eq "$(cat CLAUDE.md)" "@AGENTS.md" "case1: file is exactly the import line"
# "created" → the hash is of the file AS WRITTEN, so teardown can delete it
# only while it still hashes to that.
assert_state created "$(sciagent_sha1_file CLAUDE.md)" case1
cd ..

# ----- Case 2: present with import → exact no-op ----------------------------
mkdir case2 && cd case2
printf '@AGENTS.md\n' > CLAUDE.md
cp CLAUDE.md orig
ensure_claude_md_shim >/dev/null
assert_file_eq CLAUDE.md orig "case2: file unchanged when import already present"
assert_eq "$(grep -c '^@AGENTS\.md$' CLAUDE.md)" "1" "case2: no duplicate import line"
# No-op → NO ownership record: there is nothing for teardown to reverse, and a
# spurious "created" record here would make deactivate delete a user's file.
assert_eq "$([[ -e "$STATE" ]] && echo present || echo absent)" absent \
    "case2: no ownership record for the no-op path"
cd ..

# ----- Case 3: present without import → prepend, preserve content ----------
mkdir case3 && cd case3
cat > CLAUDE.md <<'EOF'
# Project notes

Some pre-existing user content.
EOF
# "prepended" records the hash of the content BEFORE the header goes on, so
# snapshot it here — teardown strips the header and compares the remainder.
orig_hash=$(sciagent_sha1_file CLAUDE.md)
ensure_claude_md_shim >/dev/null
assert_eq "$(head -1 CLAUDE.md)" "@AGENTS.md" "case3: import prepended at top"
assert_eq "$(grep -c '^@AGENTS\.md$' CLAUDE.md)" "1" "case3: single import line"
assert_grep '^# Project notes$' CLAUDE.md "case3: existing heading preserved"
assert_grep '^Some pre-existing user content\.$' CLAUDE.md "case3: existing prose preserved"
assert_state prepended "$orig_hash" case3
# The recorded hash must round-trip: stripping the exact "@AGENTS.md\n\n"
# header must reproduce it, or teardown warns and leaves the header behind
# forever. This is the assertion the whole ownership record exists to support.
prefix=$'@AGENTS.md\n\n'
assert_eq "$(tail -c "+$((${#prefix} + 1))" CLAUDE.md | sciagent_sha1_stream)" \
    "$orig_hash" "case3: stripping the header reproduces the recorded hash"
cd ..

pass
