#!/usr/bin/env bash
# tests/test_claude_md_shim.sh — activate.sh's ensure_claude_md_shim.
#
# Guarantees the project CLAUDE.md carries the `@AGENTS.md` import so Claude Code
# (which doesn't read AGENTS.md natively) picks up the canonical context:
#   - absent CLAUDE.md → created containing exactly `@AGENTS.md`;
#   - present WITH the import → exact no-op (no duplicate line);
#   - present WITHOUT the import → import prepended, existing content preserved.
#
# Sourced directly (smallest surface); hermetic scratch cwd via setup_tmpdir.
set -u
. "$(dirname "$0")/_lib.sh"
. "$TOOLKIT_ROOT/lib/sciagent/activate.sh"

setup_tmpdir

# ----- Case 1: absent → created with @AGENTS.md ----------------------------
mkdir case1 && cd case1
ensure_claude_md_shim >/dev/null
assert_file_exists CLAUDE.md "case1: CLAUDE.md created"
assert_grep '^@AGENTS\.md$' CLAUDE.md "case1: import present"
assert_eq "$(cat CLAUDE.md)" "@AGENTS.md" "case1: file is exactly the import line"
cd ..

# ----- Case 2: present with import → exact no-op ----------------------------
mkdir case2 && cd case2
printf '@AGENTS.md\n' > CLAUDE.md
cp CLAUDE.md orig
ensure_claude_md_shim >/dev/null
assert_file_eq CLAUDE.md orig "case2: file unchanged when import already present"
assert_eq "$(grep -c '^@AGENTS\.md$' CLAUDE.md)" "1" "case2: no duplicate import line"
cd ..

# ----- Case 3: present without import → prepend, preserve content ----------
mkdir case3 && cd case3
cat > CLAUDE.md <<'EOF'
# Project notes

Some pre-existing user content.
EOF
ensure_claude_md_shim >/dev/null
assert_eq "$(head -1 CLAUDE.md)" "@AGENTS.md" "case3: import prepended at top"
assert_eq "$(grep -c '^@AGENTS\.md$' CLAUDE.md)" "1" "case3: single import line"
assert_grep '^# Project notes$' CLAUDE.md "case3: existing heading preserved"
assert_grep '^Some pre-existing user content\.$' CLAUDE.md "case3: existing prose preserved"
cd ..

pass
