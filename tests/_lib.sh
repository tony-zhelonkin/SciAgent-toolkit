# tests/_lib.sh — shared assert helpers and tmpdir scaffold.
# Each test sources this, calls setup_tmpdir, runs assertions.
# On failure: print actual + expected + context, exit non-zero.

# shellcheck shell=bash

set -u

TOOLKIT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export SCIAGENT_TOOLKIT="$TOOLKIT_ROOT"

_TEST_NAME="${0##*/}"

assert_eq() {
    local actual="$1" expected="$2" msg="${3:-}"
    if [[ "$actual" != "$expected" ]]; then
        echo "FAIL [$_TEST_NAME] ${msg:-assert_eq}" >&2
        echo "  expected: $(printf '%q' "$expected")" >&2
        echo "  actual:   $(printf '%q' "$actual")" >&2
        exit 1
    fi
}

assert_file_eq() {
    local a="$1" b="$2" msg="${3:-}"
    if ! diff -u "$a" "$b" >/dev/null; then
        echo "FAIL [$_TEST_NAME] ${msg:-assert_file_eq}" >&2
        echo "--- expected ($b)" >&2
        echo "+++ actual   ($a)" >&2
        diff -u "$b" "$a" >&2 || true
        exit 1
    fi
}

assert_exit() {
    local expected="$1"; shift
    "$@"
    local rc=$?
    if [[ "$rc" -ne "$expected" ]]; then
        echo "FAIL [$_TEST_NAME] expected exit $expected got $rc: $*" >&2
        exit 1
    fi
}

assert_file_exists() {
    local f="$1" msg="${2:-}"
    if [[ ! -e "$f" ]]; then
        echo "FAIL [$_TEST_NAME] ${msg:-file missing: $f}" >&2
        exit 1
    fi
}

assert_symlink() {
    local f="$1" msg="${2:-}"
    if [[ ! -L "$f" ]]; then
        echo "FAIL [$_TEST_NAME] ${msg:-not a symlink: $f}" >&2
        ls -la "$(dirname "$f")" >&2 || true
        exit 1
    fi
}

assert_grep() {
    local pattern="$1" file="$2" msg="${3:-}"
    if ! grep -q "$pattern" "$file"; then
        echo "FAIL [$_TEST_NAME] ${msg:-pattern not found: $pattern in $file}" >&2
        echo "--- $file ---" >&2
        cat "$file" >&2
        exit 1
    fi
}

setup_tmpdir() {
    TMPDIR_TEST=$(mktemp -d)
    trap 'cleanup_tmpdir' EXIT
    cd "$TMPDIR_TEST"
}

cleanup_tmpdir() {
    if [[ -n "${TMPDIR_TEST:-}" && -d "$TMPDIR_TEST" ]]; then
        rm -rf "$TMPDIR_TEST"
    fi
}

pass() {
    echo "PASS [$_TEST_NAME]"
}

# build_fake_toolkit <dir>
# Create a minimal toolkit layout with fixture roles, skills, agents, commands.
# Roles defined:
#   base       — skills: s_a, s_b; agents: ag_a; commands: c_a
#   reviewer   — skills: s_b, s_c; agents: ag_b; commands: c_b (overlay)
#   alpha      — skills: s_a       (used in idempotent / max-stack tests)
build_fake_toolkit() {
    local root="$1"
    mkdir -p "$root"/{roles,skills/s_a,skills/s_b,skills/s_c,agents,commands,system-prompts,lib,bin}

    # Skills (directory format)
    echo "skill s_a"  > "$root/skills/s_a/SKILL.md"
    echo "skill s_b"  > "$root/skills/s_b/SKILL.md"
    echo "skill s_c"  > "$root/skills/s_c/SKILL.md"

    # Agents
    echo "agent ag_a" > "$root/agents/ag_a.md"
    echo "agent ag_b" > "$root/agents/ag_b.md"

    # Commands
    echo "cmd c_a"   > "$root/commands/c_a.md"
    echo "cmd c_b"   > "$root/commands/c_b.md"

    # System prompt fixture — exercised by tests that use output_style.
    # Frontmatter `name:` is the logical identifier resolved by
    # roles.sh system_prompt_path(); filename is incidental.
    cat > "$root/system-prompts/fixture-style.md" <<'EOF'
---
name: fixture-style
description: stub style for tests
---
fixture body
EOF

    # Roles
    cat > "$root/roles/base.yaml" <<EOF
name: base
description: fixture base role
skills:
  - s_a
  - s_b
agents:
  - ag_a
commands:
  - c_a
EOF
    cat > "$root/roles/reviewer.yaml" <<EOF
name: reviewer
description: fixture overlay
skills:
  - s_b
  - s_c
agents:
  - ag_b
commands:
  - c_b
EOF
    cat > "$root/roles/alpha.yaml" <<EOF
name: alpha
description: fixture alpha
skills:
  - s_a
EOF

    # Symlink lib/ and bin/ from the real toolkit so the dispatcher works.
    ln -sfn "$TOOLKIT_ROOT/lib/sciagent" "$root/lib/sciagent"
    ln -sfn "$TOOLKIT_ROOT/bin/sciagent" "$root/bin/sciagent"
}
