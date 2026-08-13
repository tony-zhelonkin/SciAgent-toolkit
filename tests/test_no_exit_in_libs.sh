#!/usr/bin/env bash
# tests/test_no_exit_in_libs.sh — structural guard for the 5.5 rule:
# library functions (lib/sciagent/*.sh) only ever `return`, never `exit`.
# Only bin/sciagent may exit. `exit` inside an awk program / subshell passed as
# a single-quoted string is allowed (it tears down the awk/subshell, not the
# caller's shell). We therefore ignore any `exit` that sits inside a
# single-quoted region, and flag only shell-level `exit` statements.
#
# Full-line `#` comments are dropped BEFORE quote counting. The quoted-region
# tracker is a single-quote parity counter over the whole file, so an apostrophe
# in ordinary prose ("what's mounted") silently flips the parity and misclassifies
# every subsequent line. That produced a false positive on roles.sh once, and it
# can equally produce a false NEGATIVE that masks a real top-level exit — the
# failure mode this test exists to catch. A comment can never BE an exit, so
# skipping comments costs no coverage.
set -u
. "$(dirname "$0")/_lib.sh"

LIB="$TOOLKIT_ROOT/lib/sciagent"

# scan_shell_exits <file>
# Emits "<lineno>: <line>" for each top-level shell `exit`, tracking single-quote
# string state across lines so awk-embedded `exit` statements are skipped.
scan_shell_exits() {
    awk '
        { line = $0 }
        # Drop full-line comments before any quote accounting: prose apostrophes
        # would otherwise flip the parity tracker below.
        line ~ /^[[:space:]]*#/ { next }
        # Toggle in_squote on every unescaped single quote in the line.
        {
            n = gsub(/'\''/, "&", line)   # count single quotes on this line
            cur = in_sq
            # A line that is currently OUTSIDE a quote and whose first token is
            # "exit" is a shell-level exit (the offender).
            if (!cur && line ~ /^[[:space:]]*exit([[:space:]]|$)/) {
                printf "%d: %s\n", NR, $0
            }
            if (n % 2 == 1) in_sq = !in_sq
        }
    ' "$1"
}

violations=()
for f in "$LIB"/*.sh; do
    while IFS= read -r match; do
        [[ -z "$match" ]] && continue
        violations+=("$(basename "$f"): $match")
    done < <(scan_shell_exits "$f")
done

if (( ${#violations[@]} > 0 )); then
    echo "FAIL [$_TEST_NAME] top-level 'exit' found in library modules:" >&2
    printf '  %s\n' "${violations[@]}" >&2
    exit 1
fi

pass
