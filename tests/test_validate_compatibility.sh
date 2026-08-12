#!/usr/bin/env bash
# tests/test_validate_compatibility.sh — `sciagent validate` enforces the
# `compatibility:` declaration grammar.
#
# `compatibility:` is the ONE frontmatter key that states what a skill needs
# from outside itself. Its whole value is that it is machine-checkable, so a
# declaration that is silently wrong is worse than no declaration at all — a
# typo'd flavour (`scaffold:` for `sciagent-scaffold:`) reads correct to a
# human and enforces nothing. Every check below is therefore a HARD failure.
#
# Cases:
#    1. no compatibility: key at all              -> exit 0 (the 83-skill floor)
#    2. valid single-flavour declaration          -> exit 0
#    3. valid four-flavour declaration            -> exit 0
#    4. unknown flavour (`scaffold:`)             -> exit 1
#    5. sciagent-toolkit naming a nonexistent lib -> exit 1
#    6. sciagent-toolkit naming an existing lib   -> exit 0
#    7. sibling-skill naming a nonexistent skill  -> exit 1
#    8. sibling-skill naming an existing skill    -> exit 0
#    9. value over 500 chars                      -> exit 1
#   10. value at exactly 500 chars                -> exit 0 (inclusive cap)
#   11. absolute sciagent-scaffold path           -> exit 1
#   12. '..'-containing sciagent-scaffold path    -> exit 1
#   13. trailing whitespace after the value       -> exit 1
#   14. ";" not followed by exactly one space     -> exit 1
#   15. "," not followed by exactly one space     -> exit 1
#   16. flavour with no items / no ": " after it  -> exit 1
#   17. present but empty                         -> exit 1
#   18. unquoted (invalid YAML, ": " in value)    -> exit 1
#   19. external-module is NOT filesystem-checked -> exit 0
#   20. sciagent-scaffold is NOT filesystem-checked (project absent) -> exit 0
#   21. sciagent-toolkit item containing "/"      -> exit 1 (escape guard)
#
# Cases 6/8/10/19/20 are the negative controls: they are the inputs that would
# still fail if a check were written too broadly (e.g. if sciagent-scaffold
# were stat()ed, or if the 500 cap were exclusive). A test suite that only
# proves things fail cannot tell an over-strict check from a correct one.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir

# _plant <compatibility-line...> — build a fresh fake toolkit whose skills/s_a
# carries the given extra frontmatter lines (may be empty). Echoes the root.
#
# The fixture gains lib/figure-style and lib/interactive-style (the two real
# helper libs symlink_create_helper_lib mounts) so the sciagent-toolkit
# existence check has something true to resolve against; build_fake_toolkit
# already supplies skills s_a, s_b and s_c for sibling-skill.
_plant() {
    local extra="$1"
    local root="$TMPDIR_TEST/tk_$RANDOM$RANDOM"
    build_fake_toolkit "$root"
    mkdir -p "$root/lib/figure-style" "$root/lib/interactive-style"
    {
        printf -- '---\n'
        printf 'name: s_a\n'
        printf 'description: A perfectly ordinary fixture skill.\n'
        [[ -n "$extra" ]] && printf '%s\n' "$extra"
        printf -- '---\n\nbody\n'
    } > "$root/skills/s_a/SKILL.md"
    printf '%s\n' "$root"
}

# _run <toolkit-root> — run validate; sets $out (combined) and $rc.
_run() {
    set +e
    out=$(SCIAGENT_TOOLKIT="$1" "$1/bin/sciagent" validate 2>&1)
    rc=$?
    set -e
}

# _expect_pass <case-label> <compatibility-line>
_expect_pass() {
    local label="$1" line="$2" root
    root=$(_plant "$line")
    _run "$root"
    if [[ "$rc" -ne 0 ]]; then
        echo "FAIL [$_TEST_NAME] $label: expected validate to PASS, got rc=$rc" >&2
        echo "  input: $line" >&2
        printf '%s\n' "$out" >&2
        exit 1
    fi
}

# _expect_fail <case-label> <compatibility-line> <grep-pattern>
# The pattern assertion is what makes this a mutation test rather than a smoke
# test: a validate that fails for some UNRELATED reason would still exit 1.
_expect_fail() {
    local label="$1" line="$2" pattern="$3" root
    root=$(_plant "$line")
    _run "$root"
    if [[ "$rc" -eq 0 ]]; then
        echo "FAIL [$_TEST_NAME] $label: expected validate to FAIL, got rc=0" >&2
        echo "  input: $line" >&2
        printf '%s\n' "$out" >&2
        exit 1
    fi
    if ! printf '%s\n' "$out" | grep -q -- "$pattern"; then
        echo "FAIL [$_TEST_NAME] $label: failed, but not for the expected reason" >&2
        echo "  input:    $line" >&2
        echo "  expected: $pattern" >&2
        printf '%s\n' "$out" >&2
        exit 1
    fi
    # The finding must name the offending skill, or a 83-skill walk is useless.
    if ! printf '%s\n' "$out" | grep -q 's_a'; then
        echo "FAIL [$_TEST_NAME] $label: finding does not name the offending skill" >&2
        printf '%s\n' "$out" >&2
        exit 1
    fi
}

# --- 1. no key at all: every pre-existing skill must stay green --------------
_expect_pass "case1 (no compatibility: key)" ""

# --- 2/3. well-formed declarations -------------------------------------------
_expect_pass "case2 (single flavour)" \
    'compatibility: "sciagent-scaffold: 02_analysis/config/analysis_config.yaml"'
_expect_pass "case3 (four flavours, multi-item)" \
    'compatibility: "sciagent-toolkit: figure-style, interactive-style; sciagent-scaffold: 02_analysis/config/analysis_config.yaml, 03_results/; sibling-skill: s_b, s_c; external-module: RNAseq-toolkit"'

# --- 4. unknown flavour ------------------------------------------------------
# The highest-value check. `scaffold:` reads correct and declares nothing.
_expect_fail "case4 (unknown flavour)" \
    'compatibility: "scaffold: 02_analysis/"' \
    'unknown flavour "scaffold"'
_expect_fail "case4b (sciagent_toolkit with underscore)" \
    'compatibility: "sciagent_toolkit: figure-style"' \
    'unknown flavour "sciagent_toolkit"'

# --- 5/6. sciagent-toolkit resolves under the toolkit's lib/ -----------------
_expect_fail "case5 (nonexistent lib)" \
    'compatibility: "sciagent-toolkit: no-such-style"' \
    "no such directory lib/no-such-style"
_expect_pass "case6 (existing lib)" \
    'compatibility: "sciagent-toolkit: interactive-style"'

# --- 7/8. sibling-skill resolves under the toolkit's skills/ -----------------
_expect_fail "case7 (nonexistent sibling)" \
    'compatibility: "sibling-skill: no-such-skill"' \
    "sibling-skill 'no-such-skill' — no such skill"
_expect_pass "case8 (existing sibling)" \
    'compatibility: "sibling-skill: s_b"'

# --- 9/10. the 500-char ceiling (mirrors quick_validate.py) ------------------
# Build a value of an exact length: the "sciagent-scaffold: " prefix is 19
# chars, so pad the single item out to the target.
_compat_of_len() {
    local target="$1" prefix="sciagent-scaffold: " pad
    pad=$(head -c "$(( target - 19 ))" < /dev/zero | tr '\0' 'x')
    printf '%s%s' "$prefix" "$pad"
}
OVER=$(_compat_of_len 501)
assert_eq "${#OVER}" "501" "case9 fixture is not 501 chars"
_expect_fail "case9 (501 chars)" \
    "compatibility: \"$OVER\"" \
    'compatibility is 501 chars (max 500)'
AT=$(_compat_of_len 500)
assert_eq "${#AT}" "500" "case10 fixture is not 500 chars"
_expect_pass "case10 (exactly 500 chars, cap is inclusive)" \
    "compatibility: \"$AT\""

# --- 11/12. sciagent-scaffold path shape ------------------------------------
_expect_fail "case11 (absolute path)" \
    'compatibility: "sciagent-scaffold: /etc/passwd"' \
    "must be repo-root-relative"
_expect_fail "case12 (parent-escaping path)" \
    'compatibility: "sciagent-scaffold: ../../etc/passwd"' \
    "must not contain a '\.\.' component"
_expect_fail "case12b ('..' in the middle)" \
    'compatibility: "sciagent-scaffold: 03_results/../../etc"' \
    "must not contain a '\.\.' component"

# --- 13. trailing whitespace -------------------------------------------------
# _fm_scalar strips a trailing quote BEFORE trailing whitespace, so this must
# be caught on the RAW line: on the parsed value the defect shows up only as a
# stray quote character.
_expect_fail "case13 (trailing whitespace after the closing quote)" \
    'compatibility: "sciagent-scaffold: 03_results/" ' \
    'has trailing whitespace'
# The stray-quote residue on its own, with NO trailing whitespace to also trip
# the check above — an unbalanced closing quote. Without this the residue check
# would be untested (case13 masks it by failing on whitespace first).
_expect_fail "case13b (stray quote in the value, no trailing whitespace)" \
    'compatibility: "sciagent-scaffold: 03_results/x""' \
    'ends in a stray quote'

# --- 14/15/16. separator and item well-formedness ---------------------------
_expect_fail "case14 (';' with no following space)" \
    'compatibility: "sciagent-scaffold: 03_results/;sibling-skill: s_b"' \
    'clauses must be separated by "; "'
_expect_fail "case15 (',' with no following space)" \
    'compatibility: "sciagent-scaffold: 03_results/,02_analysis/"' \
    'must be separated by ", "'
# A trailing-space-only value collapses under _fm_scalar (which strips the
# closing quote, then the whitespace) to the bare flavour `sciagent-scaffold:`,
# so the grammar reports the missing ": " rather than the missing item. Either
# way it is a hard failure and the message points at the right clause.
_expect_fail "case16 (flavour with no items)" \
    'compatibility: "sciagent-scaffold: "' \
    'must be followed by ": "'
_expect_fail "case16b (empty first clause, items missing mid-value)" \
    'compatibility: "sciagent-scaffold: ; sibling-skill: s_b"' \
    'has trailing whitespace'
_expect_fail "case16c (no space after the flavour colon)" \
    'compatibility: "sciagent-scaffold:03_results/"' \
    'must be followed by ": "'

# --- 17. present but empty ---------------------------------------------------
_expect_fail "case17 (empty value)" \
    'compatibility: ""' \
    'present but empty'

# --- 18. unquoted scalar -----------------------------------------------------
# Every well-formed value contains ": ", which is not a legal plain YAML
# scalar — PyYAML and any harness frontmatter parser reject the whole file.
# _fm_scalar is line-based and would read it happily, so without this check
# `sciagent validate` is the only parser in the stack that accepts it.
_expect_fail "case18 (unquoted)" \
    'compatibility: sciagent-scaffold: 03_results/' \
    'must be a quoted scalar'

# --- 19/20. the deliberate NON-checks ----------------------------------------
# external-module is unverifiable by construction: the toolkit does not ship
# the thing being named, so there is nothing truthful to stat.
_expect_pass "case19 (external-module is never resolved)" \
    'compatibility: "external-module: some-toolkit-that-does-not-exist"'
# sciagent-scaffold is shape-only. `activate` calls `cmd_validate --quiet` as
# a pre-flight, so a validate that stat()ed consumer-project paths would let a
# project's state hard-block activation. The path below exists in no toolkit
# and must still pass.
_expect_pass "case20 (sciagent-scaffold is shape-only, never stat()ed)" \
    'compatibility: "sciagent-scaffold: 02_analysis/config/analysis_config.yaml, 99_nonexistent/deeply/nested/thing.yaml"'

# --- 21. escape guard on the resolved flavours -------------------------------
# The grammar says "a directory name under lib/", so a "/" is a grammar
# violation; rejecting it is also what stops `../../etc` from being stat()ed.
_expect_fail "case21 (sciagent-toolkit item containing '/')" \
    'compatibility: "sciagent-toolkit: ../../etc"' \
    'must be a bare directory name under the toolkit lib/ directory'
_expect_fail "case21b (sibling-skill item containing '/')" \
    'compatibility: "sibling-skill: ../../etc"' \
    'must be a bare directory name under skills/'

pass
