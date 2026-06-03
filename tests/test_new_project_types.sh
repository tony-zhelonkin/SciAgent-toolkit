#!/usr/bin/env bash
# Type-conditional docs/_internal/ namespace generation + _common overlay render.
# Closes the P3<->P6 seam: software-tool gets design/benchmarks (no sessions),
# analysis gets sessions/scratch/research (no design/benchmarks), and the
# shared _common/ files land in both.

. "$(dirname "$0")/_lib.sh"

setup_tmpdir

scaffold() {
    local dir="$1" type="$2"
    "$SCIAGENT_TOOLKIT/bin/sciagent" new project "$dir" --type "$type" >/dev/null 2>&1 \
        || { echo "FAIL [$_TEST_NAME] scaffold --type $type failed" >&2; exit 1; }
}

assert_no_dir() {
    local d="$1" msg="${2:-}"
    if [[ -d "$d" ]]; then
        echo "FAIL [$_TEST_NAME] ${msg:-dir should not exist: $d}" >&2
        exit 1
    fi
}

# (a) software-tool: design/benchmarks/reasoning present, sessions absent.
tool="$TMPDIR_TEST/tool"
scaffold "$tool" software-tool
assert_file_exists "$tool/docs/_internal/design"     "software-tool missing design/"
assert_file_exists "$tool/docs/_internal/benchmarks" "software-tool missing benchmarks/"
assert_file_exists "$tool/docs/_internal/reasoning"  "software-tool missing reasoning/"
assert_no_dir      "$tool/docs/_internal/sessions"   "software-tool should not have sessions/"
# Handoff routing must resolve to design/ (P3 Step-0 reads this row).
assert_grep "session handoff.*docs/_internal/design/" "$tool/AGENTS.md" \
    "software-tool AGENTS.md missing 'session handoff -> design/' routing row"

# (b) analysis: sessions/reasoning/research present,
#     design/benchmarks absent.
# (The P3 plan named a `scratch/` dir; the implementation flattened it out. That
#  divergence is internally consistent and out of scope here, so it is not asserted.)
ana="$TMPDIR_TEST/ana"
scaffold "$ana" analysis
assert_file_exists "$ana/docs/_internal/sessions"             "analysis missing sessions/"
assert_file_exists "$ana/docs/_internal/reasoning"            "analysis missing reasoning/"
assert_file_exists "$ana/docs/_internal/research"             "analysis missing research/"
assert_no_dir      "$ana/docs/_internal/design"               "analysis should not have design/"
assert_no_dir      "$ana/docs/_internal/benchmarks"           "analysis should not have benchmarks/"

# (c) _common/ files render into both types.
for proj in "$tool" "$ana"; do
    assert_file_exists "$proj/CLAUDE.md"                    "$proj missing _common CLAUDE.md"
    assert_file_exists "$proj/.gitignore"                  "$proj missing _common .gitignore"
    assert_file_exists "$proj/docs/_internal/README.md"    "$proj missing _common docs/_internal/README.md"
done

pass
