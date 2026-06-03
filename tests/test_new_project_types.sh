#!/usr/bin/env bash
# Type-conditional docs/_internal/ namespace generation + _common overlay render.
# Covers both first-class project types: analysis and software.
# analysis  → sessions/reasoning/research present; design/benchmarks absent.
# software  → design/benchmarks/reasoning present; sessions absent.
# Both receive the shared _common/ files (CLAUDE.md, .gitignore, docs/_internal/README.md).

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

# (a) analysis: sessions/reasoning/research present, design/benchmarks absent.
# (The P3 plan named a `scratch/` dir; the implementation flattened it out. That
#  divergence is internally consistent and out of scope here, so it is not asserted.)
ana="$TMPDIR_TEST/ana"
scaffold "$ana" analysis
assert_file_exists "$ana/docs/_internal/sessions"             "analysis missing sessions/"
assert_file_exists "$ana/docs/_internal/reasoning"            "analysis missing reasoning/"
assert_file_exists "$ana/docs/_internal/research"             "analysis missing research/"
assert_no_dir      "$ana/docs/_internal/design"               "analysis should not have design/"
assert_no_dir      "$ana/docs/_internal/benchmarks"           "analysis should not have benchmarks/"

# (b) software: design/benchmarks/reasoning present, sessions absent.
sw="$TMPDIR_TEST/sw"
scaffold "$sw" software
assert_file_exists "$sw/docs/_internal/design"                "software missing design/"
assert_file_exists "$sw/docs/_internal/benchmarks"            "software missing benchmarks/"
assert_file_exists "$sw/docs/_internal/reasoning"             "software missing reasoning/"
assert_no_dir      "$sw/docs/_internal/sessions"              "software should not have sessions/"
assert_no_dir      "$sw/docs/_internal/research"              "software should not have research/"

# (c) software: packageable-library source tree present.
assert_file_exists "$sw/src"                                  "software missing src/"
assert_file_exists "$sw/tests"                                "software missing tests/"
assert_file_exists "$sw/examples"                             "software missing examples/"

# (d) software: AGENTS.md rendered and says "software" not "software-tool".
assert_file_exists "$sw/AGENTS.md"                            "software missing AGENTS.md"
if grep -q "software-tool" "$sw/AGENTS.md"; then
    echo "FAIL [$_TEST_NAME] software AGENTS.md still contains 'software-tool'" >&2
    exit 1
fi
if ! grep -q "software" "$sw/AGENTS.md"; then
    echo "FAIL [$_TEST_NAME] software AGENTS.md missing 'software' marker" >&2
    exit 1
fi

# (e) _common/ files render into both types.
for proj in "$ana" "$sw"; do
    assert_file_exists "$proj/CLAUDE.md"                    "$proj missing _common CLAUDE.md"
    assert_file_exists "$proj/.gitignore"                  "$proj missing _common .gitignore"
    assert_file_exists "$proj/docs/_internal/README.md"    "$proj missing _common docs/_internal/README.md"
done

pass
