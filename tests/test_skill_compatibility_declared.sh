#!/usr/bin/env bash
# tests/test_skill_compatibility_declared.sh — the REAL catalog, not a fixture.
#
# tests/test_validate_compatibility.sh proves the checker works against planted
# inputs. This one proves the shipped declarations are actually correct, and
# that adding them did not cost the catalog its clean bill of health:
#
#   1. `sciagent validate` exits 0 over the real toolkit (all 83 active skills).
#   2. Exactly the audited set of skills carries a `compatibility:` key —
#      neither more (a false positive re-declared) nor fewer (a declaration
#      lost to a merge).
#   3. Every declared value is within the 500-char canonical ceiling.
#   4. Every `sciagent-toolkit` item resolves to a real directory under lib/,
#      and every `sibling-skill` item to a real directory under skills/.
#      These are re-derived here from the files rather than trusted from the
#      validator, so a bug that made the validator vacuous (e.g. its
#      `compatibility:` grep silently matching nothing) cannot hide.
#
# The expected list below is the outcome of the 2026-08-11 semantic audit of
# the 15 skills flagged as "coupled". 7 of that flagged set were false
# positives and deliberately carry NO declaration; 6 skills the scan never
# flagged do. See docs/proposals/2026-08-11-offline-distribution/50_ADRs.md
# (ADR-D6) for the corrected numbers.

set -u
. "$(dirname "$0")/_lib.sh"

# --- the audited set ---------------------------------------------------------
EXPECTED=(
    annotate-bulk-rnaseq-data
    bulk-rnaseq-activity-inference
    bulk-rnaseq-gsea
    bulk-rnaseq-pathway-explorer
    cellranger-multi-to-anndata
    consensus-nmf-multirun
    decision-gate-notebook
    delegate-cli
    figure-style
    interactive-breakpoint-explorer
    iterative-peak-merging
    peak-atlas-multiome
    reasoning-trace
    scrna-cxg-host
    te-geneset-gsea
)

# Audited false positives: flagged by the static scan, cleared by the semantic
# audit. Re-declaring any of these would make the field lie in the other
# direction, so it is asserted explicitly rather than left to the count.
FALSE_POSITIVES=(
    anndatar-seurat-scanpy-conversion
    mllmcelltype-consensus-annotation
    peak-atlas-framework
    peak-atlas-unpaired
    coresh-signature-search
    scrna-pipeline-conventions
)

declare -i fail_count=0

# --- 0. the module wiring that makes the check reachable ---------------------
# `_validate_compatibility` calls `_fm_scalar` from lib/sciagent/frontmatter.sh.
# Every verb whose closure loads validate.sh must therefore also load
# frontmatter.sh — that is `validate` AND the three verbs that call
# `cmd_validate --quiet` as an internal pre-flight (activate, deactivate,
# update). Verified by hand: with frontmatter.sh dropped from those rows,
# `sciagent activate` against the real catalog dies with
# "_fm_scalar: command not found" on the first declaring skill, i.e. the
# mutation-guard is not academic. The fixture toolkits used by the activate
# tests declare no compatibility, so nothing else in the suite covers this.
while IFS= read -r row; do
    [[ "$row" == *validate.sh* ]] || continue
    if [[ "$row" != *frontmatter.sh* ]]; then
        echo "FAIL [$_TEST_NAME] bin/sciagent VERB_MODULES row loads validate.sh without frontmatter.sh:" >&2
        echo "    $row" >&2
        fail_count+=1
    fi
done < <(awk '/^declare -A VERB_MODULES=\(/,/^\)/' "$TOOLKIT_ROOT/bin/sciagent" | grep '^ *\[')

# --- 1. the whole catalog still validates ------------------------------------
set +e
out=$(SCIAGENT_TOOLKIT="$TOOLKIT_ROOT" "$TOOLKIT_ROOT/bin/sciagent" validate --quiet 2>&1)
rc=$?
set -e
if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] sciagent validate failed over the real toolkit (rc=$rc)" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# --- walk the catalog --------------------------------------------------------
declared=()
for skill_dir in "$TOOLKIT_ROOT"/skills/*/; do
    name="$(basename "$skill_dir")"
    [[ "$name" == _* ]] && continue
    file="$skill_dir/SKILL.md"
    [[ -f "$file" ]] || continue

    fm="$(awk 'NR==1 && /^---[ \t]*$/ {infm=1; next} infm && /^---[ \t]*$/ {exit} infm {print}' "$file")"
    printf '%s\n' "$fm" | grep -q '^compatibility:' || continue
    declared+=("$name")

    # Re-parse independently of lib/sciagent/frontmatter.sh: strip the outer
    # quotes only, so a mis-declared value cannot be normalised away here.
    value="$(printf '%s\n' "$fm" | awk '
        /^compatibility:/ {
            sub(/^compatibility:[ \t]*/, "")
            sub(/^"/, ""); sub(/"$/, "")
            print; exit
        }')"

    # --- 3. length ceiling ---
    if (( ${#value} > 500 )); then
        echo "FAIL [$_TEST_NAME] $name: compatibility is ${#value} chars (max 500)" >&2
        fail_count+=1
    fi

    # --- 4. toolkit-resolvable items really resolve ---
    # Split into "<flavour>: <items>" clauses on "; ", then items on ", ".
    while IFS= read -r pair; do
        [[ -z "$pair" ]] && continue
        flavour="${pair%%$'\t'*}"
        item="${pair#*$'\t'}"
        case "$flavour" in
            sciagent-toolkit)
                if [[ ! -d "$TOOLKIT_ROOT/lib/$item" ]]; then
                    echo "FAIL [$_TEST_NAME] $name: sciagent-toolkit '$item' has no lib/$item" >&2
                    fail_count+=1
                fi ;;
            sibling-skill)
                if [[ ! -d "$TOOLKIT_ROOT/skills/$item" ]]; then
                    echo "FAIL [$_TEST_NAME] $name: sibling-skill '$item' has no skills/$item" >&2
                    fail_count+=1
                fi ;;
            sciagent-scaffold)
                if [[ "$item" == /* || "$item" == *..* ]]; then
                    echo "FAIL [$_TEST_NAME] $name: sciagent-scaffold '$item' is not a clean relative path" >&2
                    fail_count+=1
                fi ;;
            external-module) ;;
            *)
                echo "FAIL [$_TEST_NAME] $name: unknown compatibility flavour '$flavour'" >&2
                fail_count+=1 ;;
        esac
    done < <(printf '%s\n' "$value" | awk '
        {
            nc = split($0, cl, "; ")
            for (i = 1; i <= nc; i++) {
                p = index(cl[i], ": ")
                if (p == 0) { print "<malformed>\t" cl[i]; continue }
                fl = substr(cl[i], 1, p - 1)
                ni = split(substr(cl[i], p + 2), it, ", ")
                for (j = 1; j <= ni; j++) print fl "\t" it[j]
            }
        }')
done

# --- 2. the declared set is exactly the audited set --------------------------
_sorted() { printf '%s\n' "$@" | sort; }
got="$(_sorted "${declared[@]}")"
want="$(_sorted "${EXPECTED[@]}")"
if [[ "$got" != "$want" ]]; then
    echo "FAIL [$_TEST_NAME] the set of skills declaring compatibility: has drifted" >&2
    echo "--- expected (audited) / +++ actual (on disk)" >&2
    diff <(printf '%s\n' "$want") <(printf '%s\n' "$got") >&2 || true
    fail_count+=1
fi

# The audited false positives must stay undeclared.
for fp in "${FALSE_POSITIVES[@]}"; do
    f="$TOOLKIT_ROOT/skills/$fp/SKILL.md"
    [[ -f "$f" ]] || continue
    if awk 'NR==1 && /^---[ \t]*$/ {infm=1; next} infm && /^---[ \t]*$/ {exit} infm' "$f" \
        | grep -q '^compatibility:'; then
        echo "FAIL [$_TEST_NAME] $fp is an audited false positive and must NOT declare compatibility:" >&2
        fail_count+=1
    fi
done

if [[ "$fail_count" -gt 0 ]]; then
    echo "FAIL [$_TEST_NAME] $fail_count compatibility problem(s) in the shipped catalog" >&2
    exit 1
fi

pass
