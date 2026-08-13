#!/usr/bin/env bash
# tests/test_template_provenance.sh
#
# templates/PROVENANCE.sha1 is what lets `link` tell "a stale copy the
# toolkit wrote" from "a file the user edited", and therefore what licenses it
# to overwrite the former. Two ways that guarantee can rot silently, both
# guarded here:
#
#   1. A managed template changes and nobody re-runs the generator. The new
#      hash is then absent from the manifest, so on the NEXT change every
#      consumer holding this version looks user-authored and gets ceded
#      instead of refreshed — the propagation defect comes back, quietly, for
#      one generation of the file.
#
#   2. A new hook template is added and never registered as managed. It
#      materializes into projects but cannot receive later refreshes.
#
#   3. A managed template gains a {{PLACEHOLDER}}. Content-provenance rests
#      entirely on the body being materialized by a plain `cp` — the moment a
#      template is rendered rather than copied, the project's bytes stop
#      matching the template's and every copy looks user-authored.
set -u
. "$(dirname "$0")/_lib.sh"

MANIFEST="$TOOLKIT_ROOT/templates/PROVENANCE.sha1"
assert_file_exists "$MANIFEST" "the provenance manifest is committed"

_hash() { sha1sum "$1" | cut -d' ' -f1; }

# The set of templates `link` materializes verbatim, derived from the hook and
# helper-shim source directories. Anything the caller materializes is checked
# here automatically, including the {{PLACEHOLDER}} invariant below.
MANAGED=()
for f in \
    "$TOOLKIT_ROOT"/templates/project/_common/.claude/hooks/*.sh.template \
    "$TOOLKIT_ROOT"/templates/project/analysis/02_analysis/helpers/*.template
do
    [[ -f "$f" ]] && MANAGED+=("${f#"$TOOLKIT_ROOT"/templates/}")
done
if [[ ${#MANAGED[@]} -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] no managed templates found — the glob is wrong" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# (1) + (2): every managed template's CURRENT hash is recorded.
# ---------------------------------------------------------------------------
for rel in "${MANAGED[@]}"; do
    cur=$(_hash "$TOOLKIT_ROOT/templates/$rel")
    if ! grep -q "^$cur  $rel\$" "$MANIFEST"; then
        echo "FAIL [$_TEST_NAME] current hash of $rel is NOT in the manifest" >&2
        echo "  hash: $cur" >&2
        echo "  fix : tools/gen-template-provenance.sh" >&2
        echo "  (until this is recorded, every consumer holding this version" >&2
        echo "   will be treated as user-authored and never refreshed)" >&2
        exit 1
    fi
done

# ---------------------------------------------------------------------------
# (3): managed templates must be copy-verbatim, never rendered.
# ---------------------------------------------------------------------------
for rel in "${MANAGED[@]}"; do
    if grep -qE '\{\{[A-Za-z0-9_]+\}\}' "$TOOLKIT_ROOT/templates/$rel"; then
        echo "FAIL [$_TEST_NAME] $rel contains a {{PLACEHOLDER}}" >&2
        echo "  Managed bodies are materialized with a plain cp, so a rendered" >&2
        echo "  template breaks content-provenance for every copy of it." >&2
        echo "  Either drop the placeholder or remove the file from the" >&2
        echo "  managed set in tools/gen-template-provenance.sh." >&2
        exit 1
    fi
done

# ---------------------------------------------------------------------------
# The generator agrees the manifest is current (also catches a hand-edit).
# Needs git history; skip loudly without it rather than fail a valid install.
# ---------------------------------------------------------------------------
if git -C "$TOOLKIT_ROOT" rev-parse --git-dir >/dev/null 2>&1; then
    if ! (cd "$TOOLKIT_ROOT" && tools/gen-template-provenance.sh --check >/dev/null 2>&1); then
        echo "FAIL [$_TEST_NAME] manifest is out of date — run tools/gen-template-provenance.sh" >&2
        exit 1
    fi
else
    echo "SKIP [$_TEST_NAME] no git history — generator --check not run" >&2
fi

# ---------------------------------------------------------------------------
# Manifest lines are well-formed: 40-hex, two spaces, a path under templates/
# that still exists. A typo'd hash would silently never match anything.
# ---------------------------------------------------------------------------
while IFS= read -r line; do
    [[ -z "$line" || "$line" == \#* ]] && continue
    if [[ ! "$line" =~ ^[0-9a-f]{40}\ \ [A-Za-z0-9_./-]+$ ]]; then
        echo "FAIL [$_TEST_NAME] malformed manifest line: '$line'" >&2
        exit 1
    fi
    p="${line#*  }"
    if [[ ! -f "$TOOLKIT_ROOT/templates/$p" ]]; then
        echo "FAIL [$_TEST_NAME] manifest names a template that no longer exists: $p" >&2
        exit 1
    fi
done < "$MANIFEST"

pass
