#!/usr/bin/env bash
# tests/test_helper_shim_coherence.sh — mount name <-> import name must agree.
#
# THE BUG THIS EXISTS TO PREVENT (shipped for two releases):
#   symlink_create_helper_lib mounts HYPHENATED directories under
#   02_analysis/helpers/ (figure-style, interactive-style). A hyphen is not
#   legal in a Python package name, so nothing can `import helpers.figure-style`.
#   The house design is therefore: hyphenated mount dir holding the contract lib
#   + an UNDERSCORED shim module next to it that projects import. figure-style
#   had its shim; interactive-style did not, and every import site in the repo
#   spelled `helpers.interactive_style.interactive_helpers` — a path that could
#   never resolve, no matter how fully the project was activated.
#
# Three coupled invariants, all derived from the source of truth rather than
# restated here (so a new mount / new import site is covered automatically):
#   (1) MOUNT -> SHIM: every directory in symlinks.sh's mount list has a
#       matching <underscored>.<ext>.template in the analysis project template
#       tree, one per language the lib actually provides.
#   (2) IMPORT -> SHIM: every `from helpers.X import` anywhere in the repo names
#       a single (dot-free) module X for which that shim template exists. A
#       dotted `helpers.a.b` spelling is a hard failure: the shim is a MODULE,
#       not a package, so it can never be traversed.
#   (3) SHIM -> MOUNT: every shim template delegates to a directory that is
#       actually in the mount list.
set -u
. "$(dirname "$0")/_lib.sh"

SYMLINKS="$TOOLKIT_ROOT/lib/sciagent/symlinks.sh"
HELPERS_TPL="$TOOLKIT_ROOT/templates/project/analysis/02_analysis/helpers"

# Source of truth for what gets mounted: the literal loop in symlinks.sh.
MOUNTS=$(sed -n 's/^[[:space:]]*for libdir in \(.*\); do$/\1/p' "$SYMLINKS")
if [[ -z "$MOUNTS" ]]; then
    echo "FAIL [$_TEST_NAME] could not parse the mount list out of $SYMLINKS" >&2
    echo "  (looked for: 'for libdir in <names>; do')" >&2
    exit 1
fi
echo "  mounts: $MOUNTS"

# ---------------------------------------------------------------------------
# (1) MOUNT -> SHIM, per language the lib provides.
# ---------------------------------------------------------------------------
declare -A SHIM_FOR_MOUNT=()
for m in $MOUNTS; do
    lib="$TOOLKIT_ROOT/lib/$m"
    assert_file_exists "$lib" "mount '$m' has no lib/ directory in this toolkit"
    under="${m//-/_}"
    SHIM_FOR_MOUNT["$under"]="$m"

    found_lang=0
    for ext in py R; do
        # Does the contract lib ship this language at all?
        compgen -G "$lib/*_helpers.$ext" >/dev/null || continue
        found_lang=1
        shim="$HELPERS_TPL/$under.$ext.template"
        if [[ ! -f "$shim" ]]; then
            echo "FAIL [$_TEST_NAME] mount '$m' ships a .$ext contract lib but has no importable shim" >&2
            echo "  expected: ${shim#"$TOOLKIT_ROOT"/}" >&2
            echo "  '$m' is hyphenated and cannot be imported directly; projects need '$under.$ext'." >&2
            exit 1
        fi
    done
    if (( found_lang == 0 )); then
        echo "FAIL [$_TEST_NAME] mount '$m' has no *_helpers.{py,R} contract lib" >&2
        exit 1
    fi
done

# ---------------------------------------------------------------------------
# (2) IMPORT -> SHIM. Scan every text source in the repo (skills, assets, docs,
#     templates, libs) for `from helpers.X import`.
# ---------------------------------------------------------------------------
sites=$(grep -rnE 'from[[:space:]]+helpers\.[A-Za-z0-9_.]+[[:space:]]+import' \
            --include='*.py' --include='*.qmd' --include='*.md' --include='*.R' \
            --include='*.template' --include='*.ipynb' \
            "$TOOLKIT_ROOT" 2>/dev/null | grep -v '/\.git/' | grep -v '/__pycache__/' || true)

if [[ -z "$sites" ]]; then
    echo "FAIL [$_TEST_NAME] found no 'from helpers.X import' sites at all — grep is broken" >&2
    exit 1
fi

n_sites=0
while IFS= read -r line; do
    [[ -n "$line" ]] || continue
    loc="${line%%:*}"
    mod=$(printf '%s\n' "$line" | sed -nE 's/.*from[[:space:]]+helpers\.([A-Za-z0-9_.]+)[[:space:]]+import.*/\1/p')
    [[ -n "$mod" ]] || continue
    n_sites=$((n_sites + 1))

    if [[ "$mod" == *.* ]]; then
        echo "FAIL [$_TEST_NAME] dotted helper import '$mod' at ${loc#"$TOOLKIT_ROOT"/}" >&2
        echo "  Shims are MODULES, not packages — 'helpers.${mod%%.*}' cannot be traversed further." >&2
        echo "  Import the flat shim: from helpers.${mod%%.*} import ..." >&2
        exit 1
    fi
    if [[ ! -f "$HELPERS_TPL/$mod.py.template" ]]; then
        echo "FAIL [$_TEST_NAME] '$mod' imported at ${loc#"$TOOLKIT_ROOT"/} but no shim ships it" >&2
        echo "  expected: templates/project/analysis/02_analysis/helpers/$mod.py.template" >&2
        echo "  (Scanning prose too is deliberate — a wrong import in a SKILL.md gets copy-pasted." >&2
        echo "   If this is a placeholder in a sentence rather than a real import, reword it.)" >&2
        exit 1
    fi
done <<< "$sites"
echo "  imports: $n_sites 'from helpers.X import' site(s), all resolvable"

# ---------------------------------------------------------------------------
# (3) SHIM -> MOUNT. A shim that points at a directory nobody mounts is dead.
# ---------------------------------------------------------------------------
n_shims=0
for shim in "$HELPERS_TPL"/*.template; do
    [[ -f "$shim" ]] || continue
    base="${shim##*/}"; base="${base%.template}"       # e.g. figure_style.py
    stem="${base%.*}"                                  # e.g. figure_style
    want="${SHIM_FOR_MOUNT[$stem]:-}"
    if [[ -z "$want" ]]; then
        continue    # a project-owned helper template, not a toolkit-lib shim
    fi
    n_shims=$((n_shims + 1))
    if ! grep -q -- "$want" "$shim"; then
        echo "FAIL [$_TEST_NAME] shim $base never references its mount dir '$want'" >&2
        exit 1
    fi
done
if (( n_shims == 0 )); then
    echo "FAIL [$_TEST_NAME] no shim template matched any mount — naming convention broke" >&2
    exit 1
fi
echo "  shims: $n_shims template(s) delegate to a real mount"

pass
