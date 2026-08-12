#!/usr/bin/env bash
# tests/test_helper_shim_propagation.sh
#
# The 02_analysis/helpers SHIM MODULES must reach a project that was already
# provisioned — not only one freshly scaffolded by `sciagent new project`.
#
# THE DEFECT. figure_style.py, figure_style.R and interactive_style.py were
# materialized exclusively by new.sh's generic template walk, which runs at
# scaffold time and never again. `sciagent activate` created the two HYPHENATED
# contract-lib mounts (02_analysis/helpers/{figure-style,interactive-style}) and
# no importable module beside them, so every skill whose `compatibility:`
# declaration names a shim path — figure-style, interactive-breakpoint-explorer,
# decision-gate-notebook — was satisfiable only in a brand-new repo. Third
# instance of one blind spot, after the hook bodies and the status line.
#
# The fix reuses the ownership discipline in lib/sciagent/ownership.sh (moved
# there out of claude_settings.sh, because an R helper is not a Claude artifact):
# create if absent, adopt a copy that already matches, refresh a copy whose
# bytes are any version this toolkit ever shipped, cede anything else.
#
# Cases (1)-(5) drive a real `sciagent activate`/`deactivate` against a
# SYNTHETIC toolkit, so every ownership branch — including "the manifest is what
# licenses the refresh" — is exercised deterministically without depending on
# this repo's git history. Cases (6) and (7) run against the REAL toolkit and the
# REAL manifest: (6) checks that what `new project` renders is byte-identical to
# what `activate` copies (the precondition content-provenance rests on), and (7)
# is the end-to-end fleet case, seeded with the oldest version of a shim that git
# says we ever shipped. Case (8) guards the module dependency closure.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir

_hash() { sha1sum "$1" | cut -d' ' -f1; }

# ---------------------------------------------------------------------------
# Synthetic toolkit: the fixture catalog + a helper-shim template tree and a
# provenance manifest we control completely. Two shims, one of which (the .py)
# has a "previously shipped" v1 and a "current" v2.
# ---------------------------------------------------------------------------
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"

# A real contract lib so symlink_create_helper_lib mounts something too — the
# shims live in the same directory as those mounts, and teardown has to cope.
mkdir -p "$FAKE/lib/figure-style"
echo "# stub R helper"  > "$FAKE/lib/figure-style/figure_helpers.R"
echo "# stub py helper" > "$FAKE/lib/figure-style/figure_helpers.py"

TPL_REL="project/analysis/02_analysis/helpers"
TPL_DIR="$FAKE/templates/$TPL_REL"
mkdir -p "$TPL_DIR"

PY_REL="$TPL_REL/figure_style.py.template"
R_REL="$TPL_REL/figure_style.R.template"

_write_py() { printf '"""figure_style shim v%s."""\nfrom helpers.figure_style import x  # v%s\n' "$1" "$1" > "$TPL_DIR/figure_style.py.template"; }
printf '## figure_style.R shim (figure-style)\nsource("figure-style/figure_helpers.R")\n' > "$TPL_DIR/figure_style.R.template"

_write_py 1
V1_BODY="$TMPDIR_TEST/py_v1.body"; cp "$TPL_DIR/figure_style.py.template" "$V1_BODY"
PY_V1=$(_hash "$V1_BODY")
_write_py 2
PY_V2=$(_hash "$TPL_DIR/figure_style.py.template")
R_CUR=$(_hash "$TPL_DIR/figure_style.R.template")

_write_manifest() {
    cat > "$FAKE/templates/PROVENANCE.sha1" <<EOF
# synthetic manifest
$PY_V2  $PY_REL
$PY_V1  $PY_REL
$R_CUR  $R_REL
EOF
}
_write_manifest

export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

STATE_DIR=".sciagent/helper_shim_state"

# ===========================================================================
# (1) An analysis project (02_analysis/ present) with NO shims: activate must
#     write every shim, byte-identical to the template, and record ownership.
# ===========================================================================
mkdir -p "$TMPDIR_TEST/p1/02_analysis"
cd "$TMPDIR_TEST/p1"

out=$("$SCIAGENT" activate base 2>&1)

assert_file_exists 02_analysis/helpers/figure_style.py "(1) .py shim materialized by activate"
assert_file_exists 02_analysis/helpers/figure_style.R  "(1) .R shim materialized by activate"
assert_eq "$(_hash 02_analysis/helpers/figure_style.py)" "$PY_V2" "(1) .py shim is byte-identical to the template"
assert_eq "$(_hash 02_analysis/helpers/figure_style.R)"  "$R_CUR"  "(1) .R shim is byte-identical to the template"
assert_eq "$(cat "$STATE_DIR/figure_style.py.sha1")" "$PY_V2" "(1) ownership record written for the .py shim"
assert_eq "$(cat "$STATE_DIR/figure_style.R.sha1")"  "$R_CUR"  "(1) ownership record written for the .R shim"
case "$out" in *"wrote: 02_analysis/helpers/figure_style.py"*) : ;;
    *) echo "FAIL [$_TEST_NAME] (1) activate did not report writing the shim: $out" >&2; exit 1 ;; esac
# Imported, never executed: the execute bit must NOT be asserted on these
# (mode `plain`), unlike the hooks and the status line.
[[ -x 02_analysis/helpers/figure_style.py ]] && {
    echo "FAIL [$_TEST_NAME] (1) shim was made executable; it is imported, not run" >&2; exit 1; }
# Re-running must be silent and change nothing (idempotent adopt path).
out_again=$("$SCIAGENT" activate base 2>&1)
case "$out_again" in *"02_analysis/helpers/figure_style"*)
    echo "FAIL [$_TEST_NAME] (1) second activate re-reported the shims: $out_again" >&2; exit 1 ;; esac
assert_eq "$(_hash 02_analysis/helpers/figure_style.py)" "$PY_V2" "(1) unchanged on re-activate"

# ===========================================================================
# (2) A project with NO 02_analysis/: activate must write no shim and must not
#     create the directory. A coordination/software repo stays untouched.
# ===========================================================================
mkdir -p "$TMPDIR_TEST/p2"
cd "$TMPDIR_TEST/p2"
"$SCIAGENT" activate base >/dev/null 2>&1
[[ -e 02_analysis ]] && { echo "FAIL [$_TEST_NAME] (2) activate created 02_analysis/ in a non-analysis repo" >&2; exit 1; }
[[ -e "$STATE_DIR" ]] && { echo "FAIL [$_TEST_NAME] (2) ownership state written with no analysis layout" >&2; exit 1; }

# ===========================================================================
# (3) THE FLEET'S STATE. A stale shim byte-identical to a version we shipped
#     earlier, with NO ownership record, must be refreshed to the current
#     template.
# ===========================================================================
mkdir -p "$TMPDIR_TEST/p3/02_analysis/helpers"
cd "$TMPDIR_TEST/p3"
cp "$V1_BODY" 02_analysis/helpers/figure_style.py
[[ -e "$STATE_DIR" ]] && { echo "FAIL [$_TEST_NAME] (3) fixture seeded an ownership record" >&2; exit 1; }

out3=$("$SCIAGENT" activate base 2>&1)
assert_eq "$(_hash 02_analysis/helpers/figure_style.py)" "$PY_V2" \
    "(3) stale-but-ours shim REFRESHED to the current template"
assert_eq "$(cat "$STATE_DIR/figure_style.py.sha1")" "$PY_V2" "(3) record updated to the refreshed content"
case "$out3" in *"refreshed: 02_analysis/helpers/figure_style.py"*) : ;;
    *) echo "FAIL [$_TEST_NAME] (3) no 'refreshed:' line: $out3" >&2; exit 1 ;; esac

# Non-vacuity: the refresh must be licensed by the manifest, not fire blindly.
# With the manifest hidden, the SAME stale shim must be left exactly alone.
mkdir -p "$TMPDIR_TEST/p3b/02_analysis/helpers"
cd "$TMPDIR_TEST/p3b"
cp "$V1_BODY" 02_analysis/helpers/figure_style.py
mv "$FAKE/templates/PROVENANCE.sha1" "$FAKE/templates/PROVENANCE.sha1.hidden"
"$SCIAGENT" activate base >/dev/null 2>&1
mv "$FAKE/templates/PROVENANCE.sha1.hidden" "$FAKE/templates/PROVENANCE.sha1"
assert_eq "$(_hash 02_analysis/helpers/figure_style.py)" "$PY_V1" \
    "(3b) without the manifest the stale shim is NOT overwritten (the fix is non-vacuous)"

# ===========================================================================
# (4) A USER-EDITED shim — content matching no version we ever shipped — is
#     kept, warned about ONCE, ceded permanently, and silent from then on.
# ===========================================================================
mkdir -p "$TMPDIR_TEST/p4/02_analysis/helpers"
cd "$TMPDIR_TEST/p4"
printf '"""my own shim."""\nMINE = True\n' > 02_analysis/helpers/figure_style.py
USER_HASH=$(_hash 02_analysis/helpers/figure_style.py)

err4=$("$SCIAGENT" activate base 2>&1 >/dev/null)
assert_eq "$(_hash 02_analysis/helpers/figure_style.py)" "$USER_HASH" "(4) user-edited shim NOT overwritten"
case "$err4" in *"will not manage it"*) : ;;
    *) echo "FAIL [$_TEST_NAME] (4) no cede warning: $err4" >&2; exit 1 ;; esac
assert_file_exists "$STATE_DIR/figure_style.py.ceded" "(4) ceded marker written"
assert_eq "$(cat "$STATE_DIR/figure_style.py.ceded")" "$USER_HASH" "(4) marker records the ceded content"
[[ -f "$STATE_DIR/figure_style.py.sha1" ]] && {
    echo "FAIL [$_TEST_NAME] (4) an ownership record was kept for a ceded file" >&2; exit 1; }

err5=$("$SCIAGENT" activate base 2>&1 >/dev/null)
case "$err5" in *"will not manage it"*)
    echo "FAIL [$_TEST_NAME] (4b) re-warned about an already-ceded shim: $err5" >&2; exit 1 ;; esac
assert_eq "$(_hash 02_analysis/helpers/figure_style.py)" "$USER_HASH" "(4b) still untouched on the second run"

# ...and deactivate leaves the user's file alone while dropping our marker.
"$SCIAGENT" deactivate >/dev/null 2>&1
assert_file_exists 02_analysis/helpers/figure_style.py "(4c) ceded shim survives deactivate"
assert_eq "$(_hash 02_analysis/helpers/figure_style.py)" "$USER_HASH" "(4c) ceded shim byte-identical after deactivate"
[[ -e "$STATE_DIR" ]] && { echo "FAIL [$_TEST_NAME] (4c) ceded marker/state dir survived deactivate" >&2; exit 1; }

# ===========================================================================
# (5) deactivate REVERSES the whole thing: every shim we wrote is removed, the
#     ownership state is gone, 02_analysis/helpers/ is pruned once empty, and
#     02_analysis/ itself — the project's own directory — is untouched.
# ===========================================================================
mkdir -p "$TMPDIR_TEST/p5/02_analysis"
cd "$TMPDIR_TEST/p5"
"$SCIAGENT" activate base >/dev/null 2>&1
assert_file_exists 02_analysis/helpers/figure_style.py "(5) shim present before deactivate"

"$SCIAGENT" deactivate >/dev/null 2>&1
[[ -e 02_analysis/helpers/figure_style.py ]] && { echo "FAIL [$_TEST_NAME] (5) shim survived deactivate" >&2; exit 1; }
[[ -e 02_analysis/helpers/figure_style.R ]]  && { echo "FAIL [$_TEST_NAME] (5) .R shim survived deactivate" >&2; exit 1; }
[[ -e "$STATE_DIR" ]] && { echo "FAIL [$_TEST_NAME] (5) ownership state survived deactivate" >&2; exit 1; }
[[ -e 02_analysis/helpers ]] && { echo "FAIL [$_TEST_NAME] (5) empty helpers/ dir left behind: $(find 02_analysis/helpers)" >&2; exit 1; }
assert_file_exists 02_analysis "(5) the project's own 02_analysis/ must never be removed"

# (5b) A user-edited shim we DID write is kept with a warning, and a sibling
#      file in helpers/ keeps the directory alive — rmdir, never rm -r.
mkdir -p "$TMPDIR_TEST/p5b/02_analysis"
cd "$TMPDIR_TEST/p5b"
"$SCIAGENT" activate base >/dev/null 2>&1
echo "# user addition" >> 02_analysis/helpers/figure_style.py
out5b=$("$SCIAGENT" deactivate 2>&1)
printf '%s\n' "$out5b" | grep -q "figure_style.py was modified" || {
    echo "FAIL [$_TEST_NAME] (5b) expected a drift warning for the edited shim" >&2
    printf '%s\n' "$out5b" >&2; exit 1; }
assert_grep '# user addition' 02_analysis/helpers/figure_style.py "(5b) the user's edit survives"
[[ -e 02_analysis/helpers/figure_style.R ]] && {
    echo "FAIL [$_TEST_NAME] (5b) the untouched shim should still have been removed" >&2; exit 1; }
assert_file_exists 02_analysis/helpers "(5b) helpers/ stays while the user's file is in it"

# ===========================================================================
# (6) THE PRECONDITION content-provenance rests on: the bytes `new project`
#     renders and the bytes `activate` copies must be THE SAME. new.sh renders
#     through sed (_subst); activate uses a plain cp. They agree only while the
#     templates carry no {{PLACEHOLDER}} — checked statically by
#     tests/test_template_provenance.sh — and while _subst keeps emitting its
#     input unchanged when there is nothing to substitute. This asserts the
#     end-to-end consequence: a freshly scaffolded repo's shims are recognised
#     as ours (adopted silently, never ceded) on its first activate.
# ===========================================================================
cd "$TMPDIR_TEST"
export SCIAGENT_TOOLKIT="$TOOLKIT_ROOT"
"$TOOLKIT_ROOT/bin/sciagent" new project "$TMPDIR_TEST/np" --type analysis >/dev/null 2>&1 \
    || { echo "FAIL [$_TEST_NAME] (6) new project --type analysis failed" >&2; exit 1; }
n_cmp=0
for tpl in "$TOOLKIT_ROOT"/templates/project/analysis/02_analysis/helpers/*.template; do
    [[ -f "$tpl" ]] || continue
    b="$(basename "${tpl%.template}")"
    assert_file_exists "$TMPDIR_TEST/np/02_analysis/helpers/$b" "(6) new project rendered $b"
    cmp -s "$TMPDIR_TEST/np/02_analysis/helpers/$b" "$tpl" || {
        echo "FAIL [$_TEST_NAME] (6) rendered $b is NOT byte-identical to its template" >&2
        echo "  content-provenance cannot recognise it, so every scaffolded repo's copy" >&2
        echo "  would be ceded to the user on its first activate." >&2
        diff <(cat "$TMPDIR_TEST/np/02_analysis/helpers/$b") <(cat "$tpl") >&2 | head -20
        exit 1; }
    n_cmp=$((n_cmp + 1))
done
[[ "$n_cmp" -ge 1 ]] || { echo "FAIL [$_TEST_NAME] (6) compared no shims — the glob is wrong" >&2; exit 1; }

# And the first activate of that scaffolded repo must adopt them silently:
# no write, no refresh, no cede warning, and a record so deactivate can reverse.
cd "$TMPDIR_TEST/np"
out6=$("$TOOLKIT_ROOT/bin/sciagent" activate base 2>&1)
case "$out6" in
    *"02_analysis/helpers/figure_style"*)
        echo "FAIL [$_TEST_NAME] (6) activate rewrote/ceded a freshly scaffolded shim: $out6" >&2
        exit 1 ;;
esac
assert_file_exists "$STATE_DIR/figure_style.py.sha1" "(6) scaffolded shim adopted (record written)"
[[ -e "$STATE_DIR/figure_style.py.ceded" ]] && {
    echo "FAIL [$_TEST_NAME] (6) a scaffolded shim was ceded instead of adopted" >&2; exit 1; }
cd "$TMPDIR_TEST"

# ===========================================================================
# (7) END-TO-END against the REAL toolkit and the REAL manifest: reproduce the
#     fleet state exactly. Seed an analysis project with the OLDEST shipped
#     version of a real shim and no ownership record, then run a real activate.
# ===========================================================================
cd "$TMPDIR_TEST"
export SCIAGENT_TOOLKIT="$TOOLKIT_ROOT"
REAL_REL="project/analysis/02_analysis/helpers/figure_style.py.template"
REAL_TPL="$TOOLKIT_ROOT/templates/$REAL_REL"
OLD_SHA=$(git -C "$TOOLKIT_ROOT" log --all --format=%H -- "templates/$REAL_REL" 2>/dev/null | tail -1)

if [[ -z "$OLD_SHA" ]] || ! git -C "$TOOLKIT_ROOT" show "$OLD_SHA:templates/$REAL_REL" >/dev/null 2>&1; then
    # Loudly, not silently: a shallow clone would otherwise make this vacuous.
    echo "SKIP [$_TEST_NAME] (7) no git history for $REAL_REL — end-to-end case not run" >&2
else
    git -C "$TOOLKIT_ROOT" show "$OLD_SHA:templates/$REAL_REL" > "$TMPDIR_TEST/old_shim.body"
    OLD_HASH=$(_hash "$TMPDIR_TEST/old_shim.body")
    CUR_HASH=$(_hash "$REAL_TPL")

    if [[ "$OLD_HASH" == "$CUR_HASH" ]]; then
        echo "SKIP [$_TEST_NAME] (7) template has only ever had one version" >&2
    else
        mkdir -p e2e/02_analysis/helpers && cd e2e
        printf '# Project AGENTS.md\n\nUser-owned content.\n' > AGENTS.md
        cp "$TMPDIR_TEST/old_shim.body" 02_analysis/helpers/figure_style.py
        # No .sciagent/helper_shim_state — exactly the fleet's state.
        [[ -e "$STATE_DIR" ]] && { echo "FAIL [$_TEST_NAME] (7) fixture seeded a record" >&2; exit 1; }

        "$TOOLKIT_ROOT/bin/sciagent" activate base >/dev/null 2>&1

        assert_eq "$(_hash 02_analysis/helpers/figure_style.py)" "$CUR_HASH" \
            "(7) real activate refreshed the stale shim to the current template"
        assert_file_exists 02_analysis/helpers/interactive_style.py \
            "(7) real activate also materializes the shim that never existed in the field"
        assert_file_exists 02_analysis/helpers/figure_style.R \
            "(7) real activate materializes the R shim"
        # The whole point of the shim: it must be importable next to its
        # hyphenated mount, which activate creates in the same directory.
        assert_symlink 02_analysis/helpers/figure-style "(7) contract-lib mount present beside the shim"

        "$TOOLKIT_ROOT/bin/sciagent" deactivate >/dev/null 2>&1
        [[ -e 02_analysis/helpers/figure_style.py ]] && {
            echo "FAIL [$_TEST_NAME] (7) real deactivate left the refreshed shim behind" >&2; exit 1; }
        cd ..
    fi
fi

# ===========================================================================
# (8) MODULE DEPENDENCY CLOSURE. ownership.sh is a third module that both
#     claude_settings.sh and symlinks.sh call into. A verb whose closure loads
#     a caller but not ownership.sh dies with "command not found" at the moment
#     it mutates — the failure mode of 4b291b2, which is invisible until a
#     specific verb runs a specific branch. Assert it structurally instead.
# ===========================================================================
BIN="$TOOLKIT_ROOT/bin/sciagent"
LIBDIR="$TOOLKIT_ROOT/lib/sciagent"

# Which lib modules actually call an ownership_* function (comments stripped)?
callers=()
for f in "$LIBDIR"/*.sh; do
    b="${f##*/}"
    [[ "$b" == "ownership.sh" ]] && continue
    if sed 's/#.*//' "$f" | grep -qE '\bownership_[a-z_]+' ; then
        callers+=("$b")
    fi
done
if [[ ${#callers[@]} -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] (8) no module calls ownership_* — the grep is wrong" >&2
    exit 1
fi

# Every VERB_MODULES entry that loads a caller must also load ownership.sh.
while IFS= read -r line; do
    verb="${line%%]*}"; verb="${verb#*[}"
    mods="${line#*\"}"; mods="${mods%\"*}"
    for c in "${callers[@]}"; do
        case " $mods " in
            *" $c "*)
                case " $mods " in
                    *" ownership.sh "*) : ;;
                    *)
                        echo "FAIL [$_TEST_NAME] (8) verb '$verb' loads $c but not ownership.sh" >&2
                        echo "  modules: $mods" >&2
                        echo "  $c calls ownership_* — this is a 'command not found' waiting to happen." >&2
                        exit 1 ;;
                esac ;;
        esac
    done
done < <(sed -n '/^declare -A VERB_MODULES=(/,/^)/p' "$BIN" | grep -E '^\s*\[[a-z]+\]="')

pass
