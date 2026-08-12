#!/usr/bin/env bash
# tests/test_lint_skill_coupling.sh — `sciagent lint --check skill-coupling`.
#
# The drift guard for skills' `compatibility:` declarations. `sciagent
# validate` checks a declaration's GRAMMAR; this check asks whether it is still
# TRUE. Half the test file is synthetic fixtures (does each rule fire?), half
# is run against the REAL shipped corpus (does it stay quiet where the
# 2026-08-11 semantic audit says it must?).
#
# Tests:
#   1.  undeclared coupling fires, with file:line
#   2.  a declared skill with the same code is silent
#   3.  a `#`-comment provenance citation is NOT evidence   <- the core design
#   4.  prose (.md) is NOT evidence
#   5.  stale declaration fires
#   6.  a mention anywhere (even a comment) acquits a declaration
#   7.  a scaffold root the toolkit does not scaffold fires
#   8.  a 01_modules/ scaffold item is told to use external-module
#   9.  undeclared sibling-skill fires only for a name that IS a skill
#   10. exit code is 0 with findings, and 0 even under --strict
#   11. --quiet suppresses
#   12. `--check all` does NOT include it
#   13. `activate` never runs it (and still activates)
#   14. REAL CORPUS: the 6 audited false positives produce no warning
#   15. REAL CORPUS: no declaring skill produces an undeclared-coupling warning
#   16. REAL CORPUS: exit code is 0
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

fail() {
    echo "FAIL [$_TEST_NAME] $1" >&2
    shift
    [[ $# -gt 0 ]] && { echo "--- output ---" >&2; printf '%s\n' "$@" >&2; }
    exit 1
}

# mkskill <name> <compatibility-value-or-empty>  — body on stdin is SKILL.md prose
mkskill() {
    local name="$1" compat="${2:-}"
    mkdir -p "$FAKE/skills/$name"
    {
        echo '---'
        echo "name: $name"
        echo "description: Fixture skill $name."
        [[ -n "$compat" ]] && echo "compatibility: \"$compat\""
        echo '---'
        cat
    } > "$FAKE/skills/$name/SKILL.md"
}

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

# (a) code coupling, NOT declared.
mkskill couple-undeclared <<'MD'
Body.
MD
mkdir -p "$FAKE/skills/couple-undeclared/scripts"
cat > "$FAKE/skills/couple-undeclared/scripts/run.py" <<'PY'
import anndata as ad
adata = ad.read_h5ad("03_results/objects/00_raw.h5ad")
PY

# (b) the same code, declared.
mkskill couple-declared "sciagent-scaffold: 03_results/objects/" <<'MD'
Body.
MD
mkdir -p "$FAKE/skills/couple-declared/scripts"
cat > "$FAKE/skills/couple-declared/scripts/run.py" <<'PY'
import anndata as ad
adata = ad.read_h5ad("03_results/objects/00_raw.h5ad")
PY

# (c) the ONLY scaffold reference is a provenance comment — the
#     peak-atlas-framework shape that produced the prior scan's false
#     positives. Both a bare and an absolute-path citation.
mkskill comment-only <<'MD'
Body.
MD
mkdir -p "$FAKE/skills/comment-only/scripts"
cat > "$FAKE/skills/comment-only/scripts/peaks.R" <<'R'
# Provenance:
#     /data2/users/Someone/proj/02_analysis/config/peak_filtration_config.yaml
#     03_results/objects/atlas.rds
library(GenomicRanges)
gr <- GRanges()
R

# (d) heavy prose coupling, no code at all — the scrna-pipeline-conventions
#     shape.
mkskill prose-only <<'MD'
Analysis scripts live under `02_analysis/stages/`, config in
`02_analysis/config/analysis_config.yaml`, and outputs land in
`03_results/<stage>/{figures,tables}/`.

```python
cfg = yaml.safe_load(open("02_analysis/config/analysis_config.yaml"))
```
MD
mkdir -p "$FAKE/skills/prose-only/references"
cat > "$FAKE/skills/prose-only/references/layout.md" <<'MD'
See `03_results/objects/01_qc.h5ad` and `00_data/raw/`.
MD

# (e) a declaration mentioned nowhere.
mkskill stale-decl "sciagent-scaffold: 03_results/unicorn/" <<'MD'
This skill says nothing about where its outputs go.
MD

# (f) a declaration mentioned only in a comment inside a script — rule 2
#     acquits broadly on purpose (any trace, anywhere, counts).
mkskill mention-in-comment "sciagent-scaffold: 03_results/griffin/" <<'MD'
Body.
MD
mkdir -p "$FAKE/skills/mention-in-comment/scripts"
cat > "$FAKE/skills/mention-in-comment/scripts/x.sh" <<'SH'
# writes into 03_results/griffin/ when the caller asks for it
echo hi
SH

# (g) scaffold root the toolkit does not scaffold.
mkskill bad-root "sciagent-scaffold: 07_elsewhere/tool.R" <<'MD'
Runs 07_elsewhere/tool.R.
MD

# (h) a sibling submodule filed as scaffold.
mkskill modules-root "sciagent-scaffold: 01_modules/RNAseq-toolkit" <<'MD'
Needs 01_modules/RNAseq-toolkit.
MD

# (i) sibling-skill: one reference to a REAL skill (s_a) and one to a
#     directory that is not a skill (the skill-creator usage-string shape).
mkskill sibling-user <<'MD'
Body.
MD
mkdir -p "$FAKE/skills/sibling-user/scripts"
cat > "$FAKE/skills/sibling-user/scripts/go.sh" <<'SH'
bash skills/s_a/scripts/helper.sh
echo "usage: package_skill.py skills/public/my-skill"
SH

run_check() { "$SCIAGENT" lint --check skill-coupling "$@" 2>&1; }

set +e
out=$(SCIAGENT_TOOLKIT="$FAKE" run_check); rc=$?
set -e

# ---------------------------------------------------------------------------
# 1-2: undeclared coupling fires with file:line; the declared twin is silent.
# ---------------------------------------------------------------------------
printf '%s\n' "$out" | grep -q "couple-undeclared: undeclared sciagent-scaffold coupling" \
    || fail "no undeclared-coupling warning for couple-undeclared" "$out"
printf '%s\n' "$out" | grep -qE "skills/couple-undeclared/scripts/run\.py:2" \
    || fail "undeclared-coupling warning does not cite file:line" "$out"
printf '%s\n' "$out" | grep -q "couple-declared:" \
    && fail "declared skill warned anyway" "$out"

# ---------------------------------------------------------------------------
# 3-4: comment lines and prose are not evidence.
# ---------------------------------------------------------------------------
printf '%s\n' "$out" | grep -q "comment-only:" \
    && fail "a #-comment provenance citation was counted as coupling" "$out"
printf '%s\n' "$out" | grep -q "prose-only:" \
    && fail "prose (.md) was counted as coupling" "$out"

# ---------------------------------------------------------------------------
# 5-6: stale declaration fires; a bare mention acquits.
# ---------------------------------------------------------------------------
printf '%s\n' "$out" | grep -q "stale-decl: stale declaration 'sciagent-scaffold: 03_results/unicorn/'" \
    || fail "no stale-declaration warning for stale-decl" "$out"
printf '%s\n' "$out" | grep -q "mention-in-comment:" \
    && fail "a declaration mentioned in a comment was called stale" "$out"

# ---------------------------------------------------------------------------
# 7-8: unscaffolded root, and 01_modules -> external-module.
# ---------------------------------------------------------------------------
printf '%s\n' "$out" | grep -q "bad-root: sciagent-scaffold '07_elsewhere/tool.R' is rooted at '07_elsewhere/'" \
    || fail "no unscaffolded-root warning for bad-root" "$out"
printf '%s\n' "$out" | grep -q "modules-root: sciagent-scaffold '01_modules/RNAseq-toolkit' names a sibling submodule" \
    || fail "no external-module hint for a 01_modules/ scaffold item" "$out"

# ---------------------------------------------------------------------------
# 9: sibling-skill fires for s_a (a real skill), not for skills/public/.
# ---------------------------------------------------------------------------
printf '%s\n' "$out" | grep -q "sibling-user: undeclared sibling-skill coupling.*skills/s_a/" \
    || fail "no undeclared sibling-skill warning for a reference to a real skill" "$out"
printf '%s\n' "$out" | grep -q "skills/public" \
    && fail "a reference to a non-skill directory was reported as sibling coupling" "$out"

# ---------------------------------------------------------------------------
# 10: exit code 0 with findings, and 0 under --strict too.
# ---------------------------------------------------------------------------
[[ "$rc" -eq 0 ]] || fail "check exited $rc with findings (must never be non-zero)" "$out"
set +e
sout=$(SCIAGENT_TOOLKIT="$FAKE" run_check --strict); src=$?
set -e
[[ "$src" -eq 0 ]] || fail "--strict promoted a skill-coupling finding to exit $src (must stay 0)" "$sout"
printf '%s\n' "$sout" | grep -q '^ERROR skill-coupling' \
    && fail "--strict turned a skill-coupling WARN into an ERROR" "$sout"

# ---------------------------------------------------------------------------
# 11: --quiet suppresses.
# ---------------------------------------------------------------------------
set +e
qout=$(SCIAGENT_TOOLKIT="$FAKE" run_check --quiet); qrc=$?
set -e
[[ "$qrc" -eq 0 ]] || fail "--quiet run exited $qrc" "$qout"
[[ -z "$qout" ]] || fail "--quiet did not suppress skill-coupling output" "$qout"

# ---------------------------------------------------------------------------
# 12: `--check all` (and the bare default) does NOT include skill-coupling.
# ---------------------------------------------------------------------------
PROJ="$TMPDIR_TEST/proj"; mkdir -p "$PROJ"
set +e
aout=$(SCIAGENT_TOOLKIT="$FAKE" "$SCIAGENT" lint --check all --project-dir "$PROJ" 2>&1)
bout=$(SCIAGENT_TOOLKIT="$FAKE" "$SCIAGENT" lint --project-dir "$PROJ" 2>&1)
set -e
printf '%s\n' "$aout" | grep -q 'skill-coupling' \
    && fail "--check all ran skill-coupling (its subject is the toolkit, not a project)" "$aout"
printf '%s\n' "$bout" | grep -q 'skill-coupling' \
    && fail "bare lint ran skill-coupling" "$bout"

# ...but the name is still accepted (not an unknown-check error).
set +e
nout=$(SCIAGENT_TOOLKIT="$FAKE" "$SCIAGENT" lint --check skill-coupling --quiet 2>&1); nrc=$?
set -e
[[ "$nrc" -eq 0 ]] || fail "--check skill-coupling rejected as unknown (exit $nrc)" "$nout"

# ---------------------------------------------------------------------------
# 13: `activate` must not run it. Its pre-flight is `validate --quiet`; a
#     drifting declaration must not appear there, nor block the mount.
# ---------------------------------------------------------------------------
APROJ="$TMPDIR_TEST/aproj"; mkdir -p "$APROJ"
(
    cd "$APROJ"
    set +e
    act=$(SCIAGENT_TOOLKIT="$FAKE" "$SCIAGENT" activate base 2>&1); arc=$?
    set -e
    [[ "$arc" -eq 0 ]] || { echo "FAIL [$_TEST_NAME] activate exited $arc with a drifting skill present" >&2
                            printf '%s\n' "$act" >&2; exit 1; }
    printf '%s\n' "$act" | grep -q 'skill-coupling' && {
        echo "FAIL [$_TEST_NAME] activate ran the skill-coupling check on its pre-flight path" >&2
        printf '%s\n' "$act" >&2; exit 1; }
    [[ -L .claude/skills/s_a ]] || { echo "FAIL [$_TEST_NAME] activate did not mount" >&2; exit 1; }
    exit 0
) || exit 1

# ---------------------------------------------------------------------------
# 14-16: THE REAL CORPUS. These are the calibration acceptance criteria.
# ---------------------------------------------------------------------------
set +e
real=$("$TOOLKIT_ROOT/bin/sciagent" lint --check skill-coupling 2>&1); realrc=$?
set -e
[[ "$realrc" -eq 0 ]] || fail "skill-coupling exited $realrc over the shipped corpus" "$real"

# The six false positives the 2026-08-11 semantic audit cleared. Each MENTIONS
# the analysis-repo layout and requires none of it; none carries a declaration.
# A warning on any of them means the comment/prose exclusion has regressed.
for fp in anndatar-seurat-scanpy-conversion mllmcelltype-consensus-annotation \
          peak-atlas-framework peak-atlas-unpaired coresh-signature-search \
          scrna-pipeline-conventions; do
    printf '%s\n' "$real" | grep -q "skill-coupling: $fp:" \
        && fail "audited false positive '$fp' produced a warning" "$real"
done

# No skill that DOES declare may be told it has undeclared coupling: that would
# mean the flavour matching is broken.
while IFS= read -r decl; do
    printf '%s\n' "$real" | grep -q "skill-coupling: $decl: undeclared" \
        && fail "declaring skill '$decl' reported as undeclared" "$real"
done < <(grep -l '^compatibility:' "$TOOLKIT_ROOT"/skills/*/SKILL.md \
         | xargs -r -n1 dirname | xargs -r -n1 basename)

pass
