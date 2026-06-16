#!/usr/bin/env bash
# tests/strand_qc/run_regression.sh — strand-split QC regression for te-gene-featurecounts.
#
# Two assertions, both run against the locked te-fc:2.0.2 container:
#   (a) Synthetic truth-table re-run. For each frozen synthetic case, run
#       featureCounts 3x (-s0/-s2/-s1) with the recorded kernel, parse the
#       .summary + per-feature counts, and assert value-equal to the frozen
#       fixtures/synthetic/expected_truth_table.tsv. Assert every derived
#       per-fragment excess in {0,+1,+2} (the mechanism gate, end-to-end).
#   (b) Regime-classifier validation. Import the runnable QC suite's SHARED
#       classifier (qc/tools/05_classify_triples.py) and assert it reproduces the
#       locked ledger tallies (AP / M / P / bilateral / excess / residuals /
#       violations) from fixtures/rcore_fixture_0019.json. This makes the regression
#       guard the runnable suite's regime logic, not a duplicated copy.
#
# The multi-GB -R CORE BAMs that produced the fixture are NOT vendored; the
# fixture's joint_counts crosstab is the frozen witness they produced, and the
# synthetic SAF + per-case SAMs regenerate the truth table inside the container.
#
# Skips gracefully (exit 0) when docker or the te-fc:2.0.2 image is unavailable,
# so the toolkit-wide tests/run-all.sh stays green on machines without docker.
# Canonical doctrine and the warning-flag taxonomy live in the toolkit docs/QC.md;
# see references/strand-split-qc.md.
set -u

IMAGE="te-fc:2.0.2"
NAME="te-gene-featurecounts/strand_qc"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAF="$HERE/fixtures/synthetic/saf/features.saf"
SAM_DIR="$HERE/fixtures/synthetic/sam"
EXPECT="$HERE/fixtures/synthetic/expected_truth_table.tsv"
FIXTURE="$HERE/fixtures/rcore_fixture_0019.json"
# Import the RUNNABLE QC suite's shared classifier (qc/tools/05_classify_triples.py) so this
# regression now also guards the runnable suite (no duplicated regime logic). A thin
# classify_triples.py shim in this dir re-exports it for back-compat.
CLASSIFY="$HERE/../../qc/tools/05_classify_triples.py"
[[ -f "$CLASSIFY" ]] || CLASSIFY="$HERE/classify_triples.py"

skip() { echo "SKIP [$NAME] $*"; exit 0; }
fail() { echo "FAIL [$NAME] $*" >&2; exit 1; }
pass() { echo "PASS [$NAME] $*"; }

command -v docker >/dev/null 2>&1 || skip "docker not on PATH — skipping strand-split regression"
docker image inspect "$IMAGE" >/dev/null 2>&1 || \
  skip "image '$IMAGE' not built — run env/build.sh to enable the strand-split regression"
command -v python3 >/dev/null 2>&1 || skip "python3 not on PATH — skipping strand-split regression"

[[ -f "$SAF" ]]      || fail "missing synthetic SAF: $SAF"
[[ -f "$EXPECT" ]]   || fail "missing expected truth table: $EXPECT"
[[ -f "$FIXTURE" ]]  || fail "missing regime fixture: $FIXTURE"
[[ -f "$CLASSIFY" ]] || fail "missing classifier: $CLASSIFY"

# version check (matches the skill's smoke-test contract).
ver="$(docker run --rm "$IMAGE" featureCounts -v 2>&1 | grep -oE 'v2\.0\.2' | head -1 || true)"
[[ "$ver" == "v2.0.2" ]] || fail "featureCounts version is '$ver', expected v2.0.2"
pass "image $IMAGE present; featureCounts $ver"

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# ---------------------------------------------------------------------------
# (a) Synthetic truth-table re-run.
# ---------------------------------------------------------------------------
# Re-run featureCounts for every (case, strand) row of the frozen table, parse
# the .summary + per-feature counts in the container, and emit one observed line
# per (case, strand) in the SAME column layout as the expected table's value
# columns. The kernel per case is read from the table (se | se_noM | pe).

# kernel name -> flag string (matches SYNTHETIC_TRUTH.md "Exact flags used").
kernel_flags() {
  case "$1" in
    se)     echo "-M -F SAF" ;;
    se_noM) echo "-F SAF" ;;
    pe)     echo "-M -F SAF -p --countReadPairs -B -C" ;;
    *)      echo "__BADKERNEL__" ;;
  esac
}

# Unique (case, sam, kernel) tuples from the expected table (skip header).
mapfile -t CASES < <(awk -F'\t' 'NR>1{print $1"\t"$2"\t"$3}' "$EXPECT" | sort -u)

OBS="$WORK/observed.tsv"
: > "$OBS"

for tuple in "${CASES[@]}"; do
  case_id="$(printf '%s' "$tuple" | cut -f1)"
  sam_base="$(printf '%s' "$tuple" | cut -f2)"
  kernel="$(printf '%s' "$tuple" | cut -f3)"
  sam="$SAM_DIR/$sam_base"
  [[ -f "$sam" ]] || fail "case '$case_id' references missing SAM: $sam"
  flags="$(kernel_flags "$kernel")"
  [[ "$flags" != "__BADKERNEL__" ]] || fail "case '$case_id' has unknown kernel '$kernel'"

  for s in 0 1 2; do
    out="$WORK/${case_id}.s${s}.counts.txt"
    log="$WORK/${case_id}.s${s}.log"
    # shellcheck disable=SC2086
    docker run --rm -u "$(id -u):$(id -g)" -v "$WORK":"$WORK" -v "$HERE":"$HERE":ro "$IMAGE" \
      featureCounts $flags -a "$SAF" -o "$out" -s "$s" -T 2 "$sam" > "$log" 2>&1 \
      || { cat "$log" >&2; fail "featureCounts failed: $case_id s$s"; }

    sm="${out}.summary"
    [[ -f "$sm" && -f "$out" ]] || fail "no output for $case_id s$s"

    # channel label per pass (reverse-stranded library convention).
    case "$s" in 0) chan=s0 ;; 1) chan=anti ;; 2) chan=sense ;; esac

    # n_fragments = distinct read names in the SAM (PE mates share a name -> 1 frag).
    nfrag="$(awk '!/^@/{print $1}' "$sam" | sort -u | wc -l | tr -d ' ')"

    # summary tallies.
    read -r asg amb nof mm sgl < <(awk -F'\t' '
      $1=="Assigned"{a=$2}
      $1=="Unassigned_Ambiguity"{b=$2}
      $1=="Unassigned_NoFeatures"{c=$2}
      $1=="Unassigned_MultiMapping"{d=$2}
      $1=="Unassigned_Singleton"{e=$2}
      END{printf "%d %d %d %d %d\n", a+0,b+0,c+0,d+0,e+0}' "$sm")

    # per-feature nonzero counts, "feat=count;feat=count" sorted, or "-".
    feats="$(awk 'NR>2 && $7>0 {print $1"="$7}' "$out" | sort | paste -sd';' -)"
    [[ -n "$feats" ]] || feats="-"

    # emit the 12 leading value columns of the expected layout (excess col added below).
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$case_id" "$sam_base" "$kernel" "s$s" "$chan" "$nfrag" \
      "$asg" "$amb" "$nof" "$mm" "$sgl" "$feats" >> "$OBS"
  done
done

# Derive per-case per-fragment excess from the observed rows and append it,
# producing the full 13-column observed table to compare to the expected table.
OBS_FULL="$WORK/observed_full.tsv"
python3 - "$OBS" "$OBS_FULL" <<'PY'
import sys
obs_in, obs_out = sys.argv[1], sys.argv[2]
rows = [l.rstrip("\n").split("\t") for l in open(obs_in)]
# group by case; columns: 0 case,3 strand,5 nfrag,6 assigned
by = {}
for r in rows:
    by.setdefault(r[0], {})[r[3]] = r
out = []
for r in rows:
    case = r[0]
    g = by[case]
    asg = {s: int(g[s][6]) for s in ("s0", "s1", "s2")}
    nfrag = int(r[5])
    total = (asg["s2"] + asg["s1"]) - asg["s0"]
    if nfrag == 0 or total % nfrag != 0:
        sys.stderr.write(f"GATE-VIOLATION {case}: excess {total} not divisible by {nfrag} fragments\n")
        sys.exit(3)
    per = total // nfrag
    if per not in (0, 1, 2):
        sys.stderr.write(f"GATE-VIOLATION {case}: per-fragment excess {per} not in {{0,1,2}}\n")
        sys.exit(3)
    out.append(r + [str(per)])
with open(obs_out, "w") as fh:
    for r in out:
        fh.write("\t".join(r) + "\n")
PY
rc=$?
[[ $rc -eq 0 ]] || fail "per-fragment excess gate failed (excess not in {0,+1,+2})"
pass "mechanism gate: every derived per-fragment excess in {0,+1,+2}"

# Compare observed (full) against the expected table body (drop header line).
EXP_BODY="$WORK/expected_body.tsv"
tail -n +2 "$EXPECT" | LC_ALL=C sort > "$EXP_BODY"
OBS_SORTED="$WORK/observed_sorted.tsv"
LC_ALL=C sort "$OBS_FULL" > "$OBS_SORTED"

if ! diff -u "$EXP_BODY" "$OBS_SORTED" > "$WORK/diff.txt"; then
  echo "----- expected vs observed diff -----" >&2
  cat "$WORK/diff.txt" >&2
  fail "synthetic truth table mismatch vs fixtures/synthetic/expected_truth_table.tsv"
fi
ncases="${#CASES[@]}"
pass "synthetic truth table reproduced for $ncases cases x {s0,s2,s1} (value-equal to frozen expected table)"

# ---------------------------------------------------------------------------
# (b) Regime-classifier validation against the frozen fixture.
# ---------------------------------------------------------------------------
python3 - "$CLASSIFY" "$FIXTURE" <<'PY'
import importlib.util, json, sys
classify_path, fixture_path = sys.argv[1], sys.argv[2]
spec = importlib.util.spec_from_file_location("classify_triples", classify_path)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
fix = json.load(open(fixture_path))
res = mod.classify(fix["joint_counts"])

# the locked ledger (UNIT 3 contract).
LOCK = {
    "AP": 169032,
    "M_asymmetric": 40362,
    "P_silent": 448490,
    "bilateral": 2478,
    "excess_observed": 506656,
    "truthtable_violations": 192,
}
bad = []
for k, exp in LOCK.items():
    if res[k] != exp:
        bad.append(f"{k}: got {res[k]!r} expected {exp!r}")
# residuals (float).
def near(a, b, tol=1e-6):
    return abs(a - b) <= tol
if not near(res["resid_excess_pct"], -25.30908545443062, 1e-4):
    bad.append(f"resid_excess_pct: got {res['resid_excess_pct']!r} expected -25.31%")
if not near(res["resid_s0Amb_pct"], 0.0):
    bad.append(f"resid_s0Amb_pct: got {res['resid_s0Amb_pct']!r} expected 0.0")
if bad:
    sys.stderr.write("CLASSIFIER MISMATCH:\n  " + "\n  ".join(bad) + "\n")
    sys.exit(4)
print(f"classifier OK: AP={res['AP']} M={res['M_asymmetric']} P={res['P_silent']} "
      f"bilateral={res['bilateral']} excess={res['excess_observed']} "
      f"resid_excess={res['resid_excess_pct']:+.2f}% resid_s0Amb={res['resid_s0Amb_pct']:.1f} "
      f"violations={res['truthtable_violations']}")
PY
rc=$?
[[ $rc -eq 0 ]] || fail "classify_triples.py did not reproduce the locked ledger from the fixture"
pass "regime classifier reproduced the locked ledger from rcore_fixture_0019.json"

echo "PASS [$NAME] all strand-split regression assertions passed"
exit 0
