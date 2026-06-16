#!/usr/bin/env bash
# run_qc.sh — master orchestrator for the strand-split TE QC suite.
#
# QC a NEW TE dataset end-to-end: run the -R CORE regime witness on a sample, prove the ledger
# closes to the read, attribute the silent pile with the -O target-revealer, gate the young set,
# and corroborate against read-independent SAF geometry. Container/bedtools tools SKIP GRACEFULLY
# (exit 0 with a clear message) when their dependency is absent, so the suite always reaches a
# consolidated verdict. Operator-invoked; not a regression (that is tests/strand_qc/).
#
# Order (04 feeds 02/05/06/07/08):
#   04 -R CORE 3-pass witness            -> joint_counts + s0/s2/s1 CORE
#   05 regime classifier (SHARED module) -> AP/M/P ledger from joint_counts
#   06 closure-table audit               -> two identities + 14-cell sum + denominator discipline
#   02 fragment-weighting check          -> matrix==nfrag (needs the s2 CORE + matrix + BAM)
#   01 Random-One BAM witness            -> 0 fragments emit >1 locus (needs the BAM)
#   07 -O silent attribution (.sh + .py) -> silent-pile young share
#   08 young assignable-evidence gate    -> conserved-fraction gate
#   09 SAF self-intersect geometry       -> geometry-vs-witness concordance
#   03 strand-invariant + per-sample meter (needs .summary files; for de_precheck handoff)
#   de_precheck                          -> 3-metric condition-linkage (needs --design)
#
# Usage:
#   run_qc.sh BAM_DIR SAF STRAND OUTDIR [--sample ID] [--design design.tsv]
#       [--matrix te_sense_s2.counts.txt] [--s0-summary F --s2-summary F --s1-summary F --gene-summary F]
#       [--s0-matrix F --s2-matrix F --s1-matrix F] [--young-regex RE]
#       [--image te-fc:2.0.2] [--threads 12] [--bam-glob '*.markdup.sorted.bam']
#
#   BAM_DIR : dir of star_salmon BAMs (one --sample is selected for the per-read witness)
#   SAF     : the grouped TE SAF (GeneID Chr Start End Strand)
#   STRAND  : library strand of the dataset (0/1/2) — recorded; the witness always runs all 3 passes
#   OUTDIR  : output directory
set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLS="$HERE/tools"
YOUNG_RE='^(L1MdT|L1MdGf|L1MdA|IAPEz)'
IMAGE="te-fc:2.0.2"
THREADS=12
BAM_GLOB='*.markdup.sorted.bam'
SAMPLE=""; DESIGN=""; MATRIX=""
S0_SUM=""; S2_SUM=""; S1_SUM=""; GENE_SUM=""
S0_MAT=""; S2_MAT=""; S1_MAT=""

if [[ $# -lt 4 ]]; then
  grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 2
fi
BAM_DIR="$1"; SAF="$2"; STRAND="$3"; OUTDIR="$4"; shift 4

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample)       SAMPLE="$2"; shift 2 ;;
    --design)       DESIGN="$2"; shift 2 ;;
    --matrix)       MATRIX="$2"; shift 2 ;;
    --s0-summary)   S0_SUM="$2"; shift 2 ;;
    --s2-summary)   S2_SUM="$2"; shift 2 ;;
    --s1-summary)   S1_SUM="$2"; shift 2 ;;
    --gene-summary) GENE_SUM="$2"; shift 2 ;;
    --s0-matrix)    S0_MAT="$2"; shift 2 ;;
    --s2-matrix)    S2_MAT="$2"; shift 2 ;;
    --s1-matrix)    S1_MAT="$2"; shift 2 ;;
    --young-regex)  YOUNG_RE="$2"; shift 2 ;;
    --image)        IMAGE="$2"; shift 2 ;;
    --threads)      THREADS="$2"; shift 2 ;;
    --bam-glob)     BAM_GLOB="$2"; shift 2 ;;
    *) echo "[run_qc] unknown arg: $1" >&2; exit 2 ;;
  esac
done

mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"

# pick the sample BAM for the per-read witness
BAM=""
if [[ -n "$SAMPLE" ]]; then
  for cand in "$BAM_DIR/$SAMPLE.markdup.sorted.bam" "$BAM_DIR/$SAMPLE.bam" "$BAM_DIR/$SAMPLE"; do
    [[ -f "$cand" ]] && { BAM="$cand"; break; }
  done
else
  # shellcheck disable=SC2086
  BAM="$(ls $BAM_DIR/$BAM_GLOB 2>/dev/null | head -1)"
fi

echo "############################################################################"
echo "# strand-split TE QC suite — $(basename "$OUTDIR")"
echo "#   BAM_DIR=$BAM_DIR  SAF=$SAF  STRAND=$STRAND"
echo "#   witness sample BAM=${BAM:-<none found>}  young-regex=$YOUNG_RE"
echo "############################################################################"

declare -A V   # tool -> verdict line (last VERDICT/SKIP line captured)
run() {        # run <label> <cmd...>; capture the final VERDICT/SKIP line
  local label="$1"; shift
  echo; echo "----- [$label] -----"
  local log="$OUTDIR/${label}.log"
  "$@" 2>&1 | tee "$log"
  local v
  # capture the tool's last GREEN/RED/SKIP verdict line; de_precheck reports its summary as
  # "GREEN: ..." / "N association(s) ... FLAGGED", so match those too.
  v="$(grep -E '^[[:space:]]*(VERDICT|SKIP|GREEN:|FAIL)|association\(s\).*FLAGGED' "$log" | tail -1)"
  V["$label"]="${v:-(no verdict line)}"
}

# 04 — the engine (feeds the rest)
if [[ -n "$BAM" ]]; then
  run 04_core_regime_witness bash "$TOOLS/04_core_regime_witness.sh" \
    --bam "$BAM" --saf "$SAF" --outdir "$OUTDIR" --image "$IMAGE" --threads "$THREADS"
else
  echo; echo "----- [04_core_regime_witness] -----"
  echo "SKIP [04] no BAM found in $BAM_DIR (glob $BAM_GLOB) — downstream per-read tools will skip"
  V["04_core_regime_witness"]="SKIP [04] no BAM found"
fi

# derive sample + joint_counts + CORE paths produced by 04
if [[ -n "$BAM" ]]; then
  sample="$(basename "$BAM")"; sample="${sample%.markdup.sorted.bam}"; sample="${sample%.bam}"
else
  sample="${SAMPLE:-sample}"
fi
JOINT="$OUTDIR/joint_counts_${sample}.json"
core_of() { local d="$1" f="$OUTDIR/$d/$(basename "$BAM").featureCounts"; [[ -f "$f" ]] || f="${f}.gz"; echo "$f"; }
F0="$(core_of s0)"; F2="$(core_of s2)"; F1="$(core_of s1)"

# 05 — shared classifier (needs joint_counts from 04)
if [[ -f "$JOINT" ]]; then
  run 05_classify_triples python3 "$TOOLS/05_classify_triples.py" "$JOINT"
  run 06_closure_audit python3 "$TOOLS/06_closure_audit.py" "$JOINT" --sample "$sample" \
    --out "$OUTDIR/closure_report_${sample}.txt"
else
  echo; echo "----- [05_classify_triples] -----"; echo "SKIP [05] no joint_counts (04 skipped)"
  V["05_classify_triples"]="SKIP [05] no joint_counts (04 skipped)"
  echo; echo "----- [06_closure_audit] -----"; echo "SKIP [06] no joint_counts (04 skipped)"
  V["06_closure_audit"]="SKIP [06] no joint_counts (04 skipped)"
fi

# 02 — fragment-weighting (needs BAM + s2 CORE + matrix)
if [[ -n "$BAM" && -f "$F2" && -n "$MATRIX" && -f "$MATRIX" ]]; then
  run 02_weighting_check bash "$TOOLS/02_weighting_check.sh" \
    --bam "$BAM" --core "$F2" --matrix "$MATRIX" --outdir "$OUTDIR" --young-regex "$YOUNG_RE"
else
  echo; echo "----- [02_weighting_check] -----"
  echo "SKIP [02] need --matrix (the s2 count matrix) + the s2 CORE from 04 + a BAM"
  V["02_weighting_check"]="SKIP [02] need --matrix + s2 CORE + BAM"
fi

# 01 — Random-One BAM witness (needs BAM)
if [[ -n "$BAM" ]]; then
  run 01_random_one_check bash "$TOOLS/01_random_one_check.sh" --bam "$BAM" --outdir "$OUTDIR"
else
  echo; echo "----- [01_random_one_check] -----"; echo "SKIP [01] no BAM"
  V["01_random_one_check"]="SKIP [01] no BAM"
fi

# 07 — -O silent attribution (.sh re-run, then .py lookup). needs 04 CORE files.
if [[ -n "$BAM" && -f "$F0" ]]; then
  run 07_silent_attribution_sh bash "$TOOLS/07_silent_attribution.sh" \
    --bam "$BAM" --saf "$SAF" --outdir "$OUTDIR" --image "$IMAGE" --threads "$THREADS"
  sOcore="$OUTDIR/sO/$(basename "$BAM").featureCounts"; [[ -f "$sOcore" ]] || sOcore="${sOcore}.gz"
  if [[ -f "$sOcore" ]]; then
    run 07_silent_attribution_py python3 "$TOOLS/07_silent_attribution.py" \
      --s0-core "$F0" --s2-core "$F2" --s1-core "$F1" --sO-core "$sOcore" \
      --outdir "$OUTDIR" --young-regex "$YOUNG_RE"
  else
    echo; echo "----- [07_silent_attribution_py] -----"; echo "SKIP [07.py] no -O CORE (07.sh skipped)"
    V["07_silent_attribution_py"]="SKIP [07.py] no -O CORE (07.sh skipped)"
  fi
else
  echo; echo "----- [07_silent_attribution] -----"; echo "SKIP [07] need 04 CORE files + a BAM"
  V["07_silent_attribution_sh"]="SKIP [07] need 04 CORE + BAM"
fi

# 08 — young assignable-evidence gate (needs 04 CORE files)
if [[ -f "$F0" && -f "$F2" && -f "$F1" ]]; then
  SILENT_SUMMARY="$OUTDIR/young_silent_summary.tsv"
  run 08_young_gate python3 "$TOOLS/08_young_gate.py" \
    --s0-core "$F0" --s2-core "$F2" --s1-core "$F1" --outdir "$OUTDIR" \
    --young-regex "$YOUNG_RE" $([[ -f "$SILENT_SUMMARY" ]] && echo --silent-summary "$SILENT_SUMMARY")
else
  echo; echo "----- [08_young_gate] -----"; echo "SKIP [08] no 04 CORE files"
  V["08_young_gate"]="SKIP [08] no 04 CORE files"
fi

# 09 — SAF self-intersect geometry (needs bedtools; skips gracefully internally)
WITNESS_ARG=""
[[ -f "$OUTDIR/young_silent_summary.tsv" ]] && WITNESS_ARG="--witness $OUTDIR/young_silent_summary.tsv"
# shellcheck disable=SC2086
run 09_geometry_overlap python3 "$TOOLS/09_geometry_overlap.py" \
  --saf "$SAF" --outdir "$OUTDIR" --young-regex "$YOUNG_RE" $WITNESS_ARG

# 03 — strand-invariant + per-sample meter (needs the four .summary files)
if [[ -n "$S0_SUM" && -n "$S2_SUM" && -n "$S1_SUM" && -n "$GENE_SUM" ]]; then
  MAT_ARGS=""
  [[ -n "$S0_MAT" && -n "$S2_MAT" && -n "$S1_MAT" ]] && \
    MAT_ARGS="--s0-matrix $S0_MAT --s2-matrix $S2_MAT --s1-matrix $S1_MAT"
  # shellcheck disable=SC2086
  run 03_strand_invariant python3 "$TOOLS/03_strand_invariant.py" \
    --s0 "$S0_SUM" --s2 "$S2_SUM" --s1 "$S1_SUM" --gene "$GENE_SUM" \
    --outdir "$OUTDIR" $MAT_ARGS
else
  echo; echo "----- [03_strand_invariant] -----"
  echo "SKIP [03] need --s0-summary --s2-summary --s1-summary --gene-summary"
  V["03_strand_invariant"]="SKIP [03] need the four .summary files"
fi

# de_precheck — 3-metric condition-linkage (needs a design matrix + the 03 QC table)
if [[ -n "$DESIGN" ]]; then
  QC_TSV="$OUTDIR/per_sample_strand_qc.tsv"
  if [[ -f "$QC_TSV" && -n "$GENE_SUM" ]]; then
    python3 "$HERE/de_precheck/build_per_sample_qc.py" --qc-tsv "$QC_TSV" --gene-summary "$GENE_SUM" \
      2>&1 | tee "$OUTDIR/de_precheck_build.log" || true
  fi
  if [[ -f "$QC_TSV" ]]; then
    run de_precheck python3 "$HERE/de_precheck/check_silentloss_vs_design.py" "$DESIGN" \
      --qc "$QC_TSV" --out "$OUTDIR/silentloss_confound_result.tsv"
  else
    echo; echo "----- [de_precheck] -----"; echo "SKIP [de_precheck] no per_sample_strand_qc.tsv (run 03 first)"
    V["de_precheck"]="SKIP [de_precheck] no per_sample_strand_qc.tsv"
  fi
fi

# ---------------------------------------------------------------------------
# consolidated verdict block
# ---------------------------------------------------------------------------
echo
echo "############################################################################"
echo "# CONSOLIDATED VERDICT — strand-split TE QC suite"
echo "#   strand=$STRAND  sample=$sample  outdir=$OUTDIR"
echo "############################################################################"
order=(01_random_one_check 02_weighting_check 03_strand_invariant 04_core_regime_witness \
       05_classify_triples 06_closure_audit 07_silent_attribution_sh 07_silent_attribution_py \
       08_young_gate 09_geometry_overlap de_precheck)
nred=0
for k in "${order[@]}"; do
  line="${V[$k]:-}"
  [[ -z "$line" ]] && continue
  printf "  %-26s %s\n" "$k" "$line"
  # RED tool, or a de_precheck design association -> count as flagged
  if [[ "$line" == *"[RED]"* || "$line" == *"FLAGGED"* ]]; then
    nred=$((nred+1))
  fi
done
echo "----------------------------------------------------------------------------"
if [[ "$nred" -eq 0 ]]; then
  echo "  SUITE VERDICT: GREEN (no RED tool; SKIPs indicate an absent dependency, not a failure)"
  echo "  Flags to carry forward: FLAG-METER-NOT-ESTIMATOR, FLAG-RESIDUAL-NOT-LOSS,"
  echo "    FLAG-ALIGNMENT-WEIGHTED, FLAG-DENOM-SCALE, FLAG-SILENTLOSS-DESIGN (see references/strand-split-qc.md + docs/QC.md)"
  exit 0
else
  echo "  SUITE VERDICT: RED ($nred tool(s) flagged) — inspect the per-tool logs in $OUTDIR"
  exit 1
fi
