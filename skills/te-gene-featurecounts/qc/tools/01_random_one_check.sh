#!/usr/bin/env bash
# 01_random_one_check.sh — BAM-witness that the genome BAM is Random-One.
#
# Generalized from kernel_trace/confirm_random_one.sh (hard-coded BAM path & sample stripped;
# now --bam <BAM> --outdir <OUT>). Operationalizes lesson 6a (KERNEL_TRACE).
#
# Under STAR --outSAMmultNmax 1 each multimapping fragment emits exactly ONE alignment locus
# (one R1 + one R2, both HI:i:1) while NH records the TRUE multiplicity (>1). Restricted to
# NH>1 fragments (the only ones that COULD emit >1 locus) we count, per readname:
#   (a) BAM lines (should be 2: R1+R2 of one locus)
#   (b) distinct emitted HI values (=1 under Random-One)
#   (c) secondary(0x100) lines (~0 — only one hit emitted)
#   (d) distinct (pos,mate) loci (<=2)
#
# Threshold: PASS iff >1-locus fragments == 0 (Random-One holds). RED if any fragment emits
# >1 locus — young high-NH families may be alignment-inflated (FLAG-ALIGNMENT-WEIGHTED).
#
# samtools is REQUIRED. If samtools is not on PATH this tool SKIPS GRACEFULLY (exit 0),
# mirroring the container-tool idiom of tests/run_regression.sh.
#
# Usage: 01_random_one_check.sh --bam <BAM> --outdir <OUT> [--samtools samtools]
set -u

NAME="01_random_one_check"
SAMTOOLS="samtools"
BAM=""
OUT=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam)      BAM="$2"; shift 2 ;;
    --outdir)   OUT="$2"; shift 2 ;;
    --samtools) SAMTOOLS="$2"; shift 2 ;;
    -h|--help)  grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "[$NAME] unknown arg: $1" >&2; exit 2 ;;
  esac
done

[[ -n "$BAM" && -n "$OUT" ]] || { echo "[$NAME] need --bam and --outdir" >&2; exit 2; }
command -v "$SAMTOOLS" >/dev/null 2>&1 || {
  echo "SKIP [$NAME] samtools ('$SAMTOOLS') not on PATH — skipping Random-One BAM witness"; exit 0; }
[[ -f "$BAM" ]] || { echo "SKIP [$NAME] BAM not found: $BAM — skipping"; exit 0; }

mkdir -p "$OUT"
sample="$(basename "$BAM")"; sample="${sample%.markdup.sorted.bam}"; sample="${sample%.bam}"
report="$OUT/random_one_${sample}.txt"

"$SAMTOOLS" view "$BAM" \
| awk -v OFS='\t' '
{
  nh=0; hi="-";
  for(i=12;i<=NF;i++){
    if($i ~ /^NH:i:/){ nh=substr($i,6)+0 }
    else if($i ~ /^HI:i:/){ hi=substr($i,6) }
  }
  if(nh<=1) next;            # only multimappers can possibly emit >1 locus
  flag=$2+0; rn=$1
  lines[rn]++
  key=rn SUBSEP hi
  if(!(key in seenhi)){ seenhi[key]=1; nhi[rn]++ }
  if(and(flag,256)) sec[rn]++
  mate = and(flag,64) ? "R1" : (and(flag,128) ? "R2" : "U")
  lockey = rn SUBSEP $3 SUBSEP $4 SUBSEP mate
  if(!(lockey in seenloc)){ seenloc[lockey]=1; nloc[rn]++ }
}
END{
  total=0; lines2=0; linesNot2=0; hi1=0; hiGt1=0; secAny=0; locOK=0; locBad=0
  for(rn in lines){
    total++
    if(lines[rn]==2) lines2++; else linesNot2++
    if(nhi[rn]==1) hi1++; else hiGt1++
    if((rn in sec)&&sec[rn]>0) secAny++
    if(nloc[rn]<=2) locOK++; else locBad++
  }
  if(total==0){ print "NH>1 fragments: 0 (no multimappers to test)"; print "MULTILOCUS=0"; exit }
  printf "NH>1 fragments (distinct readnames)        : %d\n", total
  printf "  with exactly 2 BAM lines (R1+R2)         : %d (%.4f%%)\n", lines2, 100.0*lines2/total
  printf "  with >1 distinct HI (>1 emitted locus)   : %d (%.4f%%)\n", hiGt1, 100.0*hiGt1/total
  printf "  with any secondary(0x100) line           : %d (%.4f%%)\n", secAny, 100.0*secAny/total
  printf "  with >2 distinct (pos,mate) loci (BAD)   : %d (%.4f%%)\n", locBad, 100.0*locBad/total
  printf "MULTILOCUS=%d\n", hiGt1+locBad
}' | tee "$report"

multilocus="$(awk -F= '/^MULTILOCUS=/{print $2}' "$report")"
multilocus="${multilocus:-0}"
echo
if [[ "$multilocus" -eq 0 ]]; then
  echo "VERDICT [GREEN] [$NAME] Random-One holds: 0 fragments emit >1 locus"
  exit 0
else
  echo "VERDICT [RED] [$NAME] $multilocus fragments emit >1 locus — FLAG-ALIGNMENT-WEIGHTED:" \
       "young high-NH families may be alignment-inflated; re-verify before trusting young fold-changes"
  exit 1
fi
