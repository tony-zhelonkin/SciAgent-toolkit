#!/usr/bin/env bash
# 02_weighting_check.sh — is the SENSE integer matrix fragment- or alignment-weighted?
#
# Generalized from kernel_trace/weighting_check.sh + weighting_check_0033.sh (hard-coded BAM/CORE/
# matrix/sample stripped; now via flags). Operationalizes lesson 6b + the counterfactual-risk
# pattern (KERNEL_TRACE / GAP-D5).
#
# For each fragment Assigned in the s2 (sense) CORE pass (n_targets always 1) we know its single
# GeneID. Join to the TRUE multiplicity NH from the genome BAM. Per subfamily:
#   nfrag  = distinct Assigned fragments  (fragment-weighted prediction)
#   sumNH  = sum of NH over those frags   (alignment-weighted prediction)
#   matrix = the integer count in the s2 matrix
# The matrix is fragment-weighted iff matrix == nfrag (and matrix << sumNH for high-NH young
# families). meanNH over the young set = the fold an --outSAMmultNmax all + -M run WOULD inflate by.
#
# Threshold: PASS iff matrix==nfrag for every subfamily (fragment-weighted). RED if matrix!=nfrag —
# report the realized matrix/sumNH and the would-be young meanNH fold (FLAG-ALIGNMENT-WEIGHTED).
#
# samtools is REQUIRED. SKIPS GRACEFULLY (exit 0) if samtools absent or any input missing.
#
# Usage:
#   02_weighting_check.sh --bam <BAM> --core <s2 CORE.gz> --matrix <te_sense_s2.counts.txt> \
#     --outdir <OUT> [--young-regex '^(L1MdT|L1MdGf|L1MdA|IAPEz)'] [--samtools samtools]
set -u

NAME="02_weighting_check"
SAMTOOLS="samtools"
BAM=""; CORE=""; MATRIX=""; OUT=""
YOUNG_RE='^(L1MdT|L1MdGf|L1MdA|IAPEz)'

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam)         BAM="$2"; shift 2 ;;
    --core)        CORE="$2"; shift 2 ;;
    --matrix)      MATRIX="$2"; shift 2 ;;
    --outdir)      OUT="$2"; shift 2 ;;
    --young-regex) YOUNG_RE="$2"; shift 2 ;;
    --samtools)    SAMTOOLS="$2"; shift 2 ;;
    -h|--help)     grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "[$NAME] unknown arg: $1" >&2; exit 2 ;;
  esac
done

[[ -n "$BAM" && -n "$CORE" && -n "$MATRIX" && -n "$OUT" ]] || {
  echo "[$NAME] need --bam --core --matrix --outdir" >&2; exit 2; }
command -v "$SAMTOOLS" >/dev/null 2>&1 || {
  echo "SKIP [$NAME] samtools ('$SAMTOOLS') not on PATH — skipping weighting check"; exit 0; }
for f in "$BAM" "$CORE" "$MATRIX"; do
  [[ -f "$f" ]] || { echo "SKIP [$NAME] input not found: $f — skipping"; exit 0; }
done

mkdir -p "$OUT"
sample="$(basename "$BAM")"; sample="${sample%.markdup.sorted.bam}"; sample="${sample%.bam}"
NHF="$OUT/frag_nh_${sample}.tsv"
PERSUB="$OUT/persub_nfrag_sumnh_${sample}.tsv"
MCOUNT="$OUT/matrix_count_${sample}.tsv"
JOIN="$OUT/weighting_join_${sample}.tsv"

# how to read a possibly-gzipped CORE
zcore() { case "$CORE" in *.gz) zcat "$CORE" ;; *) cat "$CORE" ;; esac; }

echo "[$NAME] [1/3] fragment->NH from BAM ..." >&2
"$SAMTOOLS" view "$BAM" \
| awk '{nh=0;for(i=12;i<=NF;i++)if($i~/^NH:i:/){nh=substr($i,6)+0;break}
        if(!($1 in seen)){seen[$1]=1; print $1"\t"nh}}' > "$NHF"

echo "[$NAME] [2/3] join CORE-assigned frags to NH, aggregate per subfamily ..." >&2
zcore \
| awk -F'\t' -v NHF="$NHF" '
  BEGIN{ while((getline l < NHF)>0){ split(l,a,"\t"); nh[a[1]]=a[2]+0 } }
  $2=="Assigned"{
     gid=$4; n=(($1 in nh)? nh[$1] : 1)
     nfrag[gid]++; sumnh[gid]+=n; if(n>maxnh[gid]) maxnh[gid]=n
  }
  END{ for(g in nfrag) printf "%s\t%d\t%d\t%d\n", g, nfrag[g], sumnh[g], maxnh[g] }' \
| sort -t$'\t' -k1,1 > "$PERSUB"

echo "[$NAME] [3/3] join with matrix counts ..." >&2
awk -F'\t' 'NR>2{print $1"\t"$NF}' "$MATRIX" | sort -t$'\t' -k1,1 > "$MCOUNT"
join -t$'\t' -1 1 -2 1 "$PERSUB" "$MCOUNT" > "$JOIN"
# columns: GeneID  nfrag  sumNH  maxNH  matrix
echo "[$NAME] joined rows: $(wc -l < "$JOIN")" >&2

# global match + young counterfactual
awk -F'\t' -v YRE="$YOUNG_RE" '
  function isyoung(g,   sf){ sf=g; sub(/:.*/,"",sf); return (sf ~ YRE) }
  {
    eqsum = ($5==$2)
    if(eqsum) eq++; else { ne++; if(ne<=10) print "  MISMATCH",$0 }
    tot++
    if(isyoung($1)){ ynf+=$2; ynh+=$3; ymx=($4>ymx?$4:ymx); ymtx+=$5 }
  }
  END{
    printf "matrix==nfrag: %d   mismatch: %d   total subfamilies: %d\n", eq, ne, tot
    if(ynf>0){
      printf "YOUNG set: nfrag=%d sumNH=%d matrix=%d  meanNH=%.2f (counterfactual alignment-weight fold)  maxNH=%d\n", \
        ynf, ynh, ymtx, ynh/ynf, ymx
    }
    print (ne==0 ? "MATCH=1" : "MATCH=0")
  }' "$JOIN" | tee "$OUT/weighting_verdict_${sample}.txt"

match="$(awk -F= '/^MATCH=/{print $2}' "$OUT/weighting_verdict_${sample}.txt")"
echo
if [[ "${match:-0}" -eq 1 ]]; then
  echo "VERDICT [GREEN] [$NAME] matrix==nfrag for every subfamily — SENSE matrix is fragment-weighted;" \
       "young families NOT alignment-inflated (the meanNH fold is the averted counterfactual)"
  exit 0
else
  echo "VERDICT [RED] [$NAME] matrix!=nfrag for some subfamily — FLAG-ALIGNMENT-WEIGHTED:" \
       "young high-NH families may be alignment-inflated up to ~meanNH-fold; re-verify the kernel"
  exit 1
fi
