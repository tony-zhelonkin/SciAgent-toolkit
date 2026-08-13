#!/usr/bin/env bash
# run_te_counting.sh — parameterized wrapper for the TE + gene featureCounts
# counting workflow, run inside the locked te-fc:2.0.2 container.
#
# Generalizes the proven 14839-DM wrapper. It:
#   1) Stages the project BAMs as SYMLINKS into a bam_fin/ dir (no copies).
#   2) Bind-mounts the symlink-TARGET dir at an IDENTICAL host=container path
#      (-v $REAL:$REAL) so the links resolve INSIDE the container, and so the
#      driver's awk parser keys output column names off the BAM basename cleanly.
#   3) Runs runFeatureCounts_TE_and_genes.sh (TE pass + gene pass + combined)
#      in te-fc:2.0.2.
#
# IMPORTANT: gene strandedness (-s) is per-library. VERIFY it (MultiQC inferred
# strandedness + the featureCounts header) for every dataset — never hardcode.
# TE pass is fixed at -s 0 (unstranded), -M, integer Random-One (no --fraction).
#
# Usage:
#   run_te_counting.sh \
#     --bam-dir   DIR    # dir containing the project's *.bam (e.g. star_salmon;
#                        #   the *.markdup.sorted.bam lean-run path is fine)
#     --gene-gtf  FILE   # nf-core filtered gene GTF
#     --te-saf    FILE   # grouped exon-subtracted TE SAF (Subfamily:Family:Class)
#     --gene-s    0|1|2  # GENE strandedness (VERIFY per library; 2=reverse)
#     --out-dir   DIR    # output base dir
#     [--threads  N]     # default 12
#     [--te-strand unstranded|sense_antisense]   # default unstranded
#     [--image    TAG]   # default te-fc:2.0.2
#     [--bam-glob GLOB]  # default '*.bam' (e.g. '*.markdup.sorted.bam')
set -euo pipefail

THREADS=12
TE_STRAND_MODE=unstranded
IMAGE=te-fc:2.0.2
BAM_GLOB='*.bam'
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() { sed -n '2,40p' "${BASH_SOURCE[0]}" >&2; exit 1; }

while (( "$#" )); do
  case "$1" in
    --bam-dir)    BAM_DIR="$2"; shift 2;;
    --gene-gtf)   GENE_GTF="$2"; shift 2;;
    --te-saf)     TE_SAF="$2"; shift 2;;
    --gene-s)     GENE_S="$2"; shift 2;;
    --out-dir)    OUT="$2"; shift 2;;
    --threads)    THREADS="$2"; shift 2;;
    --te-strand)  TE_STRAND_MODE="$2"; shift 2;;
    --image)      IMAGE="$2"; shift 2;;
    --bam-glob)   BAM_GLOB="$2"; shift 2;;
    -h|--help)    usage;;
    *) echo "Unknown arg: $1" >&2; usage;;
  esac
done

for v in BAM_DIR GENE_GTF TE_SAF GENE_S OUT; do
  [[ -z "${!v:-}" ]] && { echo "Missing required: --${v,,}" >&2; usage; }
done
command -v docker >/dev/null 2>&1 || { echo "ERROR: docker not found" >&2; exit 1; }
[[ -d "$BAM_DIR" ]]  || { echo "BAM_DIR not found: $BAM_DIR" >&2; exit 1; }
[[ -f "$GENE_GTF" ]] || { echo "GENE_GTF not found: $GENE_GTF" >&2; exit 1; }
[[ -f "$TE_SAF" ]]   || { echo "TE_SAF not found: $TE_SAF" >&2; exit 1; }
[[ "$GENE_S" =~ ^[012]$ ]] || { echo "--gene-s must be 0,1,2" >&2; exit 1; }
docker image inspect "$IMAGE" >/dev/null 2>&1 || {
  echo "ERROR: image '$IMAGE' not found. Build it: bash $HERE/../env/build.sh" >&2; exit 1; }

# Absolutize inputs.
BAM_DIR="$(cd "$BAM_DIR" && pwd)"
GENE_GTF="$(readlink -f "$GENE_GTF")"
TE_SAF="$(readlink -f "$TE_SAF")"
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"

# --- Stage BAMs as symlinks in a bam_fin/ dir; collect the real target dirs. ---
IN="$OUT/bam_fin"
mkdir -p "$IN"
shopt -s nullglob
declare -a REAL_DIRS=()
n=0
for b in "$BAM_DIR"/$BAM_GLOB; do
  [[ -e "$b" ]] || continue
  ln -sf "$(readlink -f "$b")" "$IN/$(basename "$b")"
  REAL_DIRS+=("$(dirname "$(readlink -f "$b")")")
  n=$((n+1))
done
shopt -u nullglob
(( n > 0 )) || { echo "No BAMs matching '$BAM_GLOB' in $BAM_DIR" >&2; exit 1; }
echo "Staged $n BAM symlink(s) into $IN"

# Unique real target dirs (mounted identical host=container so links resolve).
mapfile -t REAL_DIRS < <(printf '%s\n' "${REAL_DIRS[@]}" | sort -u)

# --- Build docker mount args (identical host=container paths). ---
MOUNTS=( -v "$IN":"$IN" -v "$OUT":"$OUT" -v "$HERE":"$HERE":ro )
for d in "${REAL_DIRS[@]}"; do MOUNTS+=( -v "$d":"$d":ro ); done
MOUNTS+=( -v "$(dirname "$GENE_GTF")":"$(dirname "$GENE_GTF")":ro )
MOUNTS+=( -v "$(dirname "$TE_SAF")":"$(dirname "$TE_SAF")":ro )

DRIVER="$HERE/runFeatureCounts_TE_and_genes.sh"

echo "Running counting in $IMAGE (gene -s $GENE_S, TE -s 0 -M integer, te-strand $TE_STRAND_MODE)..."
docker run --rm -u "$(id -u):$(id -g)" \
  "${MOUNTS[@]}" \
  -e TMPDIR="$OUT/tmp" \
  "$IMAGE" \
  bash -lc "mkdir -p '$OUT/tmp' && '$DRIVER' \
    -i '$IN' -o '$OUT' \
    -g '$GENE_GTF' \
    -e '$TE_SAF' \
    -S '$GENE_S' \
    -t '$THREADS' \
    --te-strand '$TE_STRAND_MODE'"

echo "Done. Outputs under: $OUT"
echo "  TE matrix:       $OUT/featurecounts_TE/te_counts_matrix.txt"
echo "  Gene matrix:     $OUT/fc_genes/count_matrices_fc/sorted_counts_matrix.txt"
echo "  Combined matrix: $OUT/combined_gene_TE_counts.tsv"
