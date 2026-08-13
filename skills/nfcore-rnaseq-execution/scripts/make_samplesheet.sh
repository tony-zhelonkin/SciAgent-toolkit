#!/usr/bin/env bash
# make_samplesheet.sh — build an nf-core/rnaseq samplesheet from a directory of
# paired Illumina FASTQs.
#
# Convention (real data, e.g. 14839-DM-0001_S1_L005_R1_001.fastq.gz):
#   <sample>_S<n>_L<lane>_R1_001.fastq.gz / ..._R2_001.fastq.gz
# The `sample` name is the basename with the trailing _S<n>_L<lane>_R<read>_001.fastq.gz
# suffix stripped. A sample spanning multiple lanes emits ONE ROW PER LANE — nf-core/rnaseq
# merges rows that share the same `sample` value (see 02-nfcore-operational-knowledge.md §F).
#
# Output: CSV on stdout with header  sample,fastq_1,fastq_2,strandedness
#   (+ a trailing seq_platform column if --seq-platform is given).
#   fastq columns are absolute paths; rows are sorted deterministically.
#
# Usage:
#   make_samplesheet.sh <fastq_dir> [strandedness] [--seq-platform PLATFORM]   > samplesheet.csv
#     strandedness defaults to `auto` (nf-core 3.x infers it). Set the VERIFIED
#     value once known — the downstream featureCounts -s must match it.
#     --seq-platform appends an optional seq_platform column (e.g. ILLUMINA),
#     written into the BAM read-group PL tag. Omit it to use the 4 required cols.
#
# Zero dependency (POSIX-ish bash + coreutils). Exits nonzero on any error.

set -euo pipefail

usage() {
  cat >&2 <<'EOF'
Usage: make_samplesheet.sh <fastq_dir> [strandedness] [--seq-platform PLATFORM]

  <fastq_dir>          directory containing paired FASTQs named
                       <sample>_S<n>_L<lane>_R1_001.fastq.gz / ..._R2_001.fastq.gz
  [strandedness]       value for the strandedness column (default: auto)
  --seq-platform PLAT  append a seq_platform column (e.g. ILLUMINA) -> BAM PL tag

Writes an nf-core/rnaseq samplesheet CSV to stdout:
  sample,fastq_1,fastq_2,strandedness[,seq_platform]
Multi-lane samples emit one row per lane (nf-core merges same-sample rows).
EOF
}

if [ "${1:-}" = "-h" ] || [ "${1:-}" = "--help" ]; then
  usage; exit 0
fi

# Parse args: positional <dir> [strandedness] plus an optional --seq-platform PLAT.
dir=""
strandedness=""
seq_platform=""
while [ "$#" -gt 0 ]; do
  case "$1" in
    --seq-platform)
      [ "$#" -ge 2 ] || { echo "ERROR: --seq-platform needs a value." >&2; usage; exit 2; }
      seq_platform="$2"; shift 2 ;;
    --seq-platform=*)
      seq_platform="${1#*=}"; shift ;;
    -*)
      echo "ERROR: unknown option: $1" >&2; usage; exit 2 ;;
    *)
      if [ -z "$dir" ]; then dir="$1"
      elif [ -z "$strandedness" ]; then strandedness="$1"
      else echo "ERROR: too many positional arguments." >&2; usage; exit 2
      fi
      shift ;;
  esac
done

if [ -z "$dir" ]; then
  echo "ERROR: <fastq_dir> is required." >&2
  usage; exit 2
fi
strandedness="${strandedness:-auto}"

if [ ! -d "$dir" ]; then
  echo "ERROR: not a directory: $dir" >&2
  exit 1
fi

# Resolve to an absolute path so fastq columns are absolute (no `cd`, keep cwd stable).
abs_dir="$(cd "$dir" && pwd)"

# Collect R1 files (sorted, deterministic). nullglob so a no-match expands to nothing.
shopt -s nullglob
r1_files=("$abs_dir"/*_R1_001.fastq.gz)
shopt -u nullglob

if [ "${#r1_files[@]}" -eq 0 ]; then
  echo "ERROR: no *_R1_001.fastq.gz files found in $abs_dir" >&2
  exit 1
fi

# Sort the R1 list for deterministic output.
IFS=$'\n' r1_files=($(printf '%s\n' "${r1_files[@]}" | sort))
unset IFS

# Build rows into an array, then sort + emit (header first).
rows=()
for r1 in "${r1_files[@]}"; do
  base="$(basename "$r1")"
  # mate path: swap R1 -> R2 in the filename
  r2_base="${base/_R1_001.fastq.gz/_R2_001.fastq.gz}"
  r2="$abs_dir/$r2_base"
  if [ ! -f "$r2" ]; then
    echo "ERROR: R1 has no matching R2: $r1 (expected $r2)" >&2
    exit 1
  fi
  # Derive sample: strip _S<n>_L<lane>_R<read>_001.fastq.gz
  # base e.g. 14839-DM-0001_S1_L005_R1_001.fastq.gz -> 14839-DM-0001
  sample="$(printf '%s\n' "$base" | sed -E 's/_S[0-9]+_L[0-9]+_R[12]_001\.fastq\.gz$//')"
  if [ "$sample" = "$base" ]; then
    echo "ERROR: filename does not match <sample>_S<n>_L<lane>_R1_001.fastq.gz: $base" >&2
    exit 1
  fi
  if [ -n "$seq_platform" ]; then
    rows+=("$sample,$r1,$r2,$strandedness,$seq_platform")
  else
    rows+=("$sample,$r1,$r2,$strandedness")
  fi
done

if [ -n "$seq_platform" ]; then
  printf 'sample,fastq_1,fastq_2,strandedness,seq_platform\n'
else
  printf 'sample,fastq_1,fastq_2,strandedness\n'
fi
printf '%s\n' "${rows[@]}" | sort
