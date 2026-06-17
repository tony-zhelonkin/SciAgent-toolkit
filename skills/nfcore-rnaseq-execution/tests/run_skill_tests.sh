#!/usr/bin/env bash
# Skill-local test entry point, auto-discovered by the toolkit's tests/run-all.sh sweep.
# Tests scripts/make_samplesheet.sh on a synthetic fixture (empty FASTQ files named per the
# Illumina convention, including a 2-lane sample). Skips GRACEFULLY (exit 0 with a notice) if
# a required tool is missing, so the bash CI never blocks on an un-bootstrapped machine.
set -euo pipefail

SKILL="$(cd "$(dirname "$0")/.." && pwd)"
GEN="$SKILL/scripts/make_samplesheet.sh"

# Graceful skip if a needed tool is absent (mirrors the mllmcelltype pattern).
for tool in bash sed sort basename mktemp; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "SKIP [nfcore-rnaseq-execution]: '$tool' not on PATH"; exit 0
  fi
done
if [ ! -x "$GEN" ] && [ ! -f "$GEN" ]; then
  echo "SKIP [nfcore-rnaseq-execution]: $GEN absent"; exit 0
fi

fail=0
note() { echo "  $1"; }
check() { # check <desc> <expected> <actual>
  if [ "$2" = "$3" ]; then note "PASS: $1"; else note "FAIL: $1 (expected [$2], got [$3])"; fail=1; fi
}

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
FQ="$TMP/data"
mkdir -p "$FQ"

# Fixture: SAMPLEA single lane; SAMPLEB spans two lanes (L005 + L006); SAMPLEC single lane.
touch "$FQ/SAMPLEA_S1_L005_R1_001.fastq.gz" "$FQ/SAMPLEA_S1_L005_R2_001.fastq.gz"
touch "$FQ/SAMPLEB_S2_L005_R1_001.fastq.gz" "$FQ/SAMPLEB_S2_L005_R2_001.fastq.gz"
touch "$FQ/SAMPLEB_S2_L006_R1_001.fastq.gz" "$FQ/SAMPLEB_S2_L006_R2_001.fastq.gz"
touch "$FQ/SAMPLEC_S3_L005_R1_001.fastq.gz" "$FQ/SAMPLEC_S3_L005_R2_001.fastq.gz"

echo "[nfcore-rnaseq-execution] testing make_samplesheet.sh"

OUT="$(bash "$GEN" "$FQ" reverse)"

# 1. Header
header="$(printf '%s\n' "$OUT" | head -n1)"
check "header" "sample,fastq_1,fastq_2,strandedness" "$header"

# 2. Row count: 4 FASTQ-pairs -> 4 data rows (+1 header = 5 lines)
nlines="$(printf '%s\n' "$OUT" | grep -c .)"
check "total lines (header + 4 rows)" "5" "$nlines"

# 3. Sample-name derivation: distinct samples = SAMPLEA, SAMPLEB, SAMPLEC
samples="$(printf '%s\n' "$OUT" | tail -n +2 | cut -d, -f1 | sort -u | paste -sd, -)"
check "distinct samples" "SAMPLEA,SAMPLEB,SAMPLEC" "$samples"

# 4. Multi-lane handling: SAMPLEB appears on TWO rows (one per lane)
b_rows="$(printf '%s\n' "$OUT" | tail -n +2 | cut -d, -f1 | grep -c '^SAMPLEB$')"
check "SAMPLEB row count (per-lane)" "2" "$b_rows"

# 5. R1/R2 pairing on SAMPLEA row: fastq_1 ends R1, fastq_2 ends R2, both absolute
a_row="$(printf '%s\n' "$OUT" | grep '^SAMPLEA,')"
f1="$(printf '%s\n' "$a_row" | cut -d, -f2)"
f2="$(printf '%s\n' "$a_row" | cut -d, -f3)"
case "$f1" in /*_R1_001.fastq.gz) p1=ok;; *) p1=bad;; esac
case "$f2" in /*_R2_001.fastq.gz) p2=ok;; *) p2=bad;; esac
check "fastq_1 absolute & R1" "ok" "$p1"
check "fastq_2 absolute & R2" "ok" "$p2"

# 6. Strandedness column carries the CLI arg
strand="$(printf '%s\n' "$a_row" | cut -d, -f4)"
check "strandedness arg" "reverse" "$strand"

# 7. Default strandedness is `auto`
def_strand="$(bash "$GEN" "$FQ" | grep '^SAMPLEA,' | cut -d, -f4)"
check "default strandedness" "auto" "$def_strand"

# 7b. --seq-platform appends a seq_platform column to header and rows
SP_OUT="$(bash "$GEN" "$FQ" auto --seq-platform ILLUMINA)"
sp_header="$(printf '%s\n' "$SP_OUT" | head -n1)"
check "seq_platform header" "sample,fastq_1,fastq_2,strandedness,seq_platform" "$sp_header"
sp_val="$(printf '%s\n' "$SP_OUT" | grep '^SAMPLEA,' | cut -d, -f5)"
check "seq_platform value" "ILLUMINA" "$sp_val"
# default (no flag) keeps the 4-column header
DEF_OUT="$(bash "$GEN" "$FQ")"
def_header="$(printf '%s\n' "$DEF_OUT" | head -n1)"
check "no seq_platform by default" "sample,fastq_1,fastq_2,strandedness" "$def_header"

# 8. Error on missing R2
ORPH="$TMP/orphan"; mkdir -p "$ORPH"
touch "$ORPH/SAMPLED_S1_L005_R1_001.fastq.gz"  # no R2
if bash "$GEN" "$ORPH" >/dev/null 2>&1; then
  note "FAIL: expected nonzero exit on missing R2"; fail=1
else
  note "PASS: errors on missing R2"
fi

# 9. Error on missing directory
if bash "$GEN" "$TMP/does_not_exist" >/dev/null 2>&1; then
  note "FAIL: expected nonzero exit on missing dir"; fail=1
else
  note "PASS: errors on missing dir"
fi

echo
if [ "$fail" -eq 0 ]; then
  echo "[nfcore-rnaseq-execution] ALL TESTS PASS"
  exit 0
else
  echo "[nfcore-rnaseq-execution] TESTS FAILED"
  exit 1
fi
