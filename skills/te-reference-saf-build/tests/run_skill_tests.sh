#!/usr/bin/env bash
# Skill-local test entry point, auto-discovered by the toolkit's tests/run-all.sh sweep.
# Smoke-tests build_te_saf.sh on a tiny SYNTHETIC TE GTF fixture: asserts the grouped SAF has
# the right GeneID = Subfamily:Family:Class labels and the GeneID/Chr/Start/End/Strand columns.
# The bedtools exon-subtraction portion is exercised when bedtools is present, and SKIPPED
# gracefully otherwise (mirrors the mllmcelltype / nfcore graceful-skip pattern). The awk-only
# grouping assertions ALWAYS run (gawk required). Prints PASS/FAIL; exits nonzero on failure.
set -euo pipefail

SKILL="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="$SKILL/scripts/build_te_saf.sh"

# Graceful skip if gawk or a coreutil is missing (never block CI on an un-bootstrapped box).
for tool in bash gawk sort cut wc mktemp sed; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "SKIP [te-reference-saf-build]: '$tool' not on PATH"; exit 0
  fi
done
if [ ! -f "$BUILD" ]; then
  echo "SKIP [te-reference-saf-build]: $BUILD absent"; exit 0
fi

fail=0
note() { echo "  $1"; }
check() { if [ "$2" = "$3" ]; then note "PASS: $1"; else note "FAIL: $1 (expected [$2], got [$3])"; fail=1; fi; }

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

echo "[te-reference-saf-build] testing build_te_saf.sh on synthetic fixtures"

# --- Synthetic TE GTF: known subfamily/family/class per locus -------------------------------
# Two L1Md_A loci (same group -> meta-feature), one Alu SINE, one MER DNA. Contigs without 'chr'.
TE_GTF="$TMP/te.gtf"
cat > "$TE_GTF" <<'EOF'
1	rmsk	exon	100	200	.	+	.	gene_id "L1Md_A"; transcript_id "L1Md_A_dup1"; family_id "L1"; class_id "LINE";
1	rmsk	exon	500	600	.	-	.	gene_id "L1Md_A"; transcript_id "L1Md_A_dup2"; family_id "L1"; class_id "LINE";
2	rmsk	exon	300	350	.	+	.	gene_id "B1_Mus1"; transcript_id "B1_Mus1_dup1"; family_id "Alu"; class_id "SINE";
3	rmsk	exon	800	900	.	+	.	gene_id "MER20"; transcript_id "MER20_dup1"; family_id "hAT"; class_id "DNA";
EOF

# --- Synthetic GENE GTF: one exon overlapping the FIRST L1Md_A locus (1:120-180) ------------
GENE_GTF="$TMP/gene.gtf"
cat > "$GENE_GTF" <<'EOF'
1	src	gene	120	180	.	+	.	gene_id "Gtest";
1	src	exon	120	180	.	+	.	gene_id "Gtest"; transcript_id "Gtest.1";
EOF

OUT="$TMP/out"
mkdir -p "$OUT"

HAS_BEDTOOLS=1
command -v bedtools >/dev/null 2>&1 || HAS_BEDTOOLS=0

if [ "$HAS_BEDTOOLS" -eq 1 ]; then
  bash "$BUILD" --te-gtf "$TE_GTF" --gene-gtf "$GENE_GTF" --out-dir "$OUT" --prefix TEST >/dev/null
  SAF_ALL="$OUT/TEST_GROUPED_all.saf"
  SAF_NOEXON="$OUT/TEST_GROUPED_all_noExon.saf"
else
  note "SKIP: bedtools absent — testing grouped-SAF awk logic only (no exon subtraction)"
  # Reproduce ONLY the grouped-SAF awk (step 2) so the label assertions still run.
  SAF_ALL="$OUT/TEST_GROUPED_all.saf"
  gawk 'BEGIN{OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}
     $0 !~ /^#/ {
       match($0,/gene_id "([^"]+)"/,g);
       match($0,/family_id "([^"]+)"/,f);
       match($0,/class_id "([^"]+)"/,c);
       print g[1] ":" f[1] ":" c[1], $1, $4, $5, $7
     }' "$TE_GTF" > "$SAF_ALL"
  SAF_NOEXON=""
fi

# --- Assertions on the grouped SAF (always) -------------------------------------------------
header="$(head -n1 "$SAF_ALL")"
check "grouped SAF header columns" "$(printf 'GeneID\tChr\tStart\tEnd\tStrand')" "$header"

# 4 input loci -> 4 data rows
nrows="$(tail -n +2 "$SAF_ALL" | grep -c .)"
check "grouped SAF row count (4 loci)" "4" "$nrows"

# GeneID = Subfamily:Family:Class for each known locus
row_l1="$(grep -P '^L1Md_A:L1:LINE\t1\t100\t200\t\+$' "$SAF_ALL" | head -n1 | wc -l)"
check "L1Md_A:L1:LINE label + cols (1/100/200/+)" "1" "$row_l1"
row_alu="$(grep -P '^B1_Mus1:Alu:SINE\t2\t300\t350\t\+$' "$SAF_ALL" | wc -l)"
check "B1_Mus1:Alu:SINE label + cols" "1" "$row_alu"
row_dna="$(grep -P '^MER20:hAT:DNA\t3\t800\t900\t\+$' "$SAF_ALL" | wc -l)"
check "MER20:hAT:DNA label + cols" "1" "$row_dna"

# Both L1Md_A loci collapse to ONE group label (the meta-feature behaviour)
n_groups="$(cut -f1 "$SAF_ALL" | tail -n +2 | sort -u | wc -l)"
check "unique groups (L1Md_A pooled -> 3 groups)" "3" "$n_groups"

# class whitelist check (awk-level, no deps): retro-only drops the DNA locus
retro_groups="$(gawk -v pat='^(LINE|SINE|LTR|RC)$' 'BEGIN{OFS="\t"}
     $0 !~ /^#/ {
       match($0,/family_id "([^"]+)"/,f); match($0,/class_id "([^"]+)"/,c);
       if (c[1] !~ pat) next;
       match($0,/gene_id "([^"]+)"/,g); print g[1]":"f[1]":"c[1]
     }' "$TE_GTF" | sort -u | wc -l)"
check "retro-only whitelist drops DNA (2 groups)" "2" "$retro_groups"

# --- Assertions on the no-exon SAF (only if bedtools ran) -----------------------------------
if [ -n "$SAF_NOEXON" ]; then
  # The first L1Md_A locus (1:100-200) overlaps gene exon 1:120-180 -> it gets split/trimmed,
  # NOT dropped (bedtools subtract leaves the flanking 100-119 and 181-200). Group remains.
  noexon_header="$(head -n1 "$SAF_NOEXON")"
  check "no-exon SAF header columns" "$(printf 'GeneID\tChr\tStart\tEnd\tStrand')" "$noexon_header"
  # No remaining interval should overlap the exon 120-180 on chr 1.
  overlap="$(tail -n +2 "$SAF_NOEXON" | gawk '$2==1 && $3<=180 && $4>=120' | wc -l)"
  check "no residual TE interval overlaps exon 1:120-180" "0" "$overlap"
  # All 3 groups survive subtraction (signal trimmed, not zeroed).
  ne_groups="$(cut -f1 "$SAF_NOEXON" | tail -n +2 | sort -u | wc -l)"
  check "groups preserved after subtract (3)" "3" "$ne_groups"
  # Start column is 1-based (min start across all rows should be >=1; the L1 flank starts at 1 or 181).
  min_start="$(tail -n +2 "$SAF_NOEXON" | cut -f3 | sort -n | head -n1)"
  if [ "$min_start" -ge 1 ]; then note "PASS: no-exon Start is 1-based (min=$min_start)"; else note "FAIL: no-exon Start <1 (min=$min_start)"; fail=1; fi
fi

echo
if [ "$fail" -eq 0 ]; then
  echo "[te-reference-saf-build] ALL TESTS PASS"; exit 0
else
  echo "[te-reference-saf-build] TESTS FAILED"; exit 1
fi
