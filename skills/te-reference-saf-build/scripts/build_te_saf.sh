#!/usr/bin/env bash
# build_te_saf.sh — build the shared TE reference SAFs (grouped + exon-subtracted no-exon).
#
# This is the OWNED, single-source-of-truth executable for the build-once TE reference
# recipe. It crystallizes the inline awk/bedtools that previously lived only in per-dataset
# /scratch READMEs (13036-DM, AdaW_eWAT_WL). See SKILL.md for the science and provenance.
#
# What it does (deterministic):
#   1) Grouped SAF from the TEtranscripts TE GTF: one row per locus, GeneID = Subfamily:Family:Class
#      (gene_id:family_id:class_id), columns  GeneID Chr Start End Strand  (1-based, TSV).
#   2) No-exon SAF: derive an exon BED FROM THE GIVEN GENE GTF (never a hardcoded sibling-project
#      path), normalize contig names, bedtools-subtract exons, verify 0 residual overlap, rebuild SAF.
#
# Inputs:
#   --te-gtf    TEtranscripts RepeatMasker GTF (e.g. GRCm39_Ensembl_rmsk_TE.gtf.gz); .gz or plain.
#   --gene-gtf  Gene GTF for the SAME build (e.g. gencode.vM37...annotation.gtf[.gz]); exon BED source.
#   --out-dir   Output directory for the SAFs.
#   --prefix    Output basename prefix (e.g. GRCm39_rmsk_TE). Produces:
#                 <prefix>_GROUPED_all.saf  and  <prefix>_GROUPED_all_noExon.saf
# Options:
#   --keep-classes REGEX   Whitelist class_id values (e.g. '^(LINE|SINE|LTR|RC)$') for a retro-only
#                          SAF. Default: keep ALL classes (the canonical decision).
#   --strip-chr yes|no     Strip a leading 'chr' from contig names on BOTH BEDs before subtract
#                          (default: yes — Ensembl mm39 has no 'chr' prefix; mismatched naming
#                          silently yields zero overlap).
#   -h | --help
#
# Dependencies: gawk (GNU awk; 3-arg match), bedtools, sort, zcat. No other deps.
#
# TODO (deferred toolkit pass): consider migrating this into TE-RNAseq-toolkit/scripts and pointing here.

set -euo pipefail

usage() { sed -n '2,40p' "$0" >&2; exit "${1:-1}"; }

TE_GTF=""; GENE_GTF=""; OUT_DIR=""; PREFIX=""
KEEP_CLASSES=""; STRIP_CHR="yes"

while [ "$#" -gt 0 ]; do
  case "$1" in
    --te-gtf)       TE_GTF="${2:-}"; shift 2;;
    --gene-gtf)     GENE_GTF="${2:-}"; shift 2;;
    --out-dir)      OUT_DIR="${2:-}"; shift 2;;
    --prefix)       PREFIX="${2:-}"; shift 2;;
    --keep-classes) KEEP_CLASSES="${2:-}"; shift 2;;
    --strip-chr)    STRIP_CHR="${2:-}"; shift 2;;
    -h|--help)      usage 0;;
    *) echo "ERROR: unknown argument: $1" >&2; usage 2;;
  esac
done

# Validate required inputs.
for v in TE_GTF GENE_GTF OUT_DIR PREFIX; do
  if [ -z "${!v}" ]; then echo "ERROR: --$(echo "$v" | tr 'A-Z_' 'a-z-') is required." >&2; usage 2; fi
done
[ -f "$TE_GTF" ]   || { echo "ERROR: --te-gtf not found: $TE_GTF" >&2; exit 1; }
[ -f "$GENE_GTF" ] || { echo "ERROR: --gene-gtf not found: $GENE_GTF" >&2; exit 1; }
case "$STRIP_CHR" in yes|no) ;; *) echo "ERROR: --strip-chr must be yes|no" >&2; exit 2;; esac

# Dependency checks.
AWK="$(command -v gawk || true)"
[ -n "$AWK" ] || { echo "ERROR: gawk (GNU awk) is required for 3-arg match()." >&2; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo "ERROR: bedtools is required." >&2; exit 1; }

mkdir -p "$OUT_DIR"
SAF_ALL="$OUT_DIR/${PREFIX}_GROUPED_all.saf"
SAF_NOEXON="$OUT_DIR/${PREFIX}_GROUPED_all_noExon.saf"

# Reader that transparently handles .gz or plain.
cat_gtf() { case "$1" in *.gz) zcat "$1";; *) cat "$1";; esac; }

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

# --- Step 2: grouped SAF (GeneID = Subfamily:Family:Class) ---------------------------------
echo "[build_te_saf] Step 2: grouped SAF -> $SAF_ALL"
cat_gtf "$TE_GTF" |
"$AWK" -v pat="$KEEP_CLASSES" 'BEGIN{OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}
     $0 !~ /^#/ {
       match($0,/gene_id "([^"]+)"/,g);      # subfamily / repName
       match($0,/family_id "([^"]+)"/,f);    # family
       match($0,/class_id "([^"]+)"/,c);     # class (LINE/LTR/SINE/DNA/RC/...)
       if (pat != "" && c[1] !~ pat) next;   # optional class whitelist
       gid = g[1] ":" f[1] ":" c[1];         # TEtranscripts-like label
       print gid, $1, $4, $5, $7
     }' > "$SAF_ALL"

n_rows=$(( $(wc -l < "$SAF_ALL") - 1 ))
n_groups=$(cut -f1 "$SAF_ALL" | tail -n +2 | sort -u | wc -l)
echo "[build_te_saf]   grouped SAF: $n_rows loci, $n_groups unique groups"

# --- Step 3: no-exon SAF via bedtools subtract --------------------------------------------
echo "[build_te_saf] Step 3: exon subtraction -> $SAF_NOEXON"

# 3a) exon BED (0-based) from the GENE GTF (derived freshly here; no sibling-project path)
cat_gtf "$GENE_GTF" |
"$AWK" 'BEGIN{OFS="\t"} $3=="exon"{print $1,$4-1,$5,".",".",$7}' > "$TMP/exons.bed"

# 3b) TE BED from the grouped SAF (0-based), GeneID carried as the BED name
"$AWK" 'BEGIN{OFS="\t"} NR>1{print $2,$3-1,$4,$1,".",$5}' "$SAF_ALL" > "$TMP/te_grouped.bed"

# 3c) normalize contig names on BOTH BEDs (default: strip leading 'chr')
if [ "$STRIP_CHR" = "yes" ]; then
  sed -E 's/^chr//' "$TMP/exons.bed"      > "$TMP/exons.norm"      && mv "$TMP/exons.norm" "$TMP/exons.bed"
  sed -E 's/^chr//' "$TMP/te_grouped.bed" > "$TMP/te_grouped.norm" && mv "$TMP/te_grouped.norm" "$TMP/te_grouped.bed"
fi

# 3d) sort both (LC_ALL=C, -k1,1 -k2,2n)
LC_ALL=C sort -k1,1 -k2,2n "$TMP/exons.bed"      -o "$TMP/exons.bed"
LC_ALL=C sort -k1,1 -k2,2n "$TMP/te_grouped.bed" -o "$TMP/te_grouped.bed"

# 3e) subtract, then VERIFY 0 residual overlap
n_overlap=$(bedtools intersect -a "$TMP/te_grouped.bed" -b "$TMP/exons.bed" -u | wc -l)
echo "[build_te_saf]   TE loci overlapping exons (pre-subtract): $n_overlap"
bedtools subtract -a "$TMP/te_grouped.bed" -b "$TMP/exons.bed" > "$TMP/te_grouped_noExon.bed"
n_residual=$(bedtools intersect -a "$TMP/te_grouped_noExon.bed" -b "$TMP/exons.bed" -u | wc -l)
if [ "$n_residual" -ne 0 ]; then
  echo "ERROR: $n_residual residual exon overlaps after subtract (expected 0)." >&2
  exit 1
fi

# 3f) rebuild SAF (Start back to 1-based)
"$AWK" 'BEGIN{OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}
     {print $4,$1,$2+1,$3,$6}' "$TMP/te_grouped_noExon.bed" > "$SAF_NOEXON"

ne_groups=$(cut -f1 "$SAF_NOEXON" | tail -n +2 | sort -u | wc -l)
echo "[build_te_saf]   no-exon SAF: $(( $(wc -l < "$SAF_NOEXON") - 1 )) loci, $ne_groups unique groups (0 residual exon overlap)"
echo "[build_te_saf] DONE."
echo "  $SAF_ALL"
echo "  $SAF_NOEXON"
