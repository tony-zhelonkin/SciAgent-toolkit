# TE + gene featureCounts counting workflow (post-nf-core, end-to-end runbook)

The operational runbook for going **from nf-core/rnaseq STAR BAMs → gene + TE subfamily
count matrices**, run inside the **locked `te-fc:2.0.2` container** (featureCounts/subread
2.0.2). SKILL.md owns the env contract and the QC gate; this file is the step-by-step
procedure. Provenance: generalized from the proven 14839-DM / 13036-DM runs
(`/scratch/nf-core/14839-DM/run_te_counting.sh` + `te_counting_record.md`). The
upstream STAR "Random-One" alignment recipe + SAF *construction* are owned by
`star-te-preprocessing` and `te-reference-saf-build`; this skill begins at the BAMs and
a pre-built SAF.

Reference paths below are the exact mm39 assets from the 14839-DM run; substitute the
project's own outdir/paths.

---

## 0. The locked container

```
IMAGE = te-fc:2.0.2        # minimal: featureCounts v2.0.2 only (no R, no Python)
```

Build once: `bash env/build.sh` (or `docker build -t te-fc:2.0.2 env/`). The image installs
the **official subread-2.0.2 Linux binary** (bioconda skipped packaging 2.0.2; the pin can
only be satisfied from the SourceForge release — same standalone binary the precedent image
`scdock-r-dev:v0.2` shipped). featureCounts v2.0.2 is the exact version 13036-DM and 14839-DM
ran (verified from raw output headers `# Program:featureCounts v2.0.2`).

> **Image lineage note:** the *legacy* `scdock-r-dev:v0.2` also contains featureCounts v2.0.2
> and was used for the original runs. The newer `scdock-r-dev:v0.5.x` images do **NOT** contain
> featureCounts/subread (only MultiQC's parser). For all new counting runs use the locked,
> minimal `te-fc:2.0.2`.

## 1. Inputs from nf-core (the STAR BAMs)

featureCounts consumes the **STAR genome BAMs** published under `--aligner star_salmon`, NOT the
Salmon quant. Depending on the nf-core run:

- **Full run** (`save_align_intermeds = true`): `star_salmon/<sample>/Aligned.out.bam` (unsorted)
  and `Aligned.sortedByCoord.out.bam` are available.
- **Lean run** (`save_align_intermeds = false`, e.g. 14839-DM): only the published
  `star_salmon/<sample>.markdup.sorted.bam` survives. **This is valid TE input** — `markdup`
  *marks* (does not remove) duplicates, and featureCounts counts them by default. Use the
  `.markdup.sorted.bam` path.

Either sorted or unsorted is acceptable — featureCounts pairs reads internally
(`METHODOLOGY.md:52-53`). Do not delete published BAMs from the outdir.

## 2. The TE SAF (pre-built, grouped, exon-subtracted)

The TE pass needs the **grouped, exon-subtracted SAF**, built once per genome build by
`te-reference-saf-build`:

```
TE_SAF = /data1/shared/ref/mouse/Ensembl/mm39/GRCm39_rmsk_TE_GROUPED_all_noExon.saf
```

The SAF `GeneID` is the group label `Subfamily:Family:Class` (e.g. `L1Md_A:L1:LINE`);
featureCounts sums all loci sharing a label into one subfamily meta-feature (~1,243 groups
for mm39). Exonic TE loci were subtracted (`bedtools subtract`) so a read counts as gene OR
TE, never both — the precondition for a valid combined matrix. (NOT the stale `_formatted`
SAF; NOT a locus-level SAF.)

## 3. BAM staging — the critical Docker pattern

featureCounts keys output column names off the BAM **basename**, and the bind mounts must let
symlinks resolve inside the container. The proven pattern (`run_te_counting.sh` automates it):

1. **Symlink** (not copy) the project BAMs into a `bam_fin/` dir.
2. Bind-mount the symlink-**target** dir(s) at an **IDENTICAL host=container path**
   (`-v $REAL:$REAL:ro`) so the links resolve to their real targets inside the container.
3. Mount `bam_fin/`, the output dir, the GTF dir, the SAF dir, and the scripts dir likewise.
4. Run as `-u $(id -u):$(id -g)` so outputs are owned by you.

With identical host=container paths the `bam_fin` symlinks resolve to their real targets and
the driver's awk parser yields clean `<sample>` column names from the basename.

## 4. Run the counting (two passes via the wrapper)

The packaged wrapper `scripts/run_te_counting.sh` does the staging + identical-path mount +
runs the vendored orchestrator `runFeatureCounts_TE_and_genes.sh` (TE pass → gene pass →
row-bind combined) inside `te-fc:2.0.2`:

```bash
scripts/run_te_counting.sh \
  --bam-dir   <STAR_BAM_DIR> \         # dir of *.bam (e.g. star_salmon)
  --bam-glob  '*.markdup.sorted.bam' \ # match the lean-run BAMs (default '*.bam')
  --gene-gtf  <GENE_GTF> \             # nf-core filtered GTF
  --te-saf    $TE_SAF \                # grouped exon-subtracted SAF (step 2)
  --gene-s    2 \                      # GENE strandedness — VERIFY per library
  --out-dir   <OUT_DIR> \
  --threads   12
                                       # --te-strand unstranded (default, -s 0) OR sense_antisense
```

The two passes (contract enforced by the vendored driver):

- **Gene pass:** GTF, `-t exon -g gene_id`, strandedness **`-S` per verified library** (14839/13036
  = `-s 2` reverse; AdaW = `-s 1` forward), multi-mappers **excluded** (featureCounts default),
  `-p --countReadPairs -B -C`. **Never assume 1 vs 2** — confirm against MultiQC inferred
  strandedness, RSeQC/Salmon, and the featureCounts header (`Strand specific : ...`).
- **TE pass:** SAF, `-M` (multi-mappers **counted** — REQUIRED under Random-One; without `-M`
  featureCounts drops the `NH>1` reads STAR keeps on the single emitted line), **NO `--fraction`**
  → **integer** Random-One counts, `-p --countReadPairs -B -C`. TE strandedness is
  **context-dependent and the field is SPLIT, NOT a fixed `-s 0`** (see "Why" below):
  - **Standalone TE-family quantification or a non-directional library →** `-s 0` (unstranded),
    the wrapper default. **Defensible (grade B): matches the dominant tool's default**
    (TEtranscripts `--stranded no`). Counts TE reads on either strand; a specificity *trade*, not
    a sensitivity gain.
  - **JOINT gene+TE matrix on a STRANDED library (more principled best-practice, grade B) →**
    count TEs at the **same strandedness as genes** and keep bidirectional biology via a
    **sense/antisense split** rather than collapsing to `-s 0`. `--te-strand sense_antisense`
    emits TE-sense (`-s 2` for a reverse lib, matched to genes) + TE-antisense (`-s 1`) matrices.
    Better FDR in the one benchmark and separates autonomous from passive TE transcription, but
    **not proven superior** for TE DE; mode-switching is **not required**.

> **Why.** Two over-claims to avoid in both directions. (1) "Always `-s 0` for bidirectional TEs"
> is mis-attributed to Teissandier 2019 (benchmarks *multimapper handling only* — "strand" appears
> once, as a fixed `-s 0` parameter; no strand-choice test). (2) "Stranded is THE field standard"
> is *also* over-claimed: the dominant tool TEtranscripts/TEcount **defaults to `--stranded no`
> (unstranded)**; best-practice (TE-Seq 2025) *recommends* stranded for directional libraries. So
> stranded-for-joint is **best-practice / mechanistic (grade B), not a benchmarked standard**, and
> the field genuinely splits. The "unstranded → better TE sensitivity" claim is **GAP** (never
> benchmarked); the only both-mode study (Savytska 2022, doi:10.3389/fgene.2022.1026847) found
> **stranded FDR (54.9%) ≤ unstranded (58.7%)**. Real TE bidirectionality is class-specific (L1-ASP/ORF0, LTR/ERV
> real; SINE/intronic largely passive), not a uniform property of TE loci.

### Exact underlying invocations

```
TE (unstranded, default):     featureCounts -M -F SAF -a <SAF> -o te_counts_raw.txt -s 0 -p --countReadPairs -B -C -T 12 <BAMs>
TE (sense, reverse lib):      featureCounts -M -F SAF -a <SAF> -o te_counts_sense_raw.txt -s 2 -p --countReadPairs -B -C -T 12 <BAMs>   # INTEGER Random-One (no --fraction)
TE (antisense, reverse lib):  featureCounts -M -F SAF -a <SAF> -o te_counts_antisense_raw.txt -s 1 -p --countReadPairs -B -C -T 12 <BAMs>   # INTEGER Random-One (no --fraction)
Gene:                         featureCounts -a <GTF> -o counts_matrix.txt -p --countReadPairs -B -C -s 2 -t exon -g gene_id -T 12 <BAMs>
```

All TE passes — primary unstranded AND the optional sense/antisense auxiliaries — are integer
Random-One (`-M`, no `--fraction`): one kernel everywhere, integer because STAR Random-One emits one
alignment/read. The fractional route (`-M --fraction` → non-integer → round()/limma-voom before
DESeq2) is a labeled **non-default alternative** (Strategy B; requires STAR `--outSAMmultNmax 100`).

(The TE `-s 0` line above is the exact 14839-DM run — a defensible standalone choice that matches
TEtranscripts' default, not retroactively wrong; for the definitive joint analysis the
sense/antisense lines matched to gene `-s 2` are the more principled best-practice, grade B.)

## 5. Outputs, QC gate, handoff

Outputs land under `<OUT_DIR>`:

| File | Path |
|---|---|
| TE subfamily matrix | `featurecounts_TE/te_counts_matrix.txt` |
| Gene matrix | `fc_genes/count_matrices_fc/sorted_counts_matrix.txt` |
| Combined (row-bind) | `combined_gene_TE_counts.tsv` |
| TE/gene raw + `.summary` | `featurecounts_TE/te_counts_raw.txt`, `fc_genes/raw_fc_output/counts_matrix.txt` |

14839-DM dims: TE 1,243 × 45; gene 78,317 × 45; combined 79,560 × 45.

**QC gate** before handing off:

- **Integer TE counts** — no fractional values (`awk 'NR>1{for(i=2;i<=NF;i++) if($i!=int($i)) exit 1}'`).
- **TE label shape** — every TE row is 3-field `Subfamily:Family:Class` (exactly 2 colons).
- **Sample order == samplesheet** — matrix columns match the samplesheet order 1:1; gene header
  == TE header.
- **No zero-libsize samples** — every per-sample gene and TE total is nonzero.
- **TE proportion = a LIBRARY-SPECIFIC sanity band, NOT a hard threshold.** TE% =
  TE_total / (gene_total + TE_total) per sample. Expect internal consistency across replicates;
  flag *wild* outliers, not an absolute number. The AdaW reference (~3.8–6.0%) is a *different*
  library/tissue, so an offset is expected. When TEs are counted **`-s 0` (unstranded) against a
  reverse-stranded gene denominator**, the mismatch inflates TE% (the gene denominator drops
  antisense/ambiguous reads the unstranded TE pass keeps): 14839-DM measured 5.5–12.5% (mean 8.7%)
  under standalone `-s 0`, internally consistent — a QC *observation*, not an error. So "TE %" is a
  QC sanity band, **not a biological transcriptome fraction**; a stranded TE recount removes the
  mismatch.

**Handoff:** the gene + TE matrices are the inputs to **`annotate-bulk-rnaseq-data`** (Ensembl→Symbol
gene annotation, `parse_te_id` TE parsing, combined annotated `DGEList`) → then DE/GSEA. Do not
perform DE here. **Joint-analysis caveats to carry with the matrix (graded options, not mandates;
see SKILL.md "Evidence & open questions"):** genes-only DESeq2 size factors (`controlGenes`) are
**grade B / contested** (TEtranscripts pools genes+TEs); the sense/antisense split is **grade B /
SQuIRE-specific**; the joint matrix is valid for **within-feature-type, across-sample DE only** —
never compare gene-vs-TE magnitude within a sample, and emit no TPM/FPKM for TE meta-features
(**grade C / mechanistic inference**, not stated in any TE primary source).

---

## Reference paths cited (14839-DM run)

| Asset | Path |
|---|---|
| BAM inputs (lean) | `/data2/nf-results/14839-DM/results_mm39_TE/star_salmon/*.markdup.sorted.bam` |
| Grouped no-exon SAF | `/data1/shared/ref/mouse/Ensembl/mm39/GRCm39_rmsk_TE_GROUPED_all_noExon.saf` |
| Gene GTF (nf-core filtered) | `<outdir>/genome/gencode.vM37.primary_assembly.annotation.filtered.gtf` |
| Vendored driver | `scripts/runFeatureCounts_TE_and_genes.sh` (gene pass → `scripts/runFeatureCounts.sh`) |
| Container | `te-fc:2.0.2` (featureCounts v2.0.2) |
