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
                                       # --te-strand unstranded (default); TE pass fixed -s 0
```

The two passes (contract enforced by the vendored driver):

- **Gene pass:** GTF, `-t exon -g gene_id`, strandedness **`-S` per verified library** (14839/13036
  = `-s 2` reverse; AdaW = `-s 1` forward), multi-mappers **excluded** (featureCounts default),
  `-p --countReadPairs -B -C`. **Never assume 1 vs 2** — confirm against MultiQC inferred
  strandedness, RSeQC/Salmon, and the featureCounts header (`Strand specific : ...`).
- **TE pass:** SAF, `-M` (multi-mappers **counted** — essential for TEs), `-s 0` (**unstranded**),
  **NO `--fraction`** → **integer** Random-One counts, `-p --countReadPairs -B -C`.

### Exact underlying invocations (from raw headers, 14839-DM)

```
TE:   featureCounts -M -F SAF -a <SAF> -o te_counts_raw.txt -s 0 -p --countReadPairs -B -C -T 12 <BAMs>
Gene: featureCounts -a <GTF> -o counts_matrix.txt -p --countReadPairs -B -C -s 2 -t exon -g gene_id -T 12 <BAMs>
```

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
  library/tissue, so an offset is expected. Crucially, the **unstranded TE pass vs reverse-stranded
  gene denominator** inflates TE% (the gene denominator drops antisense/ambiguous reads the
  unstranded TE pass keeps): 14839-DM measured 5.5–12.5% (mean 8.7%), internally consistent — a
  QC *observation*, not an error.

**Handoff:** the gene + TE matrices are the inputs to **`annotate-bulk-rnaseq-data`** (Ensembl→Symbol
gene annotation, `parse_te_id` TE parsing, combined annotated `DGEList`) → then DE/GSEA. Do not
perform DE here.

---

## Reference paths cited (14839-DM run)

| Asset | Path |
|---|---|
| BAM inputs (lean) | `/data2/nf-results/14839-DM/results_mm39_TE/star_salmon/*.markdup.sorted.bam` |
| Grouped no-exon SAF | `/data1/shared/ref/mouse/Ensembl/mm39/GRCm39_rmsk_TE_GROUPED_all_noExon.saf` |
| Gene GTF (nf-core filtered) | `<outdir>/genome/gencode.vM37.primary_assembly.annotation.filtered.gtf` |
| Vendored driver | `scripts/runFeatureCounts_TE_and_genes.sh` (gene pass → `scripts/runFeatureCounts.sh`) |
| Container | `te-fc:2.0.2` (featureCounts v2.0.2) |
