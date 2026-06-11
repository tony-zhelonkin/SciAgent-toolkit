---
name: te-gene-featurecounts
description: >-
  featureCounts TE+gene counting workflow, packaged with its own locked,
  version-pinned container (te-fc:2.0.2, featureCounts v2.0.2). Two-pass driver:
  integer Random-One TE counting (grouped Subfamily:Family:Class SAF, -M, -s 0,
  NO --fraction) plus per-library-stranded gene counting, then a row-bound
  combined matrix. Use when you have nf-core/rnaseq star_salmon BAMs (including
  the lean markdup.sorted.bam path) and need the runnable, env-locked step that
  turns them into gene + TE subfamily count matrices. For the upstream STAR
  Random-One alignment recipe use star-te-preprocessing; for building the TE SAF
  use te-reference-saf-build; for annotating the resulting matrices use
  annotate-bulk-rnaseq-data.
license: MIT
metadata:
  skill-author: SciAgent-toolkit
  last-reviewed: 2026-06-11
  version: 1.0.0
  upstream-docs: https://subread.sourceforge.net/
  scope: implementation
  category: workflow
  tier: packaged
  tags:
    - preprocessing
  requires: []
  complementary-skills:
    - te-reference-saf-build
    - star-te-preprocessing
    - annotate-bulk-rnaseq-data
  contraindications:
    - "Do not use for fractional multimapper counting. This is integer Random-One (-M, -s 0, NO --fraction); adding --fraction breaks the contract."
    - "Do not use for locus-level / copy-resolved TE quantification. The grouped SAF is subfamily-level; use SQuIRE/Telescope instead."
    - "Do not use to build the TE SAF or run STAR. The SAF is built by te-reference-saf-build and the Random-One BAMs by star-te-preprocessing; this skill begins at pre-built BAMs + SAF."
    - "Do not run featureCounts from scdock-r-dev:v0.5.x — those images lack subread. Use the locked te-fc:2.0.2 (or legacy scdock-r-dev:v0.2)."
    - "Do not perform DE, annotation, or DGEList assembly here. Hand the matrices to annotate-bulk-rnaseq-data."
---

# TE + Gene featureCounts Counting (env-locked, packaged)

## Overview

This skill is the **self-contained, runnable** step that turns nf-core/rnaseq STAR
`star_salmon` BAMs into **gene + TE subfamily count matrices**, crystallized with its own
minimal, version-pinned Docker container (`te-fc:2.0.2`, featureCounts v2.0.2). It vendors
the two-pass featureCounts driver (frozen copies) and a parameterized Docker wrapper that
handles the proven symlink-staging + identical-path bind-mount pattern. Unlike
`star-te-preprocessing` (which owns the alignment + counting *contract* and the canonical
STAR string), this skill owns the **executable counting artifact**: the locked image, the
driver, the wrapper, and the QC gate the matrices must pass.

**When to use this skill:**
- You have nf-core/rnaseq `star_salmon` BAMs (full or lean `markdup.sorted.bam`) and a
  pre-built grouped TE SAF, and need gene + TE subfamily count matrices.
- You want a reproducible, env-locked counting run (pinned featureCounts v2.0.2).
- You need the exact Docker staging recipe so `bam_fin` symlinks resolve and column names
  come out clean.

**When NOT to use this skill:**
- Building the TE SAF → `te-reference-saf-build`.
- Running STAR / producing the Random-One BAMs → `star-te-preprocessing`.
- Annotating / DGEList / DE on the matrices → `annotate-bulk-rnaseq-data`.
- Locus-level TE quantification → SQuIRE / Telescope (out of scope).

---

## Decision Tree

```
Have star_salmon BAMs + a grouped TE SAF, need count matrices?
│
├─ Want integer subfamily-level TE + gene counts, reproducibly?
│     →  THIS skill: run scripts/run_te_counting.sh in te-fc:2.0.2
│
├─ Need to build the TE SAF first?            →  te-reference-saf-build
├─ Need to produce the Random-One BAMs first? →  star-te-preprocessing
├─ Want fractional 1/n multimapper counts?    →  NOT current practice (do not add --fraction)
├─ Want locus-level / copy-resolved TE?       →  SQuIRE / Telescope (out of scope)
└─ Already have matrices, want annotation/DE? →  annotate-bulk-rnaseq-data
```

---

## The locked-env contract

- **Image:** `te-fc:2.0.2` — minimal Debian-slim + **only** featureCounts **v2.0.2**
  (subread). No R, no Python. Built from `env/Dockerfile` via `env/build.sh`.
- **Why v2.0.2 from the official binary:** v2.0.2 is the exact version 13036-DM and 14839-DM
  ran (verified from raw headers `# Program:featureCounts v2.0.2`). bioconda **skipped**
  packaging subread 2.0.2 (it jumps 2.0.1 → 2.0.3), so the pin is satisfied by installing the
  official subread-2.0.2 Linux release binary (the same standalone binary the precedent image
  `scdock-r-dev:v0.2` shipped). Build fails closed if the banner is not `v2.0.2`.
- **Image lineage:** legacy `scdock-r-dev:v0.2` also has v2.0.2 (original runs). The newer
  `scdock-r-dev:v0.5.x` images do **NOT** contain featureCounts (only MultiQC's parser). Use
  the locked `te-fc:2.0.2` for all new runs.

```bash
# Build the locked image (idempotent; aborts if < 5G free on /):
bash env/build.sh
docker run --rm te-fc:2.0.2 featureCounts -v   # -> featureCounts v2.0.2
```

---

## Quick Start

```bash
# One command does staging + identical-path mount + both passes + combined matrix,
# all inside te-fc:2.0.2:
scripts/run_te_counting.sh \
  --bam-dir   /data2/nf-results/<proj>/results_mm39_TE/star_salmon \
  --bam-glob  '*.markdup.sorted.bam' \          # lean nf-core run; default '*.bam'
  --gene-gtf  <outdir>/genome/gencode.vM37.primary_assembly.annotation.filtered.gtf \
  --te-saf    /data1/shared/ref/mouse/Ensembl/mm39/GRCm39_rmsk_TE_GROUPED_all_noExon.saf \
  --gene-s    2 \                               # GENE strandedness — VERIFY per library
  --out-dir   <OUT_DIR> \
  --threads   12
```

**Verify it worked (the QC gate — full checklist below):**

```bash
M=<OUT_DIR>/featurecounts_TE/te_counts_matrix.txt
awk 'NR>1{for(i=2;i<=NF;i++) if($i!=int($i)){print "FRACTIONAL!"; exit 1}}' "$M"  # integer
head -1 "$M"; wc -l "$M"   # ~1,243 TE meta-features (mm39), Subfamily:Family:Class rows
```

---

## Progressive Depth

### Basic Usage — the wrapper

`scripts/run_te_counting.sh` is the entry point. It:
1. Symlinks the matched BAMs into `<OUT_DIR>/bam_fin/` (no copies).
2. Bind-mounts the symlink-**target** dir(s) at an **identical host=container path**
   (`-v $REAL:$REAL:ro`) so links resolve inside the container and the driver's awk parser
   keys clean `<sample>` column names off the BAM basename.
3. Runs the vendored `runFeatureCounts_TE_and_genes.sh` (TE pass → gene pass → row-bind) in
   `te-fc:2.0.2` as `-u $(id -u):$(id -g)`.

Required args: `--bam-dir --gene-gtf --te-saf --gene-s --out-dir`. Optional: `--threads`
(12), `--te-strand` (`unstranded`), `--bam-glob` (`*.bam`), `--image` (`te-fc:2.0.2`).

### Intermediate Usage — strandedness (the error-prone variable)

The **gene** `-s` is **library-specific and must be verified per dataset** — never hardcode.
Confirm against MultiQC inferred strandedness, RSeQC/Salmon, AND the featureCounts header
(`Strand specific : reversely stranded`). 14839-DM and 13036-DM were `-s 2` (reverse, dUTP/
TruSeq); AdaW was `-s 1` (forward). The **TE** pass is **always `-s 0` (unstranded)** to avoid
halving signal from antisense TE transcription — fixed inside the driver via
`--te-strand unstranded`.

### Advanced Usage — the vendored driver (the code is the spec)

`scripts/runFeatureCounts_TE_and_genes.sh` (gene pass delegates to `scripts/runFeatureCounts.sh`)
are **frozen, vendored copies** of the TE-RNAseq-toolkit drivers — self-contained so the skill
is a runnable artifact (a deliberate reversal of the prior version-pointer ADR). Comments were
corrected vs the originals (grouped no-exon SAF; integer Random-One; `-s 0`, no `--fraction`);
**code logic is byte-identical**. Underlying invocations:

```
TE:   featureCounts -M -F SAF -a <SAF> -o te_counts_raw.txt -s 0 -p --countReadPairs -B -C -T <t> <BAMs>
Gene: featureCounts -a <GTF> -o counts_matrix.txt -p --countReadPairs -B -C -s <0|1|2> -t exon -g gene_id -T <t> <BAMs>
```

- **TE pass:** `-M` (multi-mappers counted — essential for TEs), `-s 0`, **NO `--fraction`**
  → integer Random-One. SAF `GeneID = Subfamily:Family:Class` → ~1,243 subfamily meta-features.
- **Gene pass:** multi-mappers excluded (featureCounts default), `-s` per library.
- **Combine:** row-binds gene + TE into `combined_gene_TE_counts.tsv` (valid only because
  exonic TE loci were subtracted from the SAF → no double-counting).

The full end-to-end runbook (inputs, lean BAM path, staging recipe, QC gate, handoff) lives in
**`references/te-counting-workflow.md`**.

---

## Outputs

| File | Path |
|---|---|
| TE subfamily matrix | `<OUT>/featurecounts_TE/te_counts_matrix.txt` |
| Gene matrix | `<OUT>/fc_genes/count_matrices_fc/sorted_counts_matrix.txt` |
| Combined (row-bind) | `<OUT>/combined_gene_TE_counts.tsv` |

14839-DM dims: TE 1,243 × 45; gene 78,317 × 45; combined 79,560 × 45.

---

## Verification Checklist (the QC gate)

After running, confirm before handoff:

- [ ] **Integer TE counts** — no fractional values (Random-One, no `--fraction`).
- [ ] **TE label shape** — every TE row is 3-field `Subfamily:Family:Class` (exactly 2 colons).
- [ ] **Sample order == samplesheet** — matrix columns match samplesheet order 1:1; gene header
      == TE header.
- [ ] **No zero-libsize samples** — every per-sample gene and TE total is nonzero.
- [ ] **TE proportion = a library-specific sanity band, NOT a hard threshold.** Expect internal
      consistency across replicates; flag *wild* outliers, not an absolute number. The
      unstranded-TE / stranded-gene ratio **inflates** TE% (the gene denominator drops
      antisense/ambiguous reads the unstranded TE pass keeps) — 14839-DM saw 5.5–12.5%
      (mean 8.7%), internally consistent vs the AdaW ~3.8–6.0% reference (different tissue).
- [ ] **featureCounts version** — raw headers read `# Program:featureCounts v2.0.2`.

---

## Common Pitfalls

| Symptom | Cause | Fix |
|---|---|---|
| `featureCounts: not found` in container | Used `scdock-r-dev:v0.5.x` (no subread) | Use `te-fc:2.0.2` (or legacy `scdock-r-dev:v0.2`). |
| BAM symlinks fail to open inside container | Target dir not mounted, or mounted at a different path | Mount the symlink-target dir at an **identical host=container path** (`run_te_counting.sh` does this). |
| Column names are full paths, not sample names | BAMs passed by a path the awk parser can't reduce to a basename | Use the identical-path mount + `bam_fin/` symlinks so the basename parser yields `<sample>`. |
| Fractional values in TE matrix | `--fraction` was added | Remove it; this recipe is integer `-M` Random-One. |
| Gene counts ~half expected / near zero | Wrong gene `-s` (forward vs reverse) | Set `--gene-s` from MultiQC + featureCounts header, per library. Never assume 1 vs 2. |
| Combined matrix double-counts a region | SAF still contains exonic TE loci | Use the `*_noExon.saf` from `te-reference-saf-build`. |

---

## Complementary Skills

| When you need... | Use skill | Relationship |
|---|---|---|
| Build the grouped, exon-subtracted TE SAF | `te-reference-saf-build` | Prerequisite (produces the SAF this skill consumes) |
| Produce the Random-One STAR BAMs + the alignment/counting contract | `star-te-preprocessing` | Prerequisite (produces the BAMs; owns the contract) |
| Annotate matrices, parse TE IDs, build combined DGEList | `annotate-bulk-rnaseq-data` | Next step (consumes these matrices) |
| Locus-level / copy-resolved TE quantification | SQuIRE / Telescope (external) | Alternative (out of scope) |

The canonical chain is `te-reference-saf-build` + `star-te-preprocessing` → **`te-gene-featurecounts`** → `annotate-bulk-rnaseq-data`.

---

## Resources

- **Locked image:** `env/Dockerfile` + `env/build.sh` → `te-fc:2.0.2` (featureCounts v2.0.2).
- **Vendored drivers:** `scripts/runFeatureCounts_TE_and_genes.sh`, `scripts/runFeatureCounts.sh` (frozen, comments corrected).
- **Wrapper:** `scripts/run_te_counting.sh` (staging + identical-path mount + run).
- **Runbook:** `references/te-counting-workflow.md` (end-to-end, QC gate, handoff).
- **Smoke tests:** `tests/run_skill_tests.sh` (image present, v2.0.2, synthetic integer matrix).
- **subread/featureCounts:** https://subread.sourceforge.net/ (release 2.0.2).
