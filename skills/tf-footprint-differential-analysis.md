# TF Footprint & Differential Footprint Analysis (TOBIAS + HINT-ATAC + Signac)

## Purpose

Identify transcription factor binding sites from ATAC-seq footprints and compare TF occupancy across biological conditions. Footprinting detects the physical protection of DNA by bound TFs within accessible chromatin, complementing motif accessibility (ChromVAR) with occupancy evidence.

**When to use:**
- Differential TF binding between conditions (e.g., WT vs KO)
- Validating ChromVAR TF activity with occupancy evidence
- Identifying condition-specific TF regulators from pseudobulk scATAC-seq
- Comparing TF binding across cell types

**Upstream:** Peak calling (MACS2/MACS3), pseudobulk BAM generation (ArchR/Signac), ChromVAR deviations
**Downstream:** TF regulatory networks, ChromVAR integration (quadrant analysis), GRN inference (SCENIC+)

**Key distinction:**
- **ChromVAR** measures **accessibility** (potential for binding — chromatin is open at motif sites)
- **Footprinting** measures **occupancy** (actual binding — protein physically protects DNA from Tn5)

---

## Tool Selection Guide

| Feature | **TOBIAS** (Primary) | **HINT-ATAC** (Alternative) | **Signac Footprint()** (Visualization) |
|:---|:---|:---|:---|
| **Algorithm** | Tn5 bias correction (linear regression) + footprint scoring | HMM-based detection (learns footprint states) | Normalized Tn5 insertion profiles |
| **Tn5 Bias Correction** | Best in class (`ATACorrect`) | Good (bias table) | Basic (k-mer model) |
| **Differential Analysis** | Excellent (`BINDetect` — differential binding scores + p-values) | Good (`rgt-hint differential`) | Weak (visualization only) |
| **Pseudobulk Support** | High (designed for merged BAMs) | High (designed for merged BAMs) | Medium (fragment files, slow on large data) |
| **Output** | BigWigs, volcano plots, TSV tables | BED files, line plots, statistics TXT | Seurat objects, ggplot2 |
| **Best For** | Quantitative differential analysis | Multi-condition heatmaps | Publication figures |

**Recommendation:** Use **TOBIAS** for quantitative analysis, **Signac** for visualization, **HINT-ATAC** optionally for multi-condition overview heatmaps.

---

## Installation

### TOBIAS (Python — primary)

```bash
# Create venv and install
python3.10 -m venv /opt/venvs/tobias
source /opt/venvs/tobias/bin/activate
pip install tobias

# Verify
TOBIAS --version
```

### HINT-ATAC / RGT (Python — alternative)

```bash
# Create venv and install
python3.9 -m venv /opt/venvs/rgt
source /opt/venvs/rgt/bin/activate
pip install rgt

# Download genome data (required for bias correction)
python -c "from rgt.helper import download_genome; download_genome('mm10')"

# Verify
rgt-hint footprinting --help
```

### Signac Footprint (R — visualization)

```r
# Already available in Signac
library(Signac)
library(motifmatchr)
library(JASPAR2020)        # or JASPAR2024
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
```

---

## TOBIAS Workflow (Primary — Pseudobulk)

### Overview

```
BAM + Peaks → ATACorrect → ScoreBigwig → BINDetect → Differential Binding
```

### Step 1: Tn5 Bias Correction

```bash
# Run per pseudobulk BAM (one per celltype×condition)
TOBIAS ATACorrect \
  --bam ${CELLTYPE}_${CONDITION}.bam \
  --genome mm10.fa \
  --peaks consensus_peaks.bed \
  --blacklist mm10-blacklist.v2.bed \
  --outdir ./corrected/ \
  --cores 16

# Outputs:
# - *_corrected.bw    (bias-corrected signal — USE THIS)
# - *_uncorrected.bw  (raw signal)
# - *_expected.bw     (expected Tn5 bias)
# - *_bias.bw         (estimated bias)
```

**Critical parameters:**
- `--blacklist`: ENCODE mm10 blacklist v2 (ENCFF547MET). **Mandatory** — footprinting finds "perfect" footprints in artifact regions.
- `--genome`: Must match alignment build exactly. Chromosome naming (`chr1` vs `1`) is the #1 cause of failures.
- Do NOT pre-shift reads — TOBIAS applies the +4/-5 bp Tn5 shift internally.

### Step 2: Footprint Scoring

```bash
# Generate footprint scores from corrected signal
TOBIAS ScoreBigwig \
  --signal ./corrected/${CELLTYPE}_${CONDITION}_corrected.bw \
  --regions consensus_peaks.bed \
  --output ./scores/${CELLTYPE}_${CONDITION}_footprints.bw \
  --cores 16
```

### Step 3: Motif Scanning (TFBScan)

```bash
# Scan peaks for motif occurrences
TOBIAS TFBScan \
  --motifs motifs.jaspar \
  --fasta mm10.fa \
  --regions consensus_peaks.bed \
  --outdir ./motif_hits/ \
  --cores 16
```

**Motif format:** JASPAR `.jaspar` or MEME `.meme` format. See motif database section below.

### Step 4: Differential Binding (BINDetect)

```bash
# Compare two conditions for one celltype
TOBIAS BINDetect \
  --motifs motifs.jaspar \
  --signals ./scores/cDC1_WT_footprints.bw ./scores/cDC1_Batf3KO_footprints.bw \
  --genome mm10.fa \
  --peaks consensus_peaks.bed \
  --outdir ./differential/cDC1_WT_vs_KO/ \
  --cond_names WT Batf3_KO \
  --cores 16

# Outputs:
# - bindetect_results.txt   (per-TF differential binding scores + p-values)
# - bindetect_results.xlsx  (Excel format)
# - <TF_name>/              (per-TF BED files of bound/unbound sites)
# - volcano plot PDF
```

**Key output columns in `bindetect_results.txt`:**
- `<cond1>_mean_score` / `<cond2>_mean_score`: Average footprint score per condition
- `<cond1>_<cond2>_change`: Differential binding score (positive = more bound in cond1)
- `<cond1>_<cond2>_pvalue`: Statistical significance

### Step 5: Visualization (PlotAggregate)

```bash
# Plot aggregate footprint profiles for specific TFs
TOBIAS PlotAggregate \
  --TFBS ./differential/cDC1_WT_vs_KO/IRF8/beds/IRF8_cDC1_WT_bound.bed \
  --signals ./corrected/cDC1_WT_corrected.bw ./corrected/cDC1_Batf3KO_corrected.bw \
  --output IRF8_aggregate.pdf \
  --share_y both \
  --plot_boundaries \
  --signal-on-x
```

### Batch Processing (All Contrasts)

```bash
#!/bin/bash
# Run TOBIAS for all standard contrasts

GENOME="mm10.fa"
PEAKS="consensus_peaks.bed"
BLACKLIST="mm10-blacklist.v2.bed"
MOTIFS="motifs.jaspar"
CORES=16

# Define contrasts as arrays
declare -A CONTRASTS=(
  ["I3"]="cDC1A_WT:cDC1B_WT"
  ["P1"]="cDC1_WT:cDC1_WTplusIL12"
  ["P4"]="cDC1_WT:cDC1_Batf3KO"
  ["I1"]="DC_WT:Mac_WT"
)

# Step 1-2: Correct and score all BAMs
for bam in ./bams/*.bam; do
  prefix=$(basename "$bam" .bam)
  TOBIAS ATACorrect --bam "$bam" --genome "$GENOME" --peaks "$PEAKS" \
    --blacklist "$BLACKLIST" --outdir ./corrected/ --cores "$CORES"
  TOBIAS ScoreBigwig --signal "./corrected/${prefix}_corrected.bw" \
    --regions "$PEAKS" --output "./scores/${prefix}_footprints.bw" --cores "$CORES"
done

# Step 3: Differential binding for each contrast
for contrast_id in "${!CONTRASTS[@]}"; do
  IFS=':' read -r cond1 cond2 <<< "${CONTRASTS[$contrast_id]}"
  TOBIAS BINDetect --motifs "$MOTIFS" \
    --signals "./scores/${cond1}_footprints.bw" "./scores/${cond2}_footprints.bw" \
    --genome "$GENOME" --peaks "$PEAKS" \
    --outdir "./differential/${contrast_id}_${cond1}_vs_${cond2}/" \
    --cond_names "$cond1" "$cond2" --cores "$CORES"
done
```

---

## HINT-ATAC Workflow (Alternative — Multi-Condition)

### Footprint Calling

```bash
# Per celltype pseudobulk BAM
rgt-hint footprinting \
  --atac-seq \
  --paired-end \
  --organism=mm10 \
  --output-location=./footprints/ \
  --output-prefix=${CELLTYPE} \
  ${CELLTYPE}.bam ${CELLTYPE}_peaks.narrowPeak
```

**Critical:** `--organism=mm10` requires RGT genome data to be downloaded first (for bias correction).

### Motif Matching

```bash
# Match motifs to footprint regions
rgt-motifanalysis matching \
  --organism=mm10 \
  --output-location=./motif_matching/ \
  --input-files ./footprints/${CELLTYPE}.bed
```

### Differential Footprinting (Multi-Condition)

```bash
# Compare multiple cell types simultaneously
rgt-hint differential \
  --organism=mm10 \
  --bc \
  --nc 64 \
  --mpbs-files=./motif_matching/cDC1_mpbs.bed,./motif_matching/cDC2_mpbs.bed,./motif_matching/Mac_mpbs.bed \
  --reads-files=cDC1.bam,cDC2.bam,Mac.bam \
  --conditions=cDC1,cDC2,Mac \
  --output-location=./diff_footprints/ \
  --output-prefix=celltype_comparison
```

**Outputs:**
- `*_statistics.txt`: Tag count and protection score per condition per TF
- `Lineplots/`: ATAC-seq profiles per TF per condition
- Scatter plot of TF activity dynamics

**HINT-ATAC advantage:** Multi-condition comparison produces overview heatmaps directly. TOBIAS `BINDetect` is pairwise only.

### Bias-Corrected Tracks (Genome Browser)

```bash
rgt-hint tracks \
  --bc \
  --bigWig \
  --organism=mm10 \
  ${CELLTYPE}.bam ${CELLTYPE}_peaks.narrowPeak \
  --output-prefix=${CELLTYPE}_BC
```

---

## Signac Footprint Workflow (R — Visualization)

### Setup Motif Information

```r
library(Signac)
library(Seurat)
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)

# Load Seurat object with ATAC assay
obj <- readRDS("integrated_object.rds")

# Get motif PWMs
pwm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 10090, all_versions = FALSE)  # 10090 = Mus musculus
)

# Add motif positions to object (one-time, slow)
obj <- AddMotifs(obj, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pwm)
```

### Compute and Plot Footprints

```r
# Compute footprint data for specific TFs
obj <- Footprint(
  object = obj,
  motif.name = c("IRF8", "BATF::JUN", "SPI1", "ZEB1"),
  genome = BSgenome.Mmusculus.UCSC.mm10,
  in.peaks = TRUE,      # Only motifs within peaks
  upstream = 250,
  downstream = 250
)

# Plot footprints grouped by cell type
p <- PlotFootprint(
  obj,
  features = c("IRF8", "BATF::JUN", "SPI1"),
  group.by = "celltype"
)
p + patchwork::plot_layout(ncol = 1)

# Save
ggsave("footprint_profiles.pdf", width = 8, height = 12)
```

### Interpreting Signac Footprint Plots

**Top panel (Observed):**
- Inverted "V" pattern = active TF binding (flanking enrichment + central depletion)
- Flat profile = no evidence of TF occupancy
- Deeper central valley = stronger/more frequent binding

**Bottom panel (Expected/Background):**
- Shows expected Tn5 cutting from sequence bias alone
- Compare observed vs expected to distinguish true footprints from Tn5 preference artifacts

**Comparing conditions:**
- Higher flanking signal + deeper valley = more active binding in that condition
- Similar profiles across conditions = constitutive binding (not condition-specific)

---

## ChromVAR Integration (Quadrant Analysis)

### Concept

ChromVAR (accessibility) and TOBIAS (occupancy) provide complementary TF activity evidence. Integrate via scatter plot:

| Quadrant | ChromVAR | Footprint | Interpretation |
|----------|----------|-----------|----------------|
| **Q1** | High delta | High delta | **True drivers** — chromatin is open AND TF is bound |
| **Q2** | High delta | Low delta | **Poised/pioneer** — chromatin accessible but TF not occupying (short residence time, indirect effect) |
| **Q3** | Low delta | High delta | Rare — binding in less accessible regions (pioneer factors) or artifact |
| **Q4** | Low delta | Low delta | Not active at these sites |

### Implementation in R

```r
library(tidyverse)

# Load TOBIAS results
tobias <- read_tsv("differential/cDC1_WT_vs_KO/bindetect_results.txt")

# Load ChromVAR differential results
# (Assumes you have per-motif delta deviation scores from your A3 analysis)
chromvar_diff <- read_csv("hub/A3_tf_activity_by_celltype.csv")

# Merge on TF/motif name (requires name harmonization)
combined <- inner_join(
  tobias %>% select(motif_id = output_prefix, tobias_score = WT_Batf3_KO_change,
                    tobias_pval = WT_Batf3_KO_pvalue),
  chromvar_diff %>% select(motif_id, chromvar_delta = delta_deviation),
  by = "motif_id"
)

# Quadrant plot
ggplot(combined, aes(x = chromvar_delta, y = tobias_score)) +
  geom_point(aes(color = tobias_pval < 0.05), alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ggrepel::geom_text_repel(
    data = combined %>% filter(tobias_pval < 0.05),
    aes(label = motif_id), size = 3, max.overlaps = 20
  ) +
  scale_color_manual(values = c("grey70", "red3")) +
  labs(x = "ChromVAR Delta Deviation (Accessibility)",
       y = "TOBIAS Differential Binding Score (Occupancy)",
       title = "TF Activity: Accessibility vs Occupancy") +
  theme_minimal()
```

---

## Motif Database Recommendations

### For Mouse (mm10)

| Database | Version | Mouse TFs | DC Coverage | Format | Notes |
|----------|---------|-----------|-------------|--------|-------|
| **JASPAR** | 2024 Core Vertebrates | ~800+ | Excellent (IRF8, BATF::JUN, SPI1, ETS1, RUNX) | `.jaspar`, `.meme` | Field standard, recommended |
| **HOCOMOCO** | v11 | ~400 mouse | Good (IRF, ETS, RUNX families) | `.meme`, `.pcm` | High quality, curated |
| **CIS-BP** | 2.0 | ~1200+ | Comprehensive | Custom TSV | Broadest coverage, some redundancy |
| **chromVARmotifs** | v1 | Curated | Good | R `PWMatrixList` | Used by chromVAR natively |

**Recommendation:** Use **JASPAR 2024 Core Vertebrates** for TOBIAS. If ChromVAR used a different motif set, export ChromVAR PWMs to MEME format for consistency:

```r
# Export ChromVAR motifs to MEME format for TOBIAS
library(TFBSTools)
library(universalmotif)

# Get the motif PWMs used in ChromVAR
motif_pwms <- GetMotifData(obj, slot = "pwm")  # From Signac object

# Convert to universalmotif and write MEME
motif_list <- lapply(names(motif_pwms), function(nm) {
  convert_motifs(motif_pwms[[nm]], class = "universalmotif-universalmotif")
})
write_meme(motif_list, file = "chromvar_motifs.meme", overwrite = TRUE)
```

**Critical:** Use the **same motif set** for both ChromVAR and TOBIAS to enable direct score correlation. Different PWMs for the same TF across databases yield different binding site calls.

---

## Common Pitfalls

| Issue | Cause | Fix |
|-------|-------|-----|
| TOBIAS finds 0 footprints | Chromosome naming mismatch (`chr1` vs `1`) | Ensure BAM, genome FASTA, and peaks use same naming convention |
| Spurious "perfect" footprints | Missing blacklist filtering | Always use `--blacklist mm10-blacklist.v2.bed` |
| All TFs look differentially bound | Pre-shifted BAMs + TOBIAS internal shift = double-shift | Do NOT pre-shift reads; let TOBIAS handle +4/-5 bp |
| HINT-ATAC bias correction fails | Genome data not downloaded | Run `python -c "from rgt.helper import download_genome; download_genome('mm10')"` |
| Signac Footprint() is very slow | Large fragment files for pseudobulk | Use TOBIAS for quantitative analysis; Signac for selected TFs only |
| ChromVAR vs TOBIAS disagree | Different motif databases | Use identical PWM set for both analyses |
| Low footprint signal | Insufficient sequencing depth in pseudobulk | Pool conditions per celltype to boost depth; check `.info` file read counts |
| Memory errors in ATACorrect | Whole-genome correction at once | Use `--split` flag or process chromosomes separately |

---

## Output Interpretation

### TOBIAS BINDetect

- **Differential binding score > 0**: More TF binding in condition 1 (first `--signals` argument)
- **Differential binding score < 0**: More TF binding in condition 2
- **p-value**: Based on background distribution of all binding sites for that TF
- **Bound/unbound classification**: Per-site BED files in TF-specific subdirectories

### HINT-ATAC Differential

- **Protection score**: Higher = stronger footprint (more TF occupancy)
- **Tag count score**: Number of reads around motif sites (library-size dependent)
- **Statistics file**: Per-TF activity scores across all conditions

### Signac PlotFootprint

- **Flanking enrichment**: Height of the "shoulders" around the motif center
- **Central depletion**: Depth of the "valley" at position 0 (binding site)
- **Footprint depth**: Difference between shoulder height and valley depth

---

## Performance Tips

- TOBIAS ATACorrect is the most time-consuming step (~30-60 min per BAM depending on depth)
- Use `--cores` to parallelize all TOBIAS steps
- HINT-ATAC: Use `--nc` for parallel differential analysis
- For many contrasts, run TOBIAS Steps 1-2 once per BAM, then reuse corrected/scored bigwigs for all pairwise comparisons
- Signac Footprint(): Use `in.peaks = TRUE` to restrict to peak regions (much faster)

---

## Resources

### TOBIAS
- **Publication**: Bentsen et al., Nature Communications 2020 ([DOI: 10.1038/s41467-020-18035-1](https://doi.org/10.1038/s41467-020-18035-1))
- **GitHub**: https://github.com/loosolab/TOBIAS
- **Wiki**: https://github.com/loosolab/TOBIAS/wiki
- **sc-framework (scATAC-seq)**: https://github.com/loosolab/SC-Framework

### HINT-ATAC
- **Publication**: Li et al., Genome Biology 2019 ([DOI: 10.1186/s13059-019-1642-2](https://doi.org/10.1186/s13059-019-1642-2))
- **RGT Suite**: https://reg-gen.readthedocs.io/
- **Tutorial (DC specification)**: https://reg-gen.readthedocs.io/en/latest/hint/tutorial_DC.html

### Signac
- **Footprinting vignette**: https://stuartlab.org/signac/articles/footprint
- **Signac GitHub**: https://github.com/stuart-lab/signac

### Motif Databases
- **JASPAR 2024**: https://jaspar.elixir.no/
- **HOCOMOCO v11**: https://hocomoco11.autosome.org/
- **CIS-BP 2.0**: http://cisbp.ccbr.utoronto.ca/
