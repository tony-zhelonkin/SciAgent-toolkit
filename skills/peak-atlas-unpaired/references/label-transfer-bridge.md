# The gene-activity bridge and staged CCA label transfer

Unpaired data has no shared barcodes: the scRNA cells and scATAC cells are
different cells from the same biological system. Cell identity therefore cannot
be observed on the ATAC side — it must be **transferred** from the annotated
RNA reference. This file documents the bridge that makes transfer possible and
the two-round staging used in the source pipeline.

## Why a bridge is needed

RNA features are genes; ATAC features are peaks. The two assays live in
different coordinate spaces, so there is no direct anchor between them. The
**gene-activity matrix** converts ATAC into a gene-level surrogate "RNA" assay
by summing fragments over each gene body plus its promoter. That surrogate
assay shares its feature space (gene names) with the RNA reference, so anchoring
happens there.

```r
# Provenance: S3_atac_qc_normalization.R, Step 6
gene.activities <- GeneActivity(atac)                       # Signac
atac[["RNA"]]   <- CreateAssayObject(counts = gene.activities)
atac <- NormalizeData(atac, assay = "RNA",
                      normalization.method = "LogNormalize",
                      scale.factor = median(atac$nCount_RNA))
```

`GeneActivity()` needs an `Annotation()` on the ATAC peak assay (gene models for
the build), set when the Signac object is constructed. The result is a noisy but
serviceable proxy for expression — good enough to anchor lineage identity, not
good enough to read off subtle states (which is why thresholds matter).

## The transfer call

Anchor in gene-activity space, then weight the transfer by the ATAC's own LSI so
that the prediction respects ATAC neighborhood structure:

```r
DefaultAssay(rna_ref) <- "RNA"; DefaultAssay(atac) <- "RNA"
rna_ref <- FindVariableFeatures(rna_ref, nfeatures = 3000)   # 4000 for R2
anchors <- FindTransferAnchors(
  reference = rna_ref, query = atac,
  features = VariableFeatures(rna_ref),
  reference.assay = "RNA", query.assay = "RNA",
  reduction = "cca", dims = 2:30)
pred <- TransferData(
  anchorset = anchors, refdata = rna_ref$cell_annotations,
  weight.reduction = atac[["lsi"]], dims = 2:30, k.weight = 100)
DefaultAssay(atac) <- "peaks"   # revert for downstream ATAC work
```

Key choices (provenance `S5.1a`/`S5.2a`):
- **`reduction = "cca"`** — canonical correlation analysis finds the shared
  axes of variation between the two assays; the right anchoring method for
  cross-modality transfer in gene-activity space.
- **`weight.reduction = atac[["lsi"]]`** — the transfer weights come from the
  ATAC LSI, not the (noisy) gene-activity assay, so predictions vary smoothly
  over the ATAC manifold.
- **`dims = 2:30`** — LSI component 1 typically captures sequencing depth, so it
  is dropped (start at 2), as in standard Signac LSI usage.
- **`k.weight = 100`** — anchor weighting; higher is stricter.

## Two rounds: coarse, then refined

A single transfer of fine labels onto noisy gene activity is unreliable.
Staging it stabilizes the result:

| Round | Reference | nfeatures | Labels | Threshold | Output column | Drives |
|---|---|---|---|---|---|---|
| **R1** | whole RNA reference | 3000 | COARSE lineage (recode fine labels via a `coarse_map`, demote off-target to "Other") | `>= 0.50` | `r1_highconf` | Strategy A |
| **R2** | focused RNA reference (target lineages only) | 4000 | REFINED subtype (transferred directly) | `>= 0.55` | `r2_highconf` | Strategy B |

R2 runs on the **R1-labeled ATAC** — coarse lineage first, then resolve
subtypes within. The reference is narrowed to the lineages of interest so the
fine-label transfer is not distracted by irrelevant cell types. The stricter
0.55 threshold reflects that fine labels carry more transfer error than coarse
ones.

## LowConf — keeping bad labels out of peaks

A cell is **LowConf** when its `prediction.score.max` is below the round's
threshold OR its predicted label is not in the accepted target set. LowConf
cells are **excluded from label-resolved peak calling**:

```r
atac_A <- subset(atac, r1_highconf != "LowConf")   # Strategy A input
atac_B <- subset(atac, r2_highconf != "LowConf")   # Strategy B input
```

Why: a cell whose lineage is uncertain would smear its fragments into the wrong
group's pseudobulk, contaminating that lineage's peaks. The label-free Strategy
C is exempt — it groups only by condition and uses no transferred labels, so
LowConf cells still contribute there (and only there).

## Pitfalls

| Symptom | Cause | Fix |
|---|---|---|
| Anchors near zero / nonsense transfer | Anchored on the peak assay, not gene activity | Set `DefaultAssay` to the gene-activity ("RNA") assay on BOTH objects before `FindTransferAnchors` |
| LSI component 1 dominates, depth artifacts | `dims` started at 1 | Use `dims = 2:30` (drop the depth component) |
| Fine labels look random | Transferred fine labels in one shot on whole reference | Stage: R1 coarse on whole ref, then R2 refined on a focused ref over the R1-labeled ATAC |
| Lineage peaks contaminated | LowConf cells included in Strategy A/B | `subset(... != "LowConf")` before label-resolved calling |

## See also

- `references/abc-strategies.md` — how the transferred labels drive A and B.
- `peak-atlas-framework/references/macs3-scatac-params.md` — the CallPeaks
  wrapper and 501bp/QC filters that consume these groupings.
- `scripts/label_transfer_strategy.R` — `add_gene_activity_bridge`,
  `transfer_labels`.
