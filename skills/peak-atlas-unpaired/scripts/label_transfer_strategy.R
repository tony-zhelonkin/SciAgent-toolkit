# label_transfer_strategy.R - Unpaired bridge + staged CCA label transfer +
#                             per-strategy (A/B/C) Signac peak calling.
#
# Provenance (reproduced faithfully from the DC_Dictionary foundation layer):
#   GeneActivity bridge       -> foundation/S3_atac_qc_normalization.R (Step 6)
#   R1 coarse transfer        -> foundation/S5.1a_integration_r1.R
#   R2 refined transfer       -> foundation/S5.2a_integration_r2.R
#   Strategy A (per celltype) -> foundation/S5.1b_call_strategy_a_peaks.R
#   Strategy B (subtype x cond)-> foundation/S5.2b_call_strategy_b_peaks.R
#   Strategy C (label-free)   -> foundation/S4.1c_call_strategy_c_peaks.R
#
# The unpaired regime has no shared barcodes between RNA and ATAC. Identity is
# transferred from an annotated RNA reference through a gene-activity surrogate
# assay (the only shared feature space). Genome/path specifics from the source
# project (mouse mm10/mm39, effective.genome.size 1.87e9, macs3 path) are
# function arguments with sensible defaults - override for your build.

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(GenomicRanges)
  library(IRanges)
})

# ---------------------------------------------------------------------------
# 1. Bridge: gene-activity surrogate RNA assay on the ATAC pool
#    (provenance: S3_atac_qc_normalization.R Step 6)
# ---------------------------------------------------------------------------

#' Build the gene-activity bridge assay on an ATAC object
#'
#' Sums fragments over gene bodies + promoters into a gene-level "RNA" assay -
#' the shared feature space anchoring RNA and ATAC. Requires an Annotation()
#' on the ATAC assay (set during Signac object construction).
#'
#' @param atac Seurat ATAC object with a fragments-backed peak assay.
#' @param assay_name Name for the new gene-activity assay (default "RNA").
#' @return atac with the gene-activity assay added and LogNormalized.
add_gene_activity_bridge <- function(atac, assay_name = "RNA") {
  gene.activities <- GeneActivity(atac)
  atac[[assay_name]] <- CreateAssayObject(counts = gene.activities)
  atac <- NormalizeData(
    object = atac,
    assay = assay_name,
    normalization.method = "LogNormalize",
    scale.factor = median(atac[[paste0("nCount_", assay_name)]][, 1])
  )
  atac
}

# ---------------------------------------------------------------------------
# 2. Staged CCA label transfer (RNA -> ATAC in gene-activity space)
#    (provenance: S5.1a_integration_r1.R, S5.2a_integration_r2.R)
# ---------------------------------------------------------------------------

#' Transfer labels from an annotated RNA reference onto the ATAC pool
#'
#' Anchors in gene-activity space (reduction = "cca"), weights the transfer by
#' the ATAC LSI reduction, then marks cells below `prediction_threshold` (or
#' outside `target_labels`) as "LowConf" so they can be excluded from
#' label-resolved peak calling.
#'
#' @param rna_ref Annotated RNA Seurat object; labels in `label_col`.
#' @param atac    ATAC Seurat object with a gene-activity assay + "lsi".
#' @param label_col Column in rna_ref@meta.data holding the labels.
#' @param target_labels Character vector of accepted labels. Predictions not in
#'   this set are demoted to "LowConf". NULL accepts any predicted label.
#' @param label_map Optional named vector recoding fine RNA labels into coarse
#'   labels (R1 use). NULL transfers labels directly (R2 use).
#' @param prediction_threshold Min prediction.score.max to keep a label
#'   (R1 default 0.50; pass 0.55 for the stricter R2 round).
#' @param dims LSI/CCA dimensions (default 2:30).
#' @param k_weight TransferData k.weight (default 100).
#' @param nfeatures Variable features in the reference (R1 3000, R2 4000).
#' @param prefix Metadata prefix for the prediction columns ("r1"/"r2").
#' @return atac with prediction metadata and a `<prefix>_highconf` label column.
transfer_labels <- function(rna_ref, atac, label_col = "cell_annotations",
                            target_labels = NULL, label_map = NULL,
                            prediction_threshold = 0.50, dims = 2:30,
                            k_weight = 100, nfeatures = 3000, prefix = "r1") {
  DefaultAssay(rna_ref) <- "RNA"
  DefaultAssay(atac) <- "RNA"   # gene activity = shared space for anchoring

  rna_ref <- FindVariableFeatures(rna_ref, nfeatures = nfeatures)
  var_features <- VariableFeatures(rna_ref)

  anchors <- FindTransferAnchors(
    reference = rna_ref, query = atac,
    features = var_features,
    reference.assay = "RNA", query.assay = "RNA",
    reduction = "cca", dims = dims
  )

  predictions <- TransferData(
    anchorset = anchors,
    refdata = rna_ref[[label_col]][, 1],
    weight.reduction = atac[["lsi"]],
    dims = dims, k.weight = k_weight
  )

  DefaultAssay(atac) <- "peaks"   # revert for downstream ATAC work

  pred <- as.data.frame(predictions)
  colnames(pred) <- gsub("\\.", "_", colnames(pred))

  # Optional coarse recode (R1); otherwise transfer the predicted id directly
  resolved <- if (!is.null(label_map)) {
    dplyr::recode(pred$predicted_id, !!!label_map, .default = "Other")
  } else {
    pred$predicted_id
  }

  # High-confidence gate: score >= threshold AND label in target set
  keep <- pred$prediction_score_max >= prediction_threshold
  if (!is.null(target_labels)) keep <- keep & resolved %in% target_labels
  highconf <- ifelse(keep, resolved, "LowConf")

  pred_prefixed <- pred
  colnames(pred_prefixed) <- paste0(prefix, "_", colnames(pred_prefixed))
  pred_prefixed[[paste0(prefix, "_label")]] <- resolved
  pred_prefixed[[paste0(prefix, "_highconf")]] <- highconf

  atac <- AddMetaData(atac, metadata = pred_prefixed)
  atac
}

# ---------------------------------------------------------------------------
# 3. Shared per-group post-process (identical across A/B/C)
#    (provenance: S5.1b / S5.2b / S4.1c - the "Process to 501bp" blocks)
# ---------------------------------------------------------------------------

#' Normalize a single group's MACS peaks to 501bp and apply CALL-stage filters
#'
#' summit -> 501bp resize -> standard chromosomes -> drop blacklist ->
#' drop N-rich (N fraction >= n_max).
#'
#' @param peaks_gr GRanges from Signac::CallPeaks (summit in mcols$peak).
#' @param blacklist GRanges blacklist (empty GRanges = skip).
#' @param bsgenome BSgenome object for the N-content filter.
#' @param std_chr Standard chromosomes to keep (default mouse chr1..19,X,Y).
#' @param n_max N-fraction cutoff; peaks at/above are dropped (default 0.1).
#' @return GRanges of 501bp summit-centered, filtered peaks.
process_group_peaks <- function(peaks_gr, blacklist, bsgenome,
                                 std_chr = paste0("chr", c(1:19, "X", "Y")),
                                 n_max = 0.1) {
  if ("peak" %in% names(mcols(peaks_gr))) {
    summit_pos <- start(peaks_gr) + mcols(peaks_gr)$peak
    peaks_gr <- GRanges(seqnames = seqnames(peaks_gr),
                        ranges = IRanges(start = summit_pos, width = 1),
                        strand = "*")
  }
  peaks_gr <- resize(peaks_gr, width = 501, fix = "center")
  peaks_gr <- keepSeqlevels(peaks_gr, std_chr, pruning.mode = "coarse")
  if (length(blacklist) > 0) {
    peaks_gr <- peaks_gr[!overlapsAny(peaks_gr, blacklist)]
  }
  peak_seqs <- Biostrings::getSeq(bsgenome, peaks_gr)
  n_content <- Biostrings::letterFrequency(peak_seqs, "N") / 501
  peaks_gr[n_content < n_max]
}

# ---------------------------------------------------------------------------
# 4. Call peaks for one strategy across its groups, then merge
#    (provenance: S5.1b / S5.2b / S4.1c - identical structure, different group)
# ---------------------------------------------------------------------------

#' Call MACS3 peaks per group for one strategy and merge to a scored union
#'
#' For each group with >= min_cells cells, calls Signac::CallPeaks WITH
#' group.by (critical - without it MACS sees all fragments and every group
#' returns identical peaks), post-processes to 501bp, then merges across groups
#' with reduce(min.gapwidth = 50) and scores each merged peak by the number of
#' groups that support it (support count).
#'
#' @param atac ATAC object (already subset to high-confidence cells for A/B).
#' @param group_col Metadata column defining the groups (e.g. "r1_highconf"
#'   for A, "cluster_condition" for B, "orig.ident" for label-free C).
#' @param blacklist,bsgenome Passed to process_group_peaks().
#' @param min_cells Minimum cells per group (A/C = 100, B = 20).
#' @param effective_genome_size MACS effective genome size (mouse 1.87e9).
#' @param macs_path Path to the macs3 binary.
#' @param outdir,strategy Output dir and strategy tag ("A"/"B"/"C").
#' @return GRanges union with `support_count`, `strategy`, `groups` metadata.
call_strategy_peaks <- function(atac, group_col, blacklist, bsgenome,
                                min_cells = 100,
                                effective_genome_size = 1.87e9,
                                macs_path = "macs3",
                                outdir = tempdir(), strategy = "A") {
  grp <- atac[[group_col]][, 1]
  counts <- table(grp)
  groups <- names(counts)[counts >= min_cells]
  if (length(groups) == 0) stop("No group has >= ", min_cells, " cells in ", group_col)

  per_group <- list()
  for (g in groups) {
    peaks_gr <- tryCatch(
      CallPeaks(
        object = atac,
        group.by = group_col,   # CRITICAL: triggers Signac fragment filtering
        idents = g,
        macs2.path = macs_path,
        outdir = outdir,
        name = paste0(strategy, "_", make.names(g)),
        format = "BEDPE",
        effective.genome.size = effective_genome_size,
        additional.args = "-q 0.01 --call-summits --nolambda --keep-dup all",
        cleanup = FALSE
      ),
      error = function(e) { warning("CallPeaks failed for ", g, ": ", e$message); NULL }
    )
    if (is.null(peaks_gr) || length(peaks_gr) == 0) next
    pg <- process_group_peaks(peaks_gr, blacklist, bsgenome)
    mcols(pg)$group <- g
    mcols(pg)$strategy <- strategy
    per_group[[g]] <- pg
  }
  if (length(per_group) == 0) stop("No peaks called for strategy ", strategy)

  # Merge across groups (reduce gapwidth=50), re-resize to 501bp
  union_gr <- reduce(do.call(c, base::unname(per_group)), min.gapwidth = 50)
  union_gr <- resize(union_gr, width = 501, fix = "center")

  # Support count = number of groups overlapping each merged peak (the score)
  support_matrix <- sapply(per_group, function(p) overlapsAny(union_gr, p))
  mcols(union_gr)$support_count <- rowSums(support_matrix)
  mcols(union_gr)$groups <- apply(support_matrix, 1, function(r)
    paste(names(per_group)[r], collapse = ","))
  mcols(union_gr)$score_max <- mcols(union_gr)$support_count
  mcols(union_gr)$strategy <- strategy
  union_gr
}

# ---------------------------------------------------------------------------
# Canonical order (unpaired pipeline):
#   atac <- add_gene_activity_bridge(atac)
#   atac <- transfer_labels(rna_whole,    atac, target_labels = coarse, prefix = "r1",
#                           label_map = coarse_map, prediction_threshold = 0.50)
#   atac <- transfer_labels(rna_focused,  atac, target_labels = refined, prefix = "r2",
#                           prediction_threshold = 0.55, nfeatures = 4000)
#   atac_A <- subset(atac, r1_highconf != "LowConf")
#   atac_B <- subset(atac, r2_highconf != "LowConf")
#   atac_B$cluster_condition <- paste0(atac_B$r2_highconf, "_", atac_B$orig.ident)
#   A <- call_strategy_peaks(atac_A, "r1_highconf",      bl, bsg, min_cells = 100, strategy = "A")
#   B <- call_strategy_peaks(atac_B, "cluster_condition",bl, bsg, min_cells = 20,  strategy = "B")
#   C <- call_strategy_peaks(atac,   "orig.ident",       bl, bsg, min_cells = 100, strategy = "C")  # label-free
# ---------------------------------------------------------------------------
