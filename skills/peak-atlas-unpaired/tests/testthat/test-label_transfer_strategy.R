# Contract for the bridge + transfer + per-strategy calling
# (scripts/label_transfer_strategy.R). These depend on Seurat / Signac::CallPeaks /
# BSgenome, so they stay documented scaffolds with mocks. The notes ARE the spec.
source_script("label_transfer_strategy.R")

test_that("process_group_peaks emits 501bp peaks on standard chromosomes only", {
  skip_scaffold("the N-content filter needs a BSgenome; mock Biostrings::getSeq or pass a toy genome")
  # peaks <- GRanges with a 'peak' summit-offset col, some on chrUn / chrM
  # out <- process_group_peaks(peaks, blacklist = GRanges(), bsgenome = toy_bsgenome)
  # expect_true(all(width(out) == 501L))
  # expect_true(all(as.character(seqnames(out)) %in% paste0("chr", c(1:19, "X", "Y"))))
})

test_that("call_strategy_peaks sets group.by, gates on min_cells, scores by #groups", {
  skip_scaffold("mock Signac::CallPeaks (mock_call_peaks) so MACS3 is not required")
  # With a mock CallPeaks returning a known per-group GRanges:
  #   - groups with < min_cells are skipped (B uses 20; A/C use 100)
  #   - if NO group clears min_cells -> the function stop()s
  #   - support_count of a merged peak == number of contributing groups
  #   - cross-group merge uses reduce(min.gapwidth = 50) then resize to 501bp
})

test_that("transfer_labels demotes low-confidence / off-target predictions to LowConf", {
  skip_scaffold("needs Seurat refs + FindTransferAnchors; mock TransferData output")
  # Given a mocked predictions frame:
  #   - prediction_score_max < prediction_threshold -> '<prefix>_highconf' == 'LowConf'
  #   - predicted label not in target_labels         -> 'LowConf'
  #   - label_map recodes fine -> coarse labels (R1); NULL passes labels through (R2)
})

test_that("add_gene_activity_bridge adds a normalized gene-level assay", {
  skip_scaffold("needs Signac GeneActivity() + an Annotation()-bearing ATAC object")
  succeed()  # placeholder so the scaffold is a valid, skipped test
})
