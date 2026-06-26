# call_peaks_multistrategy.R
#   Paired-multiome consensus peak calling: per-group pseudobulk + 50/50
#   pseudo-replicate split + MACS3 + reproducibility filter, then the
#   four-strategy cross-strategy reconciliation (CALL -> MERGE -> TIER).
#
# Provenance (reproduced faithfully):
#   /data2/users/JCRLab/JBader/JBader_scHFD/02_analysis/
#     S2_6b_atac_peak_recalling_multistrategies.R
#       create_pseudobulk_bedpe(), call_macs3_peaks(),
#       call_peaks_with_pseudoreps(), call_peaks_by_strategy(),
#       the four-strategy combine + support + iterative merge.
#   Pseudo-rep helpers (split_pseudoreplicates, filter_by_pseudoreps,
#   read_narrowpeak) reproduced from 01_scripts/R/peak_utils.R.
#
# The shared primitives are SOURCED from peak-atlas-framework, not re-defined:
#   iterative_overlap.R  -> clusterGRanges, convergeClusterGRanges
#   support_voting.R     -> calculate_strategy_support, add_adjusted_score
#   normalize_width.R    -> normalize_to_501bp
#   blacklist.R          -> load_blacklist, remove_blacklist_peaks
#
# Project-specific genome / MACS path / chromosome set are FUNCTION ARGUMENTS
# with mouse mm39 defaults (the source project). For human, pass
# macs_genome = "hs", std_chroms = paste0("chr", c(1:22, "X", "Y")), and an
# hg38 blacklist BED.

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(IRanges)
})

# ---- source the framework primitives -----------------------------------------
# Adjust FRAMEWORK_SCRIPTS to wherever peak-atlas-framework/scripts lives.
FRAMEWORK_SCRIPTS <- Sys.getenv(
  "PEAK_ATLAS_FRAMEWORK_SCRIPTS",
  unset = "../../peak-atlas-framework/scripts"
)
for (f in c("iterative_overlap.R", "support_voting.R",
            "normalize_width.R", "blacklist.R")) {
  src <- file.path(FRAMEWORK_SCRIPTS, f)
  if (file.exists(src)) source(src) else
    warning("Framework script not found (source it manually): ", src)
}

# ==============================================================================
# Configuration (source: S2_6b config block)
# ==============================================================================
MIN_CELLS_PER_GROUP      <- 50    # need 2x for two pseudo-reps
MAX_CELLS_PER_GROUP      <- 1000  # downsample large groups for speed
MIN_FRAGS_PER_PSEUDOREP  <- 5000  # skip a rep below this many fragments
PEAK_EXTEND              <- 250   # summit +/- 250 = 501bp

# scATAC-optimized MACS3 params (Corces & Granja / ENCODE); see framework ref.
MACS3_PARAMS <- "--nomodel --shift -100 --extsize 200 --call-summits --keep-dup all -q 0.01"

# ==============================================================================
# Pseudo-rep helpers (reproduced from peak_utils.R)
# ==============================================================================

#' Split cells into two pseudo-replicates by a random 50/50 split.
split_pseudoreplicates <- function(cells, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- length(cells)
  n_rep1 <- floor(n / 2)
  rep1_idx <- sample(n, n_rep1)
  list(rep1 = cells[rep1_idx], rep2 = cells[-rep1_idx])
}

#' Keep only rep1 peaks that overlap a rep2 peak (PRIMARY noise filter).
filter_by_pseudoreps <- function(peaks_rep1, peaks_rep2,
                                 min_overlap = 1, verbose = TRUE) {
  if (length(peaks_rep1) == 0 || length(peaks_rep2) == 0) {
    if (verbose) message("    Empty pseudo-rep, returning empty GRanges")
    return(GRanges())
  }
  reproducible_mask <- overlapsAny(peaks_rep1, peaks_rep2, minoverlap = min_overlap)
  if (verbose) {
    message(sprintf("    Rep1: %d, Rep2: %d, Reproducible: %d (%.1f%% retained)",
                    length(peaks_rep1), length(peaks_rep2), sum(reproducible_mask),
                    100 * sum(reproducible_mask) / max(length(peaks_rep1), 1)))
  }
  peaks_rep1[reproducible_mask]
}

#' Read a MACS3 narrowPeak file into GRanges (0-based -> 1-based).
read_narrowpeak <- function(file) {
  if (!file.exists(file)) stop("File not found: ", file)
  df <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(df) <- c("chr", "start", "end", "name", "score", "strand",
                    "signalValue", "pValue", "qValue", "peak")
  GRanges(
    seqnames = df$chr,
    ranges = IRanges(start = df$start + 1, end = df$end),
    strand = "*",
    score = df$score, signalValue = df$signalValue,
    pValue = df$pValue, qValue = df$qValue,
    peak = df$peak  # summit offset from start
  )
}

# ==============================================================================
# Pseudobulk + MACS3 (reproduced from S2_6b)
# ==============================================================================

#' Write a pseudobulk BEDPE for `cells` by filtering per-sample fragment files.
#' Strips a barcode prefix of the form SAMPLE_<barcode>-<n>. Returns n fragments.
create_pseudobulk_bedpe <- function(obj, cells, output_file, atac_assay = "ATAC") {
  frag_objs <- Fragments(obj[[atac_assay]])
  if (length(frag_objs) == 0) stop("No fragment files found in ", atac_assay, " assay")

  cell_df <- data.frame(full_barcode = cells, stringsAsFactors = FALSE)
  cell_df$sample <- sub("_[ACGTN]+-[0-9]+$", "", cells)
  cell_df$core_barcode <- sub("^.*_([ACGTN]+-[0-9]+)$", "\\1", cells)

  if (file.exists(output_file)) unlink(output_file)
  system(sprintf("touch '%s'", output_file))

  for (i in seq_along(frag_objs)) {
    frag_path <- GetFragmentData(frag_objs[[i]], slot = "path")
    if (!file.exists(frag_path)) {
      warning(sprintf("Fragment file not found: %s", frag_path)); next
    }
    # Infer the sample name from the 10x outs/ path; fall back to the parent dir.
    path_parts <- strsplit(frag_path, "/")[[1]]
    outs_idx <- which(path_parts == "outs")
    frag_sample <- if (length(outs_idx) > 0) path_parts[outs_idx[1] - 1]
                   else basename(dirname(dirname(frag_path)))

    sample_cells <- cell_df$core_barcode[cell_df$sample == frag_sample]
    if (length(sample_cells) == 0) next

    barcode_file <- tempfile(pattern = paste0("barcodes_", frag_sample, "_"),
                             fileext = ".txt")
    writeLines(sample_cells, barcode_file)
    # Expand each fragment by its duplicate count ($5) into a 1bp-anchored BEDPE.
    cmd <- sprintf(
      "zcat '%s' | grep -v '^#' | awk 'NR==FNR{barcodes[$1]=1; next} $4 in barcodes {for(c=1;c<=$5;c++) print $1\"\\t\"$2\"\\t\"$2+1\"\\t\"$1\"\\t\"$3-1\"\\t\"$3}' '%s' - >> '%s'",
      frag_path, barcode_file, output_file
    )
    system(cmd, wait = TRUE)
    unlink(barcode_file)
  }
  as.numeric(system(sprintf("wc -l < '%s'", output_file), intern = TRUE))
}

#' Call MACS3 on a BEDPE and return the narrowPeak as GRanges.
call_macs3_peaks <- function(bedpe_file, output_name, macs_dir,
                             macs3_path = "macs3", macs_genome = "mm",
                             macs_params = MACS3_PARAMS) {
  peaks_file <- file.path(macs_dir, paste0(output_name, "_peaks.narrowPeak"))
  macs_cmd <- sprintf(
    "%s callpeak -t '%s' -g %s -f BEDPE -n '%s' --outdir '%s' %s 2>&1",
    macs3_path, bedpe_file, macs_genome, output_name, macs_dir, macs_params
  )
  system(macs_cmd, intern = TRUE)
  if (!file.exists(peaks_file)) {
    warning(sprintf("MACS3 did not produce output for %s", output_name)); return(NULL)
  }
  read_narrowpeak(peaks_file)
}

# ==============================================================================
# Per-cluster pseudo-replicate caller (core Corces/Granja routine)
# ==============================================================================

#' Call reproducible peaks for one cell group with pseudo-replicate filtering.
#' @return GRanges of reproducible 501bp peaks (with group/strategy/n_cells), or NULL.
call_peaks_with_pseudoreps <- function(obj, cells, group_name, strategy_name,
                                       macs_dir, pseudorep_dir,
                                       macs3_path = "macs3", macs_genome = "mm",
                                       atac_assay = "ATAC", seed = 1234567) {
  message(sprintf("    Processing %s (%d cells)...", group_name, length(cells)))

  if (length(cells) < MIN_CELLS_PER_GROUP * 2) {
    message(sprintf("      Skipping: need %d cells for pseudo-reps, have %d",
                    MIN_CELLS_PER_GROUP * 2, length(cells)))
    return(NULL)
  }
  if (length(cells) > MAX_CELLS_PER_GROUP) {
    set.seed(seed)
    cells <- sample(cells, MAX_CELLS_PER_GROUP)
    message(sprintf("      Downsampled to %d cells", MAX_CELLS_PER_GROUP))
  }

  output_base <- gsub("[^A-Za-z0-9_]", "_", paste0(strategy_name, "_", group_name))

  # 1. Split 50/50 into pseudo-reps.
  set.seed(seed)
  pseudoreps <- split_pseudoreplicates(cells)
  message(sprintf("      Pseudo-reps: %d + %d cells",
                  length(pseudoreps$rep1), length(pseudoreps$rep2)))

  # 2. Pseudobulk BEDPE per rep.
  bedpe_rep1 <- file.path(pseudorep_dir, paste0(output_base, "_rep1.bedpe"))
  bedpe_rep2 <- file.path(pseudorep_dir, paste0(output_base, "_rep2.bedpe"))
  n_frags_rep1 <- create_pseudobulk_bedpe(obj, pseudoreps$rep1, bedpe_rep1, atac_assay)
  n_frags_rep2 <- create_pseudobulk_bedpe(obj, pseudoreps$rep2, bedpe_rep2, atac_assay)
  message(sprintf("      Fragments: rep1=%s, rep2=%s",
                  format(n_frags_rep1, big.mark = ","),
                  format(n_frags_rep2, big.mark = ",")))

  if (n_frags_rep1 < MIN_FRAGS_PER_PSEUDOREP || n_frags_rep2 < MIN_FRAGS_PER_PSEUDOREP) {
    message(sprintf("      Skipping: insufficient fragments (min: %d)", MIN_FRAGS_PER_PSEUDOREP))
    unlink(c(bedpe_rep1, bedpe_rep2)); return(NULL)
  }

  # 3. MACS3 per rep.
  peaks_rep1_raw <- call_macs3_peaks(bedpe_rep1, paste0(output_base, "_rep1"),
                                     macs_dir, macs3_path, macs_genome)
  peaks_rep2_raw <- call_macs3_peaks(bedpe_rep2, paste0(output_base, "_rep2"),
                                     macs_dir, macs3_path, macs_genome)
  unlink(c(bedpe_rep1, bedpe_rep2))
  if (is.null(peaks_rep1_raw) || is.null(peaks_rep2_raw)) {
    message("      MACS3 failed for one or both pseudo-reps"); return(NULL)
  }

  # 4. Normalize each rep to 501bp on the summit (framework helper).
  peaks_rep1 <- normalize_to_501bp(peaks_rep1_raw, extend = PEAK_EXTEND)
  peaks_rep2 <- normalize_to_501bp(peaks_rep2_raw, extend = PEAK_EXTEND)

  # 5. PRIMARY NOISE FILTER: keep rep1 peaks overlapping rep2.
  reproducible_peaks <- filter_by_pseudoreps(peaks_rep1, peaks_rep2, verbose = TRUE)
  if (length(reproducible_peaks) == 0) {
    message("      No reproducible peaks found"); return(NULL)
  }

  mcols(reproducible_peaks)$group    <- group_name
  mcols(reproducible_peaks)$strategy <- strategy_name
  mcols(reproducible_peaks)$n_cells  <- length(cells)
  reproducible_peaks
}

#' Call peaks for every group under one grouping variable (one strategy).
#' @return named list of GRanges (one per group that passed the guardrails).
call_peaks_by_strategy <- function(obj, grouping_var, strategy_name,
                                   macs_dir, pseudorep_dir,
                                   macs3_path = "macs3", macs_genome = "mm",
                                   atac_assay = "ATAC") {
  message(sprintf("\n=== Strategy %s: calling peaks by %s ===", strategy_name, grouping_var))
  groups <- unique(obj@meta.data[[grouping_var]])
  groups <- groups[!is.na(groups)]

  peak_list <- list()
  for (grp in groups) {
    grp_str <- as.character(grp)
    cells <- colnames(obj)[obj@meta.data[[grouping_var]] == grp]
    peaks_grp <- call_peaks_with_pseudoreps(
      obj, cells, grp_str, strategy_name,
      macs_dir = macs_dir, pseudorep_dir = pseudorep_dir,
      macs3_path = macs3_path, macs_genome = macs_genome, atac_assay = atac_assay
    )
    if (!is.null(peaks_grp) && length(peaks_grp) > 0) peak_list[[grp_str]] <- peaks_grp
  }
  total_peaks <- sum(vapply(peak_list, length, integer(1)))
  message(sprintf("  Strategy %s complete: %d groups, %s reproducible peaks",
                  strategy_name, length(peak_list), format(total_peaks, big.mark = ",")))
  peak_list
}

# ==============================================================================
# Four-strategy reconciliation: MERGE within, support BEFORE final merge, TIER.
# ==============================================================================

#' Reconcile the four strategies into one fixed-width consensus atlas.
#'
#' @param all_strategy_peaks Named list of length 4 (RNA, ATAC, WNN, CellType),
#'   each itself a named list of per-group GRanges (output of call_peaks_by_strategy).
#' @param blacklist GRanges of blacklist regions (load via blacklist.R).
#' @param BSgenome BSgenome object for the chromosome-boundary filter.
#' @param std_chroms Standard chromosomes to keep (default mouse).
#' @return GRanges consensus atlas with score/adjusted_score/n_strategies/strategies/peak_id.
reconcile_strategies <- function(all_strategy_peaks, blacklist, BSgenome,
                                 std_chroms = paste0("chr", c(1:19, "X", "Y"))) {

  # 1. Merge WITHIN each strategy (rank by raw score), after QC filtering.
  strategy_merged <- list()
  for (strat in names(all_strategy_peaks)) {
    strat_peaks <- do.call(c, base::unname(all_strategy_peaks[[strat]]))
    if (length(strat_peaks) == 0) next
    strat_peaks <- keepSeqlevels(strat_peaks, std_chroms, pruning.mode = "coarse")
    strat_peaks <- remove_blacklist_peaks(strat_peaks, blacklist)
    sm <- convergeClusterGRanges(strat_peaks, by = "score", decreasing = TRUE)
    mcols(sm)$strategy_source <- strat
    strategy_merged[[strat]] <- sm
    message(sprintf("  %s: %s -> %s peaks", strat,
                    format(length(strat_peaks), big.mark = ","),
                    format(length(sm), big.mark = ",")))
  }

  # 2. Combine the four merged sets.
  combined <- do.call(c, base::unname(strategy_merged))

  # 3. Boundary filter (drop peaks running past a chromosome end).
  chrom_lengths <- seqlengths(BSgenome)
  pc <- as.character(seqnames(combined))
  valid <- !is.na(chrom_lengths[pc]) & start(combined) >= 1 & end(combined) <= chrom_lengths[pc]
  combined <- combined[valid]

  # 4. Support BEFORE the final merge (else "bridging" => everything reads n=4).
  combined <- calculate_strategy_support(combined, strategy_merged)
  message("Strategy support distribution (BEFORE merge):")
  print(table(mcols(combined)$n_strategies))

  # 5. Support-weighted score boost, then the final iterative-overlap merge.
  combined <- add_adjusted_score(combined)  # adjusted = score*(1+0.5*(n-1))
  final <- convergeClusterGRanges(combined, by = "adjusted_score",
                                  decreasing = TRUE, verbose = TRUE)

  # Recompute support on the final survivors for reporting (informational).
  final <- calculate_strategy_support(final, strategy_merged)

  mcols(final)$peak_id <- paste0(seqnames(final), ":", start(final), "-", end(final))
  names(final) <- mcols(final)$peak_id
  final
}

# ------------------------------------------------------------------------------
# Sketch of a full run (the four strategies, then reconcile):
#   blacklist <- load_blacklist(bed_path = "mm39-blacklist.bed")
#   common <- list(macs_dir = macs_dir, pseudorep_dir = prep_dir,
#                   macs3_path = "/path/to/macs3", macs_genome = "mm")
#   s_rna  <- do.call(call_peaks_by_strategy, c(list(obj, "RNA_clusters",  "RNA"),      common))
#   s_atac <- do.call(call_peaks_by_strategy, c(list(obj, "ATAC_clusters", "ATAC"),     common))
#   s_wnn  <- do.call(call_peaks_by_strategy, c(list(obj, "WNN_clusters",  "WNN"),      common))
#   s_ct   <- do.call(call_peaks_by_strategy, c(list(obj, "refined_cell_type","CellType"), common))
#   atlas  <- reconcile_strategies(
#               list(RNA = s_rna, ATAC = s_atac, WNN = s_wnn, CellType = s_ct),
#               blacklist, BSgenome.Mmusculus.UCSC.mm39)
# ------------------------------------------------------------------------------
