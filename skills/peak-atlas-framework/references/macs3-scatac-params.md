# scATAC MACS3 parameters and the CALL-stage QC filters

The CALL stage turns a per-group pseudobulk into summit-centered fixed-width peaks. This file documents the exact MACS3 invocation, why each flag matters for scATAC, the 501bp normalization, and the chromosome / blacklist / N-content filters that clean the raw calls.

## Canonical scATAC MACS3 call

The pipeline calls MACS3 directly on a per-pseudorep BEDPE:

```bash
macs3 callpeak \
  -t <pseudobulk.bedpe> \
  -g mm \
  -f BEDPE \
  --nomodel \
  --shift -100 \
  --extsize 200 \
  --call-summits \
  --keep-dup all \
  -q 0.01
```

(Provenance: `MACS3_PARAMS` in `02_analysis/S2_6b_atac_peak_recalling_multistrategies.R` and `01_scripts/R/README_peak_calling.md`. Use `-g hs` and a human genome for human data.)

### Why each flag

| Flag | Why for scATAC |
|---|---|
| `-f BEDPE` | Pseudobulk fragments are paired-end intervals; BEDPE preserves both Tn5 ends. |
| `--nomodel` | ATAC fragments are too short / shifted to fit MACS's shifting model — skip model building. |
| `--shift -100 --extsize 200` | Re-center signal on the Tn5 cut site: shift each read 100bp upstream and extend 200bp, so the pileup peaks at the open-chromatin insertion point (±100bp), not the read body. |
| `--call-summits` | Reports the summit position per peak — the anchor for 501bp normalization. Without it you cannot center peaks faithfully. |
| `--keep-dup all` | In single-cell pseudobulk, "duplicate" fragments are genuine independent insertions from different cells, not PCR artifacts. Discarding them throws away real signal. |
| `-q 0.01` | FDR cutoff. MACS3 is intentionally permissive here; downstream pseudo-rep + consensus filtering, not the q-value, is the real noise filter. |

A widely-used alternative ATAC shift is `--shift -75 --extsize 150` (Buenrostro/ENCODE bulk convention; see `iterative-peak-merging`). The `-100 / 200` pairing used here centers a slightly wider Tn5 window; both are valid, just keep one consistent across all groups.

## 501bp summit-centered fixed width

After calling, every peak is re-centered on its summit and resized to a fixed width with `scripts/normalize_width.R::normalize_to_501bp` (summit ± `extend`, default `extend=250` → 501bp).

Why fixed-width 501bp:
- **chromVAR / TFBSTools standard.** Equal-width windows make per-peak accessibility directly comparable and keep motif scanning unbiased by peak width.
- **Summit-centered.** The summit is the Tn5 insertion maximum; centering keeps the strongest base at the window center so overlap resolution (the MERGE stage) compares like with like.
- **Enables iterative overlap.** Uniform widths make the winner-take-all overlap removal well-defined.

## QC filters applied to the raw calls

Run these before / during MERGE (provenance: `S2_6b_atac_peak_recalling_multistrategies.R` and `createIterativeOverlapPeakSet.R`):

1. **Standard-chromosome filter.** Keep only `chr1`-`chr19`, `chrX`, `chrY` for mouse (`chr1`-`chr22`, `chrX`, `chrY` for human). Drop scaffolds / `_` contigs / `chrM`.
   ```r
   std_chroms <- paste0("chr", c(1:19, "X", "Y"))
   peaks <- keepSeqlevels(peaks, std_chroms, pruning.mode = "coarse")
   ```
2. **Chromosome-boundary (cliff) filter.** Drop peaks whose extended window runs past a chromosome end (`subsetByOverlaps(peaks, chrom_sizes, type = "within")`).
3. **ENCODE blacklist removal.** `scripts/blacklist.R::remove_blacklist_peaks` (subsetByOverlaps invert). The source project lifts the ENCODE mm10 blacklist v2 to mm39; `load_blacklist(bed_path=...)` is parameterized for any build.
4. **N-content filter.** Drop peaks whose sequence is mostly unknown bases — the source uses a strict `N <= percent` with `percent = 0` (any N). A relaxed variant drops peaks with N fraction `>= 0.1`. Implement via `Biostrings::letterFrequency` on `getSeq(BSgenome, peaks)`.

## Signac CallPeaks wrapper (unpaired pseudobulk)

The unpaired project calls peaks through Signac rather than shelling out to MACS3 directly:

```r
peaks <- Signac::CallPeaks(
  object = obj,
  group.by = "transferred_celltype",        # set the grouping!
  effective.genome.size = 1.87e9,            # mouse mm39 effective size
  additional.args = "-q 0.01 --call-summits --nolambda --keep-dup all"
)
```

Notes:
- `--nolambda` is added for **sparse pseudobulk**: with few cells the local-lambda background model is unstable, so MACS uses a global background instead.
- `effective.genome.size = 1.87e9` is the mappable mouse genome; use `2.7e9` for human (`-g hs`).
- Always set `group.by` — see the pitfall "MACS returns identical peaks across groups."

Use the direct MACS3 call when you control the pseudobulk/pseudorep splitting yourself (the multiome pipeline); use the Signac wrapper when Signac already manages fragments and grouping (the unpaired pipeline). Both feed the same MERGE stage.

## See also

- `references/iterative-overlap-merge.md` — what happens to these calls next.
- `snapatac2-atac-preprocessing`, `signac-chromatin-analysis` — full calling/fragment mechanics.
