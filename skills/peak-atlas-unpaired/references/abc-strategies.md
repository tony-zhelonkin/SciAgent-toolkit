# Strategies A / B / C in depth

The crux of the unpaired atlas is calling peaks three complementary ways over
the SAME ATAC cells, then letting the framework's support-voting and merge
combine them. Each strategy answers a different question and has a different
robustness to label-transfer error. Provenance: `S5.1b` (A), `S5.2b` (B),
`S4.1c` (C).

## The three groupings

| Strategy | Grouping column | Min cells | Question it answers | Transfer-error robustness |
|---|---|---|---|---|
| **A** | coarse cell type (`r1_highconf`) | 100 | What are the lineage-identity peaks (cDC1 vs cDC2 vs Mac vs Mono)? | partial — coarse labels are stable |
| **B** | refined subtype × condition (`r2_highconf` + `orig.ident`) | **20** | What are the rare condition-specific interaction peaks? | low — fine labels are noisiest |
| **C** | condition only (`orig.ident`) — **label-free** | 100 | What is the accessible landscape per condition, at maximum power? | full — uses NO transferred labels |

All three call peaks the same way (the shared Signac wrapper) and post-process
the same way (501bp + chromosome/blacklist/N filters + `reduce(gapwidth=50)`).
Only the grouping column and minimum-cell threshold differ.

## Strategy A — lineage identity (coarse)

Groups by the R1 coarse label after dropping LowConf cells. With `min = 100`
cells/group it is statistically solid. Each group is peak-called, post-processed
to 501bp, and the groups are merged; the merged peak's **score is its support
count** (`n_celltypes`) — how many of the coarse types contain it:

- support 4 (all types) = shared / core accessibility,
- support 1 = lineage-specific peak.

```r
atac_A <- subset(atac, r1_highconf != "LowConf")
A <- call_strategy_peaks(atac_A, "r1_highconf", bl, bsg,
                         min_cells = 100, strategy = "A")
```

## Strategy B — refined subtype × condition (the rare-biology lens)

Groups by `paste0(r2_highconf, "_", orig.ident)` after dropping LowConf cells.
This is where condition-specific rare biology lives — e.g. a subtype that only
appears in a knockout condition. The **minimum is lowered to 20** cells
(from the usual 50) precisely to keep such rare groups:

> a biologically critical subtype with ~43 cells in one condition and ~24 in
> another would be erased by a 50-cell floor; 20 keeps it.

```r
atac_B <- subset(atac, r2_highconf != "LowConf")
atac_B$cluster_condition <- paste0(atac_B$r2_highconf, "_", atac_B$orig.ident)
B <- call_strategy_peaks(atac_B, "cluster_condition", bl, bsg,
                         min_cells = 20, strategy = "B")
```

B's support count (`n_groups`) and the derived `n_conditions` feed several tier
rules (`references/tier-stratification.md`), so B is the strategy that most
shapes the final atlas's rare-biology content. It is also the least robust to
transfer error — which is exactly why it is hedged by C.

## Strategy C — label-free condition pseudobulk (the hedge)

Groups by condition (`orig.ident`) ONLY. No cell-type labels are used, so
**Strategy C is immune to label-transfer error**: even if the RNA→ATAC transfer
were entirely wrong, C still recovers the per-condition accessible landscape at
maximum statistical power (every cell in the condition contributes). It is the
safety net that catches real peaks A and B might miss because of mislabeling.

```r
C <- call_strategy_peaks(atac, "orig.ident", bl, bsg,
                         min_cells = 100, strategy = "C")   # LowConf cells included
```

C is also the only strategy that legitimately includes LowConf cells, because it
never relies on their labels.

## The group.by bug (why it matters)

`Signac::CallPeaks(idents = g)` alone does NOT subset fragments — it only
subsets cells for cell selection, while MACS is still handed every fragment
file. The consequence is that **every group returns identical peaks**. Setting
`group.by` triggers Signac's fragment filtering so MACS sees only the group's
fragments:

```r
CallPeaks(object = atac, group.by = group_col, idents = g, ...)   # correct
CallPeaks(object = atac, idents = g, ...)                          # BUG: identical peaks
```

This was a real bug fixed in the source pipeline; it is also the framework
pitfall "MACS returns identical peaks across groups."

## What each strategy contributes (read it off the UpSet)

Plot an UpSet / 3-way Venn of A vs B vs C peaks:
- the **A∩B∩C core** is high-confidence shared accessibility,
- **C-only** peaks are the transfer-error hedge paying off (real signal A/B
  missed),
- **B-only** peaks are candidate rare biology — promising but the riskiest, to
  be gated by the external-correlation test before trusting (Tier 0c).

## See also

- `references/label-transfer-bridge.md` — where A's and B's labels come from.
- `references/tier-stratification.md` — how support counts become tiers.
- `peak-atlas-framework/references/support-voting.md` — `n_strategies` and the
  compute-before-merge rule the framework applies across A/B/C.
- `scripts/label_transfer_strategy.R` — `call_strategy_peaks`,
  `process_group_peaks`.
