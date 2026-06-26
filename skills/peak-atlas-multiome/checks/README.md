# checks/

Runnable gates that exit non-zero on failure. The structural and FRiP gates are
shared and live in `peak-atlas-framework/checks/` (`check_peak_atlas.R`,
`check_frip_retention.R`); this child adds only the multiome-specific marker gate.

| Check | Asserts |
|---|---|
| `check_marker_retention.R` | Rare-marker promoter-peak retention of the filtered atlas vs the original is `>= 0.95` (default), and no marker set loses ALL its promoter peaks. Default panels: HSC (Gata2/Procr/Hlf), pDC (Tcf4/Siglech/Bst2), Neutrophil (S100a8/S100a9/Cebpe). Override with a markers file; set `MARKER_ENSDB` for non-mouse builds. |

Usage: `Rscript checks/check_marker_retention.R original.rds filtered.rds [min_ratio] [markers.txt]`
