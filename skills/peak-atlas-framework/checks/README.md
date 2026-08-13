# checks/

Standalone R validators for a candidate peak atlas. Each exits non-zero on
failure so it can gate a pipeline.

| Check | Asserts |
|---|---|
| `check_peak_atlas.R` | No overlapping peaks; every peak exactly the fixed width (default 501bp); zero blacklist overlaps (when supplied); valid coordinates. Reproduces `validate_peak_atlas()`. |
| `check_frip_retention.R` | Filtered median FRiP >= threshold (default 0.90) of the original median. |

Usage:

```bash
Rscript check_peak_atlas.R atlas.rds [blacklist.bed] [expected_width]
Rscript check_frip_retention.R original_frip.rds filtered_frip.rds [min_ratio]
```

See `references/validation-battery.md` for the full validation battery and the
suggested gate order.
