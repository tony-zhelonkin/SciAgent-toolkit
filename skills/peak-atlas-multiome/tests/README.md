# peak-atlas-multiome — test scaffolds

Non-brittle test scaffolds for the paired-multiome peak-calling logic. They
assert **invariants / properties** on synthetic fixtures and are guarded by
`skip_scaffold()` until wired up in an R env. See
`peak-atlas-framework/tests/README.md` for the full testing philosophy (it is
the shared spine; this skill follows the same rules).

## What is covered

| Script / check | Key invariants tested |
|---|---|
| `apply_primary_rescue_filter.R` | adaptive threshold `= clamp(rate*n, floor, cap)`; PRIMARY keeps a peak clearing its bar in ANY cluster; **rare-type protection** (adaptive keeps a peak a global threshold would drop); RESCUE keeps 4-strategy peaks; `final = primary | rescue`; monotone in `floor` |
| `call_peaks_multistrategy.R` | pseudo-rep split partitions cells ~50/50, deterministic; reproducibility filter keeps only rep1∩rep2; narrowPeak 0→1 based; reconciliation non-overlapping + support set; **group.by required** (mock-MACS contract) |
| `checks/check_marker_retention.R` | retention ≥ 0.95 → exit 0; below, or a panel dropping to zero → exit 1 |

The MACS3 / Signac `CallPeaks` paths are **mocked** (`mock_call_peaks`) so the
pure pseudo-rep / reconciliation logic is testable without external tools; the
gene-coordinate (EnsDb) parts of marker retention are documented as TODO
fixtures.

## Running

```bash
bash tests/run_skill_tests.sh
# Set PEAK_ATLAS_FRAMEWORK_SCRIPTS if the framework lives elsewhere; the helper
# defaults it to ../peak-atlas-framework/scripts so sourced primitives resolve.
# (call_peaks_multistrategy.R resolves the same path relative to its OWN file
# when the env var is unset, so this is belt-and-braces, not a requirement.)
```
