# Test fixtures — SYNTHETIC and DE-IDENTIFIED

Every file here is **hand-authored / synthetic**. None is a copy, sample, or derivative of any
real dataset, and **none can reconstruct the biology of any real dataset**. They exist solely to
exercise code branches and catch regressions in offline (no-API) tests.

## What was deliberately kept out

- **No identity / condition metadata.** No cell types, lineages, donors, treatments,
  temperatures, sample IDs, barcodes, GEO/Synapse accessions, or project paths.
- **Abstract axes and categories.** Evidence axes are generic (`proliferation`, `stress_response`,
  `apoptosis`, `activation`); categories are `catA/catB/catC`; programs are `P1/P2/P3/P9` with
  placeholder genes (`G1`, `GA`, `RPL1`, …). These encode no real program structure.
- **Numbers are chosen to hit branches**, not measured: each cluster targets a specific render path
  (axis cap, companion force-render, artifact-category flag, unmapped program, signature
  de-saturation, present/absent flag token), not a biological state.

## Files

| File | Purpose |
|---|---|
| `markers_synth.csv` | 6 clusters × generic public gene symbols (textbook housekeeping/cell-cycle genes — zero PII). Marker identity is irrelevant offline. |
| `evidence_synth.csv` | 6 clusters of synthetic binned-z panel + program/signature/flag columns for `PanelEvidenceProvider`. |
| `programs_synth.csv` | 3 synthetic program-decode rows (one `Technical` artifact); program `P9` is intentionally absent → `(unmapped)`. |
| `cassette_response.json` | Synthetic `interactive_consensus_annotation`-shaped dict for the no-API end-to-end test; includes a 1-vs-1 split (cluster 3) and labels covering exact/synonym/novel harmonize branches. |
| `axes_synth.csv` | 6 synthetic per-cell two-axis rows covering every `reconcile` branch. |
| `profile_custom.yaml` | Minimal custom (join-enabled, closed-vocab) profile for `load_profile` + `reconcile` CLI tests. |
