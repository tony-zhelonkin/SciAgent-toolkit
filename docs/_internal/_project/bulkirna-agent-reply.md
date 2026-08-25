# Reply owed to the bulkiRNA agent — it is not unblocked, and its premise is wrong

**Date:** 2026-08-24. The agent asked for one instruction to unblock three skill
files. The honest answer has two parts, and the second contradicts its analysis.

## Part 1 — no, the three files are unchanged

Verified on disk at `715feb83`:

| Item | Skill | State |
|---|---|---|
| C5 | `coresh-signature-search` | All four scripts still ship: `coresh_batch.R`, `extract_gene_loadings.R`, `symbols_to_entrez.R`, `validate_coresh_install.R` |
| C6 | `bulk-rnaseq-gsea` | `references/msigdb.md:7` still says "the RNAseq-toolkit wrapper around `clusterProfiler::GSEA()`" |
| C6 | `annotate-bulk-rnaseq-data` | `SKILL.md:71` still reads `RNAseq-toolkit **v0.2.0**` |

**Zero skills mention bulkiRNA.** One part did move: the
`01_scripts/RNAseq-toolkit/...` citations no longer name any single mount. See
part 3.

## Part 2 — the gate is not only a sequencing decision

The agent wrote: *"Phase 3 already discharged the prerequisite — a skill can reach
an installed bulkiRNA in the image."* **Measured, that is false today.**

| Layer | bulkiRNA | `coresh_*` |
|---|---|---|
| `scdock-r-dev:v0.5.10` — **7 of 8 live containers** | **absent entirely** | no |
| `scdock-r-dev:v0.5.13` — 1 live container | 0.4.0, 64 exports | no |
| `scdock-r-dev:v0.5.14` — built, **running nowhere** | 0.6.0, 79 exports | yes |
| source `HEAD` | 1.1.0.9000, 59 exports | yes |

`packageVersion("bulkiRNA")` errors with "there is no package called 'bulkiRNA'" on
`v0.5.10`. So `coresh_search()` reaches **0 of 8** containers. Rewriting C5 to name
it would replace a stale instruction with an unrunnable one.

The owner has since **deferred container recreation** until scio and bulkiRNA are
both parked and agreed. So the instruction the agent asked for cannot be given yet,
and the reason is a delivery fact rather than a preference.

## Part 3 — what the floor must be when it does unblock

The safe set is the **58-function intersection of `0.6.0` and source `HEAD`** —
everything except `gs_master_columns`. Naming anything outside it buys one image
version and then breaks:

- **21 legacy exports** are live in the pinned image and already deleted upstream.
  `run_gsea`, `load_reference_db`, `normalize_gsea_results`,
  `download_gatom_references` and the `gsea_*` plot family are among them. The
  catalog currently names ten. **Do not route to any of them.**
- `gs_master_columns` is upstream-only.

## What did change, and why it is a waypoint

`01_scripts/RNAseq-toolkit/` was recorded by phase 01 as existing in zero projects.
That was wrong: it exists in `12868-EH`, and `Gama_Vivian_DRP1_bulkRNAseq` uses a
third spelling, `01_Scripts/`. Seven projects use `01_modules/`. So no single
spelling is correct for the fleet, and the citations now read `<toolkit>/`, defined
per file as wherever the project vendors the tree. When the package lands these
lines name a function and the vendored path leaves the catalog.

## Two corrections to pass back

1. **The submodule is no longer untouchable.** The standing "do not touch that
   submodule" instruction was superseded: this repository has had 24 commits this
   session, and the owner has authorised fleet writes. Collision risk is now the
   only real coordination concern, which is what the agent already identified.
2. **The repository is being renamed to `scio`.** Its own paths
   (`scbio-docker/toolkits/SciAgent-toolkit/skills/...`) change. Cite skills by
   name, or through `.claude/skills/` and `.agents/skills/`, rather than by vendor
   path — `lint --check harness-links` now reports the vendor path as a read path.

## Task
`#47` phase 02 stays blocked, gated on `#49`, which the owner deferred. `#50`, the
reverse pass for the unclaimed `de_*` family, is gated the same way.
