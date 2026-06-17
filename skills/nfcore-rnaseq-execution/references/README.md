# references/

Deep-dive topic guides loaded on demand from `SKILL.md`. One file per topic (e.g. `parameters.md`, `advanced-usage.md`, `benchmarks.md`). Aim for ≤300 lines per file; anything longer should be split by sub-topic. Referenced from SKILL.md as `see references/<topic>.md`.

## Contents

- **`nfcore-options.md`** — nf-core/rnaseq options cheatsheet: launch CLIs (modern `nf-core pipelines launch` + classic `nextflow run`), samplesheet columns (incl. `seq_platform`/`seq_center` and multi-row read concatenation), `strandedness=auto` inference + MultiQC verification, fq lint behaviour, two-step `--save_align_intermeds` → `--skip_alignment` BAM reprocessing, and optional rRNA removal (incl. the GPU `ribodetector --use_gpu_ribodetector` path). Loaded on demand from SKILL.md.
- **`dataset-record-template.md`** — the reusable per-dataset provenance template. Every nf-core/rnaseq run instantiates a filled copy (project ID, pipeline `-r` + revision hash, verbatim `nextflow run` command, `te_star.config` snapshot, samplesheet path, VERIFIED strandedness, SAF/featureCounts invocation, software versions, output paths, QC, per-dataset deviations). First line points back to the canonical recipe (`star-te-preprocessing` + `nfcore-rnaseq-execution`). Durable home: beside results on `/data1` (source of truth) + Obsidian vignette mirror (ADR-E = Both).

> Note: the canonical TE STAR args and `te_star.config` are **owned by `star-te-preprocessing`**, not duplicated here. This skill points to them.
