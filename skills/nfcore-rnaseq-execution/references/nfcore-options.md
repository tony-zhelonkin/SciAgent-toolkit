# nf-core/rnaseq options cheatsheet

Operational options for `nf-core/rnaseq` (current latest `-r 3.26.0`; prior 13036/AdaW runs used `-r 3.20.0`). Always pin `-r` explicitly — changing the version can change STAR/tool defaults, so the TE recipe in `star-te-preprocessing` must be re-validated against whichever version a given run pins (decided at Phase V, not here).

Authoritative docs: https://nf-co.re/rnaseq/3.26.0

## Launch CLIs (both valid)

```bash
# Modern interactive launcher (prompts/validates params, writes nf-params.json):
nf-core pipelines launch nf-core/rnaseq -r 3.26.0

# Classic direct invocation (still fully supported; used by the dataset records):
nextflow run nf-core/rnaseq -r 3.26.0 -profile docker --input ... --outdir ...
```

## Samplesheet columns

Required (4): `sample,fastq_1,fastq_2,strandedness`.

Optional:
- `seq_platform`, `seq_center` — written into BAM read-group `PL`/`CN` tags. For a uniform value across all samples, use params `--seq_platform ILLUMINA` / `--seq_center <name>` instead of per-row columns.
- Reprocessing workflow only: `genome_bam`, `transcriptome_bam`, `percent_mapped`.

Behaviour:
- **Same `sample` value on multiple rows → the pipeline CONCATENATES the raw reads** before analysis. This is exactly what `make_samplesheet.sh` relies on: it emits one row per lane sharing the sample id, so multi-lane libraries are merged automatically.
- Spaces in a `sample` value are auto-converted to underscores.

## Strandedness = auto (the default the generator emits)

`strandedness: auto` makes the pipeline subsample 1M reads, infer strand with Salmon, and propagate the result. Thresholds: `--stranded_threshold 0.8`, `--unstranded_threshold 0.1`.

Verification (this is THE mechanism — do not "assume"): the inferred call appears in MultiQC under the **"Strandedness checks"** section, cross-checking Salmon vs RSeQC with pass/fail. Capture the inferred value from MultiQC into the dataset record; the downstream featureCounts gene `-s` must match it.

## fq lint

Runs by default at the start and after every FASTQ-manipulating step; a lint error stops the workflow. Tune with `--extra_fqlint_args`. The paired-read-name check (`P001`) is prone to false failure and is **disabled by default** via `--disable-validator P001`. If a specific validator misfires on usable FASTQs, extend `--extra_fqlint_args` to disable it.

## Two-step BAM reprocessing (cheap quantification re-runs — useful for TE)

TE counting is a post-alignment step, so re-aligning to tweak quantification is wasteful. Instead:

1. Run once with `--save_align_intermeds`. This publishes the BAMs and writes `<outdir>/samplesheets/samplesheet_with_bams.csv`.
2. Reprocess quantification/QC without re-aligning:
   ```bash
   nextflow run nf-core/rnaseq -r <same -r> -profile docker \
     --input <outdir>/samplesheets/samplesheet_with_bams.csv \
     --skip_alignment ...
   ```

Caveats: only works for BAMs produced by THIS pipeline, and you must not mix aligner types — `star_salmon` BAMs reprocess as `star_salmon`.

## rRNA removal (optional, OFF by default)

Enable with `--remove_ribo_rna --ribo_removal_tool {sortmerna|ribodetector|...}`. SortMeRNA default DB is `smr_v4.3_fast_db`.

GPU-accelerated path (this machine has a GPU):
```bash
--remove_ribo_rna --ribo_removal_tool ribodetector --use_gpu_ribodetector
```
`--use_gpu_ribodetector` auto-applies `--gpus all` for Docker; requires the NVIDIA driver and x86_64.

Comparability note: the prior TE runs did **not** remove rRNA. Keep it off for cross-dataset comparability unless you deliberately change it (and record the deviation).
