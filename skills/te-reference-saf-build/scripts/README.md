# scripts/

Executable helpers the skill invokes directly. Self-documenting (usage header at the top).

## Contents

- **`build_te_saf.sh`** — the OWNED, single-source-of-truth build script for the TE reference SAFs. Inputs: `--te-gtf` (TEtranscripts RepeatMasker GTF), `--gene-gtf` (gene GTF for the exon BED, same build), `--out-dir`, `--prefix`. Produces `<prefix>_GROUPED_all.saf` (grouped, `GeneID = Subfamily:Family:Class`) and `<prefix>_GROUPED_all_noExon.saf` (exon-subtracted, 0-overlap-verified). Derives the exon BED FROM `--gene-gtf` (never a hardcoded sibling-project path), normalizes contig names, deterministic. `--keep-classes REGEX` optionally whitelists classes (default: keep all). Deps: `gawk`, `bedtools`, coreutils. Smoke-tested by `tests/run_skill_tests.sh`.
