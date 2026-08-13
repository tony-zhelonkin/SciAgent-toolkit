# references/

Deep-dive topic guides loaded on demand from `SKILL.md`. One file per topic. Aim for ≤300 lines per file. Referenced from SKILL.md as `see references/<topic>.md`.

## Contents

- **`reference-build-record-template.md`** — the reusable per-genome-build provenance template. Every TE reference (re)build instantiates a filled copy: genome build, TE-GTF source URL + date + md5, gene GTF used for the exon BED, exact `build_te_saf.sh` command, output paths + md5s, class-composition choice, verification (incl. the mandatory 0 residual-exon-overlap check), and the datasets that consume the build. Durable home: beside the artifacts under `/data1/shared/ref/<species>/<provider>/<build>/`.

> The TE-counting *usage* runbook and the per-dataset run record are owned by `star-te-preprocessing` / `nfcore-rnaseq-execution`, not duplicated here. This skill owns the build-once construction.
