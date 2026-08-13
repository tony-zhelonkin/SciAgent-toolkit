# references/

Deep-dive and canonical artifacts loaded on demand from `SKILL.md`.

Contents:
- **`te_star.config`** — the owned canonical per-project nf-core/rnaseq Groovy config (below).
- **`te-counting-workflow.md`** — the end-to-end TE-counting runbook (post-nf-core): STAR BAM
  staging → TEtranscripts GTF → grouped `_noExon` SAF (`bedtools subtract`) → two-pass
  `runFeatureCounts_TE_and_genes.sh` → gene + TE subfamily matrices → QC gate → handoff to
  `annotate-bulk-rnaseq-data`. Generalized from the per-dataset `/scratch/nf-core/13036-DM/README.md`.
  Points to (does not copy) the TE-RNAseq-toolkit v2.0.0 driver.

## `te_star.config`

The **single owned canonical copy** of the per-project nf-core/rnaseq Groovy config
for TE-compatible runs. This skill is the single source of truth for this file
(SSoT policy, artifact `04-knowledge-management-proposal.md` §3): every other
location (per-dataset run dirs) holds a *frozen snapshot* labelled
"from star-te-preprocessing v1.0.0", never an editable copy.

**Provenance:** copied byte-for-byte from the 13036-DM / AdaW_eWAT_WL canonical
runs, 2026-06-08. It carries references (mm39 / GENCODE vM37), BAM retention
(`save_align_intermeds = true`), the Docker `NXF_UID/GID` trick, `/data2` temp
redirection, and resource caps. **The TE-specific STAR behavior lives in the
`--extra_star_align_args` string in `SKILL.md`, NOT in this config.**
