# TE counting workflow → moved to the `te-gene-featurecounts` skill

The post-nf-core **TE + gene featureCounts counting workflow**, together with its own
locked, version-pinned container, now lives in the **`te-gene-featurecounts`** packaged
skill. That skill is the self-contained, runnable artifact: the pinned image
**`te-fc:2.0.2`** (featureCounts **v2.0.2**), the vendored two-pass driver
(`runFeatureCounts_TE_and_genes.sh` → `runFeatureCounts.sh`), the parameterized Docker
wrapper that does the `bam_fin` symlink-staging + identical-path bind-mount, and the QC
gate the matrices must pass.

See:
- `skills/te-gene-featurecounts/SKILL.md` — the env-locked contract, runbook, QC gate.
- `skills/te-gene-featurecounts/references/te-counting-workflow.md` — the detailed end-to-end runbook.

`star-te-preprocessing` (this skill) still owns the alignment + counting **contract** and the
canonical STAR "Random-One" string; **SAF construction** is owned by `te-reference-saf-build`.
The handoff chain is `te-reference-saf-build` + `star-te-preprocessing` → `te-gene-featurecounts`
→ `annotate-bulk-rnaseq-data`.
