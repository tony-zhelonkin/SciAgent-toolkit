# Phase 01 — the three-way inventory

**Repo:** scio · **Read-only. Produce a table, change no skill.** · Read
`00_INDEX.md` first.

## Why read-only

Two of the three versions in play disagree about whether the functions exist. A
rewrite authored against the source tree would name an API that no live container
has. Establish the ground truth first; phase 02 acts on it.

## What to produce

One row per **claim**, not per skill. A claim is any place a skill names a script
path, a `source()` call, a vendored-toolkit version, or a package function.

| Column | Meaning |
|---|---|
| `skill` | directory name |
| `file:line` | where the claim is |
| `claim` | the path, function or version named, verbatim |
| `source_api` | the export in the bulkiRNA source tree that covers it, or `—` |
| `installed_api` | the same, checked against what a container actually has |
| `verdict` | `superseded` · `superseded-but-undelivered` · `stands` · `unrelated` |

`superseded-but-undelivered` is the row that matters: the package absorbed it,
and the reader's container does not have the absorbing function yet.

## How to check each layer

**Source tree** — `/data1/users/antonz/pipeline/bulkiRNA`, `1.1.0.9000` at
`ca99d8a`:

```bash
grep '^export' /data1/users/antonz/pipeline/bulkiRNA/NAMESPACE | sed 's/export(//;s/)//' | sort
```

**Installed** — R is **not on the host**. Use a running container, and check more
than one image version, because the fleet spans two:

```bash
docker ps --format '{{.Names}}\t{{.Image}}' | grep dev-core     # v0.5.10 and v0.5.13 are both live
docker exec <container> bash -lc 'Rscript -e "
  cat(as.character(packageVersion(\"bulkiRNA\")), \"\n\");
  cat(sort(getNamespaceExports(\"bulkiRNA\")), sep=\"\n\")"'
```

Measured 2026-08-21 in `jr-mc_devcontainer-dev-core-1` (`scdock-r-dev:v0.5.13`):
**bulkiRNA 0.4.0, 64 exports, no `coresh_*` layer and no `bulkirna_api`.** Check a
`v0.5.10` container too — if it differs, the inventory needs a column per image
version and phase 02 needs a floor.

Do **not** trust `bulkirna_api()`'s `superseded_by` column: `R/api.R:130` fills it
with `NA`, and it is absent from 0.4.0 entirely.

**The claims:**

```bash
cd <toolkit>/skills
grep -rn "RNAseq-toolkit\|clusterProfiler\|01_scripts/\|source(" \
  annotate-bulk-rnaseq-data bulk-rnaseq-gsea bulk-rnaseq-pathway-explorer \
  coresh-signature-search gatom-metabolomic-predictions iterative-peak-merging \
  nfcore-rnaseq-execution peak-atlas-framework peak-atlas-multiome \
  star-te-preprocessing te-gene-featurecounts te-geneset-gsea \
  te-reference-saf-build
```

## Three questions the table must answer

1. **Which of the thirteen are actually in scope?** The three ATAC skills
   (`iterative-peak-merging`, `peak-atlas-framework`, `peak-atlas-multiome`)
   probably name `01_scripts/` for an unrelated reason. Say so and drop them
   rather than carrying them forward.
2. **What does the TE track need?** Five skills belong to `TE-RNAseq-toolkit`, a
   different package. Record whether it has an exported API at all and stop
   there — its audit is not this plan.
3. **What is the floor?** Phase 02 has to pick one: write against the source API
   and accept that the text is ahead of every container until the image pin
   moves, or write against the installed API and revisit. Recommend one, with the
   version numbers behind it.

## Also record, in the other direction

Fifteen exports no skill mentions — `coresh_chunks`, `coresh_match`,
`coresh_sets`, `coresh_convergence`, `gs_coregulation`, `gs_leading_edge`,
`gs_ranks`, `gs_score`, `de_volcano`, `de_volcano_grid`, `de_pca`, `de_pca_3d`,
`de_md_plot`, `de_bfc_plot`, `filter_confounder_genes`. For each: is this
capability an existing skill should route to, a skill that does not exist yet, or
correctly invisible? One line each.

## Deliverable

The table, the three answers, the reverse pass, and a one-paragraph
recommendation on the floor. Write it to `01_inventory.report.md` beside this
brief. Change nothing under `skills/`.
