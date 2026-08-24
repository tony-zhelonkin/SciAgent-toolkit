# Phase 01 — the three-way inventory, measured

**Run 2026-08-21 on the host, read-only. No file under `skills/` was changed.**

> **Correction, 2026-08-24.** Row 2 of §3 says `01_scripts/RNAseq-toolkit/`
> "exists in **zero** projects" and that "every project uses `01_modules/`".
> Both are wrong. A survey of all nine projects holding the tree found **three
> spellings**: `01_modules/` in 7, `01_scripts/` in `12868-EH` — with all six
> cited paths present — and `01_Scripts/` in
> `GVDRP1_prj/Gama_Vivian_DRP1_bulkRNAseq`. The verdict `false today` is wrong
> for row 2; the defect is a **layout migration in progress**, so no single
> spelling is correct for the fleet. Rows 5 and 6 stand as measured.
>
> Phase 02a acted on the bad measurement first (`fff5ff8`) and was corrected in
> `91c3869`: the citations now name no mount at all.

The brief asked for three layers. There are **four**, and the fleet is split down
the middle between two of them. That split, not staleness, is the finding.

---

## 1. The four API layers

| Layer | Version | Exports | `coresh_*` | `bulkirna_api` |
|---|---|---|---|---|
| `scdock-r-dev:v0.5.10` — **7 of 8 live containers** | **absent** | **0** | no | no |
| `scdock-r-dev:v0.5.13` — 1 live container (`jr-mc`) | 0.4.0 | 64 | no | no |
| `scdock-r-dev:v0.5.14` — **built 18 h ago, running nowhere** | 0.6.0 | 79 | **yes** | **yes** |
| source `HEAD` `ca99d8a` | 1.1.0.9000 | 59 | yes | yes |

Measured by `getNamespaceExports("bulkiRNA")` inside each container, and by
`git show <sha>:NAMESPACE` for the pinned commit.

Two corrections to `00_INDEX.md`, which was written from the source tree and one
container:

- **bulkiRNA is not installed at all on `v0.5.10`.** `packageVersion` errors with
  "there is no package called 'bulkiRNA'". The INDEX said "eight of nine
  containers run v0.5.10 while VERSION says v0.5.14" and inferred a version gap;
  the gap is an absence.
- **`BULKIRNA_SHA = e42c2de1` is tag `0.6.0`, not 0.4.0 and not HEAD.** It has 79
  exports — more than either neighbour. The `v0.5.14` image that carries it is
  already built.

### The surface did not just grow, it was refactored

```
0.4.0  (64)  ⊂  0.6.0  (79)          0.6.0 = 0.4.0 + 15 new
0.6.0  (79)  −  21 legacy  + 1  =  HEAD (59)
```

- **0.6.0 adds 15**: the whole `coresh_*` layer, `bulkirna_api`,
  `gene_to_entrez`/`entrez_to_gene`, `gatom_download_refs`, `gsdb_coresh`,
  `gs_coregulation`, `filter_confounder_genes`, `bulkirna_stochastic`.
- **HEAD removes 21 legacy names** that 0.4.0 and 0.6.0 both export:
  `run_gsea`, `run_gsea_analysis`, `load_reference_db`, `list_reference_dbs`,
  `normalize_gsea_results`, `gsea_barplot`, `gsea_dotplot`, `gsea_dotplot_facet`,
  `gsea_running_sum_plot`, `plot_all_gsea_results`, `save_gsea_log`,
  `empty_gsea_tibble`, `list_to_term2gene`, `parse_gmx`, `parse_mitoxplorer`,
  `download_gatom_references`, `convert_human_to_mouse`, `create_MD_plot`,
  `create_standard_volcano`, `custom_minimal_theme_with_grid`, `filter_by_size`.
- **HEAD adds 1**: `gs_master_columns`.

**This is the trap that decides phase 02.** The catalog already names **ten of
those twenty-one** legacy functions:

| Legacy name the catalog uses | Files | In 0.4.0 / 0.6.0 | In HEAD |
|---|---|---|---|
| `normalize_gsea_results` | 4 | yes | **gone** |
| `load_reference_db` | 3 | yes | **gone** |
| `run_gsea` | 3 | yes | **gone** |
| `download_gatom_references` | 2 | yes | **gone** |
| `gsea_barplot`, `gsea_dotplot`, `gsea_running_sum_plot` | 2 each | yes | **gone** |
| `run_gsea_analysis`, `parse_gmx`, `parse_mitoxplorer` | 1 each | yes | **gone** |

The naive modernisation — prefix the existing names with `bulkiRNA::` — compiles
against the pinned image and is **already wrong against the package's own HEAD**.
`gatom-metabolomic-predictions:101` is the clean example: it sources
`load_reference_db.R`, the tempting fix is `bulkiRNA::download_gatom_references()`,
and the correct target is `gatom_download_refs`, which does not exist in 0.4.0.

---

## 2. The fleet is split, and each half breaks the other half's instruction

| Project | Image | bulkiRNA | `01_modules/RNAseq-toolkit` |
|---|---|---|---|
| Meta-Aging | v0.5.10 | absent | **present** |
| Meta-Aging/14616-DM | v0.5.10 | absent | **present** |
| Meta-Aging/14782-DM | v0.5.10 | absent | **present** |
| 14839-DM-cGAS | v0.5.10 | absent | **present** (+ TE toolkit) |
| DC-nexus | v0.5.10 | absent | **present** |
| DC-nexus/DC_mouse_cancer | v0.5.10 | absent | **present** |
| STING-JR | v0.5.10 | absent | **present** |
| **JR-MC** | **v0.5.13** | **0.4.0** | **absent** |

The two conditions are perfectly complementary. Vendored tree present ⇔ package
absent.

- The skills' current `source("01_modules/RNAseq-toolkit/…")` instructions are
  **true and executable in 7 of 8 projects**, and dead in JR-MC.
- A bulkiRNA rewrite would be **executable in 1 of 8**, and only partly — 0.4.0
  has no `coresh_*` layer at all.

The vendored tree is `v0.2.0-9-g752481f`, and every cited file is really there
(`scripts/General/`, `scripts/GSEA/GSEA_processing/run_gsea.R`). So the skills'
"RNAseq-toolkit v0.2.0" claim is **accurate**, not stale. That is the opposite of
what the origin report assumed.

**Consequence: no single version of the text is correct for the fleet as it
stands.** This is not a documentation-lag problem that a rewrite fixes. It is a
half-finished migration, and the skill is being asked to paper over it.

---

## 3. The claim inventory

141 grep hits across the thirteen skills. Grouped by claim, because the same
claim repeats:

| # | Claim | Where | Verdict |
|---|---|---|---|
| 1 | `source("01_modules/RNAseq-toolkit/scripts/General/{io_helpers,annotate_genes,dge_helpers,provenance}.R")` | `annotate-bulk-rnaseq-data` SKILL:71, gene-annotation:22,64–68, te-annotation:49,70–71 | **superseded-but-undelivered** — `read_counts_matrix`, `read_metadata`, `annotate_genes`, `build_dge`, `write_session_provenance`, `ensure_dir` exist in all three package versions; delivered to 1 of 8 |
| 2 | `01_scripts/RNAseq-toolkit/…` (10 hits, `bulk-rnaseq-gsea` ×5 files) | custom-db:96, master-tables:405,407, msigdb:318,319,399, visualization:435–437, SKILL:196 | **false today, independent of bulkiRNA** — every project uses `01_modules/`. This path exists in **zero** projects. Two spellings of the same mount in one skill family |
| 3 | `source(file.path(DIR_TOOLKIT, "GSEA/GSEA_processing/run_gsea.R"))` | `bulk-rnaseq-gsea` SKILL:77, msigdb:51,328 | **superseded-but-undelivered** → `gs_test`. Note `run_gsea` also survives as a *legacy export* in 0.4.0/0.6.0 — do not route to it |
| 4 | `clusterProfiler::GSEA(…)`, `library(clusterProfiler)` (34 hits) | `bulk-rnaseq-gsea` ×4 files, `coresh-…-bridge:125`, `bulk-rnaseq-pathway-explorer:20` | **stands** — `clusterProfiler` and `fgsea` are installed in v0.5.10 *and* v0.5.13. `gs_test` wraps it; the trap lines (`eps = 0`, `p.adjust` vs `qvalue`) stay regardless |
| 5 | `source("02_analysis/helpers/normalize_gsea.R")`, `…/pathway_utils.R` (6 hits) | master-tables:56,104–105, custom-db:107 | **false today** — `02_analysis/helpers/` is empty in 7 of 8 projects; the eighth holds six different files. Neither file exists anywhere |
| 6 | `source("02_analysis/config/config.R")` (6 hits) | `bulk-rnaseq-gsea` SKILL:76, msigdb:50, master-tables:104, visualization:61–62 | **false in 7 of 8** — present only in 14616-DM |
| 7 | `load_reference_db()` + `download_gatom_references()` | `gatom-metabolomic-predictions:99,101` | **superseded-but-undelivered, and the obvious target is a trap** → `gatom_download_refs`/`gatom_refs`. `gatom` the package is also **absent from v0.5.10**, so this skill is doubly gated |
| 8 | `scripts/{coresh_batch,extract_gene_loadings,symbols_to_entrez,validate_coresh_install}.R` | `coresh-signature-search`, self-contained, relative paths | **superseded-but-undelivered** → `coresh_search`, `coresh_loadings`/`gsdb_coresh`, `gene_to_entrez`, `coresh_validate`. All four are **0.6.0-only** — delivered to **0 of 8** containers |
| 9 | `01_modules/TE-RNAseq-toolkit/R/*.R` + "**v0.1.0**" (34 `TE-RNAseq-toolkit` hits) | `annotate-bulk-rnaseq-data`, `te-geneset-gsea` | **stands, version label two majors stale** — see §5 |
| 10 | "TE-RNAseq-toolkit **v2.0.0** / **v2.0.1**" | `star-te-preprocessing:68,156`, `te-reference-saf-build:224` | **stands** — actual is `v2.0.3-5`. Same repo labelled `v0.1.0` and `v2.0.x` inside one catalog |
| 11 | `<TE-RNAseq-toolkit / bulkRNAseq_pipeline_scripts vX.Y.Z>` | `nfcore-rnaseq-execution/references/dataset-record-template:79` | **stands** — a fill-in placeholder, not a claim |
| 12 | `01_scripts/R_scripts/createIterativeOverlapPeakSet.R` (5 hits) | `iterative-peak-merging` | **unrelated** — a consumer analysis script. Belongs to #8 (`scripts/`→`stages/`) |
| 13 | `01_scripts/R/peak_utils.R` provenance comments; intra-skill `source()`/`sys.source()` | `peak-atlas-framework`, `peak-atlas-multiome` | **unrelated** — provenance to `JBader_scHFD`, and self-referential sourcing |

The brief's four verdicts needed a fifth. **`false today`** — rows 2, 5, 6 — is a
claim that names a path existing in no project, and it has nothing to do with
bulkiRNA. Twenty-two hits. It is the same defect class as everything else in the
handoff's table: *an instruction naming something the mechanism does not
guarantee.* It is also the cheapest thing in this report to fix, and fixing it
does not wait on any image decision.

---

## 4. Q1 — the ATAC trio is out

**Yes, drop all three.** `iterative-peak-merging`, `peak-atlas-framework`,
`peak-atlas-multiome` matched the grep for two innocent reasons: provenance
comments citing `/data2/users/JCRLab/JBader/JBader_scHFD/01_scripts/R/peak_utils.R`
(the project the skill was distilled from), and `source()`/`sys.source()` calls
that load a *sibling file inside the same skill*. No RNA package, no vendored
toolkit. Thirteen skills in scope becomes ten.

## 5. Q2 — the TE track needs a version reconciliation, not an API migration

`TE-RNAseq-toolkit` has **no `NAMESPACE` and no `DESCRIPTION`** — it is not an R
package, so it exports nothing and `source()` is the only route available. No
audit against an API is possible.

What it does have is a version problem the catalog owns:

- Source and the one project vendoring it are both at **`v2.0.3-5-gb78935f`**.
- `annotate-bulk-rnaseq-data` and `te-geneset-gsea` call it **`v0.1.0`**.
- `star-te-preprocessing` and `te-reference-saf-build` call it **`v2.0.0`/`v2.0.1`**.

I checked whether the `v0.1.0` claims survived the jump. **They do.** All four
cited files are present at `v2.0.3`; `create_te_genesets.R:109,175` still carries
the hardcoded `source("01_modules/TE-RNAseq-toolkit/R/te_utils.R")` that
`te-annotation.md`'s T1b note documents; the `gs_name`/`gene_symbol` column
contract still holds (`create_te_genesets.R:85–87`). So the substance is right and
only the label is wrong — a label fix, not a rewrite. Recorded and stopped, per
the brief.

## 6. Q3 — the floor

**Recommendation: the floor is the 58-function intersection of `0.6.0` and source
`HEAD`, and phase 02 does not start until the fleet is on `v0.5.14`.**

Two parts, and the second is the load-bearing one.

**The set.** `0.6.0 ∩ HEAD` is 58 of 59 current exports — everything but
`gs_master_columns`. A skill written against that set is correct in the pinned
image *and* correct against the package as it stands today, so it survives drift
in both directions. Naming anything outside it costs one image version of
correctness and then breaks: the 21 legacy exports are live in the pin and
already deleted upstream, and `gs_master_columns` is upstream-only. Ten of those
21 are in the catalog right now (§1), so this is not hypothetical.

**The gate.** Writing against any bulkiRNA version is wrong for 7 of 8 projects
today, because the package is not installed there. The image that fixes this is
built and running nowhere. So the ordering is not "rewrite, then hope the image
follows" — it is:

```
recreate containers on v0.5.14  ──▶  the floor becomes real  ──▶  phase 02
```

I considered the alternative — have the skill probe
(`requireNamespace("bulkiRNA")`) and document both routes. It is the house rule
("probe, don't assert") and it is what I would do if the migration were going to
sit half-finished for weeks. I am recommending against it as the *default*
because a two-route skill doubles the text in the one place text is scarce, and
it teaches the vendored `source()` pattern to every reader for as long as it
stands. The probe is the right answer only if §7's decision comes back "the fleet
stays split".

**What phase 02 can do before the gate, safely:** rows 2, 5, 6 and 10 — the 22
`false today` hits and the TE version labels. Those are wrong under every image
version and their fix does not name a package function.

---

## 7. The decision this hands back

Phase 02 as written in `02_pilot.md` assumed a rewrite was executable. It is not,
yet. The gate is one owner decision:

**Recreate the eight dev-core containers on `scdock-r-dev:v0.5.14`?** The image
is built. Doing it makes the 0.6.0 API real everywhere, retires the vendored
`01_modules/RNAseq-toolkit` submodule from seven projects, and unblocks phase 02
for all five bulkiRNA-scope skills at once. Not doing it leaves the fleet split
and phase 02 has to ship the probe-and-branch form instead.

This is the two-hop propagation problem in a third place. Hop 1 is not the
toolkit submodule here, it is the **image tag** — and it has the same property
the trace measured: a correct artifact exists and no mechanism delivers it.

---

## 8. The reverse pass — capability the catalog never learned

Fifteen 0.6.0-and-HEAD exports that no skill names. Confirmed by
`grep -rlw` across `skills/`: **the catalog contains zero occurrences of the
string "bulkiRNA"**, and names only six floor functions at all
(`annotate_genes`, `read_counts_matrix`, `read_metadata`, `build_dge`,
`write_session_provenance`, `format_pathway_name`) — and those only because
RNAseq-toolkit's scripts happened to use the same names.

| Export | Verdict |
|---|---|
| `coresh_chunks`, `coresh_match`, `coresh_sets` | plumbing under `coresh_search` — correctly invisible |
| `coresh_convergence` | a **missing trap doc**: "did the sweep converge" is exactly the judgement `coresh-signature-search` should teach. Route to it |
| `gs_coregulation` | new capability, no skill. Candidate for a skill, not a paragraph |
| `gs_leading_edge`, `gs_ranks`, `gs_score` | belong in `bulk-rnaseq-gsea` — leading-edge extraction is a routine question it currently cannot answer |
| `de_volcano`, `de_volcano_grid`, `de_md_plot`, `de_bfc_plot` | no DE-visualisation skill exists. The gap the catalog is least aware of |
| `de_pca`, `de_pca_3d` | same; QC-adjacent |
| `filter_confounder_genes` | a methodology decision with a wrong default — wants a decision gate, not a wrapper mention |
| `bulkirna_stochastic` | reproducibility contract; belongs in whatever skill sets seeds |

**Two skills' worth of uncovered capability, and one missing trap.** Worth its
own plan; do not fold it into phase 02.

## 9. Also confirmed

- `bulkirna_api()`'s `superseded_by` is `rep(NA_character_, …)` at `R/api.R:130`,
  and `bulkirna_api` is absent from 0.4.0 entirely. The mapping in §3 was
  authored by reading and will have to be re-authored next time. Recommendation 1
  in `00_INDEX.md` §6 stands and is now better evidenced: **populating that
  column turns this report into a join.**
- The vendored `RNAseq-toolkit/MIGRATION.md` is about git branching strategy, not
  the bulkiRNA consolidation. There is no migration guide.
- `/data1/users/antonz/pipeline/RNAseq-toolkit/` is empty. The source of truth for
  the vendored tree is the submodule remote, not the host checkout.

## 10. New defects, for the task list

1. **`01_scripts/` vs `01_modules/`** — `bulk-rnaseq-gsea` names a mount point
   that exists in no project, across five files, ten hits.
2. **`02_analysis/helpers/{normalize_gsea,pathway_utils}.R`** — sourced in three
   places, exists nowhere.
3. **`02_analysis/config/config.R`** — sourced in four places, exists in 1 of 8.
4. **`TE-RNAseq-toolkit` labelled `v0.1.0` and `v2.0.x` in one catalog**; actual
   `v2.0.3`.
5. **The `de_*` visualisation family has no skill** (§8).

1–4 are `false today`, fixable now, and blocked on nothing.

---

## 11. The mount-shape check, run early — deletion is safe

`02_pilot.md` asks for this grep before removing
`coresh-signature-search/scripts/*.R`. It is read-only, so I ran it here rather
than let phase 02 start blind.

**One live dependency, and it does not reach through the mount.**
`/data2/users/JCRLab/13403-YD-analysis/13403-YD_Christina` uses all three scripts
in earnest (`remediation-03-coresh-origin.md` documents an mmu CoReSh sweep and
says `extract_gene_loadings.R` "is used verbatim"). It sources
`02_Analysis/helpers/R/coresh/` — its own copy. No `.R` or `.qmd` anywhere in the
fleet sources `.claude/skills/…/scripts/` or
`01_Modules/SciAgent-toolkit/skills/…/scripts/`.

**Verdict: the scripts can be deleted outright.** No shims needed.

But the copy count is the finding. That one project holds **four md5-identical
copies** of the same 269-line kernel:

```
e34d1dc…  01_Modules/SciAgent-toolkit/skills/coresh-signature-search/scripts/   (the mount)
e34d1dc…  02_Analysis/helpers/R/coresh/                                         (what actually runs)
e34d1dc…  01_Modules/.ref/dc-hum-verse-tooling/…/coresh-slice/lib/              (a reference checkout)
e34d1dc…  <toolkit>/skills/coresh-signature-search/scripts/                     (the source)
```

Identical, so nothing has drifted yet — which is luck, not a mechanism.
bulkiRNA's own plan notes the same kernel was md5-identical across four checkouts
before extraction. Five homes for one function is precisely what `coresh_search`
exists to end, and it is the strongest argument in this report for making the
migration real rather than documenting around it.

### Fleet finding, incidental but worth the task list

That project is bound **the pre-#37 way**: `.claude/skills` is a real directory
holding **41 per-skill symlinks** (of 85 skills), each an absolute container path
`/workspaces/13403-YD_Christina/01_Modules/…` — dangling on the host, and 44
skills simply unreachable. Its pin is `5e5347e`. It has no running container, so
it is dormant rather than broken-in-use, and it is a sweep (#30) target: the
single-directory-symlink form plus relative targets fixes both defects at once.

The casing is also non-standard — `01_Modules/`, `02_Analysis/`. Every path claim
in §3 is false there for a reason that has nothing to do with bulkiRNA.
