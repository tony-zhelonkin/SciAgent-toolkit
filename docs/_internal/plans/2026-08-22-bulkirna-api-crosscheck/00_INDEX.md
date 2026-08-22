# Plan INDEX — the catalog does not know bulkiRNA exists

**Date:** 2026-08-22 · **Scope:** IN — cross-check every skill that names a
vendored RNAseq/TE toolkit path against the installed bulkiRNA API, and delete
the implementation depth the package has absorbed / OUT — writing bulkiRNA code,
bumping the image pin, the TE track's own package audit, the fleet sweep

**Origin.** An agent working in `RNAseq-toolkit` reported three owner-gated skill
files (`coresh-signature-search`, `bulk-rnaseq-gsea`,
`annotate-bulk-rnaseq-data`) and correctly said the gate was a decision, not a
technical dependency. Verifying it surfaced a larger and more awkward picture.

> **Phase 01 ran 2026-08-21 — read `01_inventory.report.md`, and read it before
> §1 below.** Two premises in this INDEX were measured wrong. bulkiRNA is not
> version-behind on `v0.5.10`, it is **absent** (7 of 8 live containers); and
> `BULKIRNA_SHA = e42c2de1` is tag **`0.6.0` with 79 exports**, in an image that
> is already built. The report supersedes §1's version table and re-scopes §4.

---

## 1. What the cross-check found

**Thirteen skills name a vendored toolkit path or clusterProfiler. Zero mention
bulkiRNA.** A 59-export package absorbed most of what they describe and the
catalog cannot see it.

```
annotate-bulk-rnaseq-data   bulk-rnaseq-gsea   bulk-rnaseq-pathway-explorer
coresh-signature-search     gatom-metabolomic-predictions
iterative-peak-merging      nfcore-rnaseq-execution
peak-atlas-framework        peak-atlas-multiome
star-te-preprocessing       te-gene-featurecounts
te-geneset-gsea             te-reference-saf-build
```

**Three versions are live at once**, and this is the finding that shapes the
whole plan:

| Layer | State | Has `coresh_*` / `bulkirna_api`? |
|---|---|---|
| bulkiRNA source tree | `1.1.0.9000` at `ca99d8a` (`v1.0.0-1`), 59 exports | yes |
| **installed in the running image** | **`0.4.0`, 64 exports** | **no** |
| what the skills claim | `RNAseq-toolkit` **v0.2.0** script paths | n/a |

The image pins by commit — `BULKIRNA_SHA = e42c2de1…` in
`scbio-docker/docker/base/R/install_core.R:178` — so the API a project sees is a
property of its **image version**, and eight of nine live containers run
`v0.5.10` while `scbio-docker/VERSION` says `v0.5.14`.

**So rewriting a skill to call `coresh_search()` would be wrong for every project
on the current image.** The function does not exist there. This is the
`delegate-cli` lesson arriving in a second place: *a skill must not assert an API
that version-drifts underneath it.* Whatever these skills end up saying has to
survive being read inside a container one or two image versions behind.

## 2. The principle to apply

**A skill owns *when* and *why*, plus the traps that cannot be expressed in code.
The package owns *how*.** Where the API absorbed the how, the skill deletes it
and routes to the verb. Where the trap is knowledge the package does not hold —
the CORESH Entrez-integer trap, species mismatch, `eps = 0` truncating p-values —
it stays. That is the non-regenerable half of the asset and it does not shrink.

Test for every line: *does the package now do this? Then the line is a
re-implementation and goes. Does the package silently get it wrong unless the
reader knows something? Then the line is the point of the skill and stays.*

## 3. The mapping, as far as reading establishes it

High-confidence supersessions, source-tree API:

| Skill ships / cites | Absorbed by |
|---|---|
| `coresh-signature-search/scripts/coresh_batch.R` | `coresh_search` |
| `…/extract_gene_loadings.R` | `coresh_loadings`, `gsdb_coresh` |
| `…/symbols_to_entrez.R` | `gene_to_entrez`, `entrez_to_gene` |
| `…/validate_coresh_install.R` | `coresh_validate` |
| GSEA execution via `clusterProfiler::GSEA()` | `gs_test` |
| `references/master-tables.md`, `normalize_gsea.R` | `gs_to_master`, `gs_master_columns`, `gs_validate_master`, `gs_stat_types` |
| `references/visualization.md` | `gs_plot_{bar,dot,heatmap,running}`, `theme_bulki` |
| `references/msigdb.md` | `gsdb_msigdb`, `gsdb_load` |
| `references/custom-db.md` | `gsdb_from_file`, `gsdb_register`, `gsdb_list`, `gsdb_info` |
| `RNAseq-toolkit v0.2.0 scripts/General/{io_helpers,annotate_genes,dge_helpers,provenance}.R` | `read_counts_matrix`, `read_metadata`, `annotate_genes`, `build_dge`, `write_session_provenance`, `ensure_dir` |
| `gatom-metabolomic-predictions:101` `source(.../load_reference_db.R)` + `download_gatom_references()` | `gatom_download_refs`, `gatom_refs`, `gatom_de`, `gatom_genes`, `gatom_module`, `gatom_save_html` |

Also unclaimed by any skill: `coresh_chunks`, `coresh_match`, `coresh_sets`,
`coresh_convergence`, `gs_coregulation`, `gs_leading_edge`, `gs_ranks`,
`gs_score`, `de_volcano`, `de_volcano_grid`, `de_pca`, `de_pca_3d`, `de_md_plot`,
`de_bfc_plot`, `filter_confounder_genes`. Capability the catalog never learned
about — worth a pass in the other direction.

**The registry that should have made this mechanical is empty.** `bulkirna_api()`
returns a `superseded_by` column and `R/api.R:130` fills it with
`rep(NA_character_, …)`. The schema exists; the data does not. So this mapping is
authored by reading, and it will have to be re-authored next time.

## 4. Phases

| # | Ships | Repo | Parallel? |
|---|---|---|---|
| 01 | the three-way inventory: skill claim × source API × installed API, per skill | scio | first |
| 02 | rewrite the three highest-confidence skills as the pilot | scio | after 01 |
| 03 | the remaining ten, triaged by 01 | scio | after 02 |

`01_inventory.md` and `02_pilot.md` carry the briefs. Phase 03 is deliberately
unwritten: phase 01 decides whether ten more skills need this or whether some of
those thirteen name `01_scripts/` for an unrelated reason (the three ATAC skills
probably do).

## 5. Sequencing — this is the part that matters

**Land phases 01–02 before the fleet sweep (#30), or they wait for the next
delivery event.** Today's delegation trace measured that cost exactly: a correct
fix sat in this repository for **13 h 25 min** while the project that needed it
ran the stale copy, because nothing re-pins a vendored submodule mid-session. A
skill fix that misses the sweep misses 25 projects.

Collision risk is real but small: today's work touched
`skills/reasoning-trace` — none of the thirteen. Check `git log --oneline -5`
before starting.

## 6. Recommendations back to the RNAseq-toolkit side

Not this repo's work, but this plan is the evidence for both:

1. **Populate `superseded_by` in `bulkirna_api()`.** The package knows what it
   absorbed; scio can only guess. A populated registry turns this cross-check
   from an afternoon of reading into a join, and it is the right home for the
   mapping — the same "one home per claim" rule this toolkit runs on.
2. **The image pin is the real gate.** `coresh_*` and `bulkirna_api` reach a
   project only when `BULKIRNA_SHA` advances and the image is rebuilt and the
   containers recreated. Deciding to rewrite the skills implies deciding to move
   that pin; otherwise the skills describe an API their readers do not have.

## 7. Out of scope

- Writing or changing bulkiRNA code.
- Bumping `BULKIRNA_SHA` or building an image.
- The TE track (`TE-RNAseq-toolkit`, five skills) needs its own API audit against
  its own package; phase 01 records what it finds and stops there.
- The fleet sweep.
