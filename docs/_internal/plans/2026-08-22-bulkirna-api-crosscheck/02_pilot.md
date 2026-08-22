# Phase 02 — rewrite three skills as the pilot

**Repo:** scio · **Blocked by:** phase 01's floor decision · Read `00_INDEX.md`
§2 for the principle.

> **Phase 01 answered this brief's two open questions on 2026-08-21.** Read
> `01_inventory.report.md` §6 and §11 before starting.
>
> 1. **The floor decision landed on option 3 below — wait for the image.** The
>    package is not installed at all on `v0.5.10`, which is 7 of 8 live
>    containers, so no bulkiRNA text is executable for them. `v0.5.14` carries
>    the 0.6.0 API and is already built. **This phase is blocked on recreating
>    the containers, not on a scbio-docker code change.**
> 2. **The mount-shape grep is done. Deleting `scripts/*.R` is safe** — the one
>    live consumer sources its own copy, and no file in the fleet reaches through
>    the mount. No shims needed.
>
> Phase 01 also split this phase in two. **Phase 02a below is unblocked and
> should ship now**; the rewrite waits.
>
> The floor set, when it unblocks, is `0.6.0 ∩ HEAD` — 58 functions. **Never name
> one of the 21 legacy exports** (`run_gsea`, `load_reference_db`,
> `normalize_gsea_results`, `download_gatom_references`, the `gsea_*` plots …):
> they are live in the pinned image and already deleted upstream, so writing them
> buys one image version and then breaks.

---

## Phase 02a — the false claims, unblocked

Twenty-two claim hits name paths that exist in **no project**, under every image
version. They are wrong today, their fix names no package function, and they gate
on nothing. Ship them independently of everything above.

1. **`01_scripts/RNAseq-toolkit/…` → `01_modules/RNAseq-toolkit/…`** — ten hits
   across five `bulk-rnaseq-gsea` files (`custom-db:96`, `master-tables:405,407`,
   `msigdb:318,319,399`, `visualization:435,436,437`, `SKILL:196`). All eight live
   projects use `01_modules/`. One skill family, two spellings of one mount.
2. **`source("02_analysis/helpers/normalize_gsea.R")` and `…/pathway_utils.R`** —
   `master-tables:56,104,105`, `custom-db:107`. Neither file exists in any
   project; `02_analysis/helpers/` is empty in seven of eight.
3. **`source("02_analysis/config/config.R")`** — `SKILL:76`, `msigdb:50`,
   `master-tables:104`, `visualization:61,62`. Present in 1 of 8 (14616-DM).
   Either stop asserting it or state it as a project-local convention the reader
   supplies.
4. **`TE-RNAseq-toolkit` version labels** — `annotate-bulk-rnaseq-data` and
   `te-geneset-gsea` say `v0.1.0`; `star-te-preprocessing:68,156` and
   `te-reference-saf-build:224` say `v2.0.0`/`v2.0.1`; the repo is at `v2.0.3-5`.
   Phase 01 verified every cited file, the hardcoded-`source()` trap
   (`create_te_genesets.R:109,175`) and the `gs_name`/`gene_symbol` contract all
   survive the jump — so this is a **label correction only**. Do not touch the
   substance, and do not start the TE track's own audit.

Same gates as below. This is the `docs-layout` lesson in the catalog: an
instruction naming something the mechanism does not guarantee.

---

## The three

The subset the RNAseq-toolkit agent named, and the highest-confidence rows in the
inventory:

1. **`coresh-signature-search`** — ships four scripts under `scripts/`, all
   superseded: `coresh_batch.R` → `coresh_search`, `extract_gene_loadings.R` →
   `coresh_loadings`/`gsdb_coresh`, `symbols_to_entrez.R` →
   `gene_to_entrez`/`entrez_to_gene`, `validate_coresh_install.R` →
   `coresh_validate`. The SKILL.md names them at nine places.
2. **`bulk-rnaseq-gsea`** — four reference files describing a wrapper around
   `clusterProfiler::GSEA()` plus `01_scripts/RNAseq-toolkit/...` paths. The
   `gs_*` and `gsdb_*` layers cover execution, master tables, plots and database
   registration.
3. **`annotate-bulk-rnaseq-data`** — `SKILL.md:71` pins **RNAseq-toolkit v0.2.0**
   by name and lists four script files, all now exports.

## What to keep, explicitly

Deleting the wrong half destroys the asset. **Keep every trap.** Named examples
that must survive verbatim in substance:

- the **CORESH Entrez-integer trap** and the assertion that catches it
  (`coresh-signature-search` ~line 193) — the package cannot warn you that your
  symbols silently matched nothing;
- **species mismatch** and the ortholog remedy (~199);
- **`eps = 0`** or p-values truncate at ~1e-4 (`bulk-rnaseq-gsea:173`);
- the **`p.adjust` vs `qvalue`** divergence between GSEA runs
  (`references/master-tables.md:396`) — unless phase 01 shows `gs_stat_types`
  resolves it, in which case route to that and say so;
- timing figures that set expectations (~10–20 s variance-only, 2–5 min with
  p-values).

## Shape after the rewrite

*When to reach for this, why, the traps, and which verb does it.* Not how to
implement it. A skill that used to teach a procedure becomes a skill that teaches
a judgement and names a function.

Frontmatter `description` stays ≤ 350 chars and must stop promising what the file
no longer contains — `bulk-rnaseq-gsea`'s currently advertises "MSigDB execution
via clusterProfiler/fgsea".

## The deletion — checked, and clear

Removing `coresh-signature-search/scripts/*.R` changes **mount shape**, so phase
01 ran the fleet grep first. Result in `01_inventory.report.md` §11:

- One live consumer, `13403-YD_Christina`, runs all three scripts in earnest —
  from its own copy at `02_Analysis/helpers/R/coresh/`.
- **Nothing in the fleet sources the mount path** (`.claude/skills/…/scripts/` or
  `01_Modules/SciAgent-toolkit/skills/…/scripts/`).

**Delete them outright. No shims.** Do re-run the grep if the sweep lands first —
it changes what is bound where.

Worth carrying into the rewrite's prose: that one project holds four
md5-identical copies of the same 269-line kernel. Five homes for one function is
the condition `coresh_search` exists to end.

## The floor — phase 01 chose option 3

It found worse than "absent from the live image": the whole package is absent from
`v0.5.10`, and `coresh_*` reaches **0 of 8** containers. Option 3 it is. The
options are kept below because option 1 becomes the answer if the owner decides
the fleet stays split.

1. Say the version the verb requires, and give the one-line check the reader can
   run: `"coresh_search" %in% getNamespaceExports("bulkiRNA")`. **Probe, do not
   assert** — the same rule `delegate-cli` now follows for codex flags.
2. Keep both paths, package-first, with the script path named as the fallback for
   older images. Costs length, and length is what this plan is trying to reduce.
3. Wait for the image pin. Then this phase blocks on a scbio-docker change and
   should say so rather than shipping text nobody can run.

Do not invent a version-detection helper for scio to ship. That is a fourth
delivery mechanism for one package's problem.

## Gates

```bash
bash tests/run-all.sh                     # 65 passing, 0 failing
bin/scio lint --check toolkit --strict    # clean; descriptions ≤ 350 chars
bash tests/test_mount_layout_identity.sh  # if any scripts/ file moved
bash tests/test_skill_link_integrity.sh
```

## Do not

- Do not write or change bulkiRNA code, and do not touch `BULKIRNA_SHA`.
- Do not touch the five TE skills or the three ATAC skills.
- Do not delete a trap because the API "probably handles it". Verify or keep.
- Do not run the fleet sweep. But do finish before it — see `00_INDEX.md` §5.
