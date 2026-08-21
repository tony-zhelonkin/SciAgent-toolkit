# 2026-06-17 — Annotate-skill refactor: replicate the reference gene+TE annotation style, reproducibly

**Author:** architect (SciAgent-toolkit)
**Status:** design — no code written yet
**Scope of this doc:** design + ordered execution plan for refactoring the gene+TE annotation
skills (`annotate-bulk-rnaseq-data`) and their sibling-repo helpers (RNAseq-toolkit,
TE-RNAseq-toolkit), plus scaffolding a new `te-geneset-gsea` skill.

---

## 1. Problem statement

The `annotate-bulk-rnaseq-data` skill is a **thin router**: its two reference docs cite
helper functions that live in two *external* sibling repos (RNAseq-toolkit, TE-RNAseq-toolkit),
not in SciAgent-toolkit. The router prose, the sibling helper code, and the user's canonical
reference annotation script have drifted apart. The router cannot currently reproduce the
user's annotation style, and several citations point at code paths that do not match the spec
(or do not exist at all).

The user's canonical reference is `13036/02_Analysis/00_annotate_data.R` (a single self-contained
Rscript). The refactor must make the skill + helpers reproduce that script's *outputs and
semantics* generically (config-driven, not 13036-hardcoded), while preferring to **repoint/wire**
existing correct code over rewriting it.

### Ground-truth corrections to the gap brief (trust these)

- **No `R/` dir in RNAseq-toolkit.** It is a script library; gene code lives under
  `scripts/General/{io_helpers,annotate_genes,dge_helpers}.R`. The gap text's "RNAseq-toolkit R/
  factory fns" phrasing refers to the **TE**-RNAseq-toolkit `R/` dir, not RNAseq-toolkit.
- **No literal `SYMBOL` column in TE-RNAseq-toolkit `R/`.** The TE geneset gene slot column is
  named `gene_symbol`, and it equals `subfamily` only when `use_subfamily_as_gene = TRUE`
  (the default). The spec's "SYMBOL=subfamily" is conceptually right but not literal.
- **Two `te_utils.R` files exist — thin (`scripts/`) vs rich (`R/`).** The skill currently sources
  the THIN `scripts/te_utils.R` (verified present, 940 bytes): its `build_te_annotation()` emits
  `Symbol, Ensembl, subfamily, family, class, type` (TE label in `Symbol`, `Ensembl=NA`,
  `type="TE"`). The RICH `R/te_utils.R::build_te_annotation()` emits
  `feature_id, is_te, subfamily, family, class, replication_type`. The repoint (T1) is thin→rich,
  not the removal of a dangling reference (`scripts/te_utils.R` is real). `feature_type` is NOT
  emitted by either `build_te_annotation`; it is added by `create_combined_dge()` /
  `annotate_dge_features()` (verified `R/create_combined_dge.R:178`).
- **Version tags.** SciAgent-toolkit's own tags are `v3.0.0`, `v3.1.0` only — it has no `.gitmodules`
  and does not vendor the sibling toolkits. The annotate skill's `RNAseq-toolkit v2.0.0` /
  `TE-RNAseq-toolkit v2.0.1`/`v2.0.0` pins are *phantom*. The real released tags are **RNAseq-toolkit
  `v0.2.0`** and **TE-RNAseq-toolkit `v0.1.0`** (verified: annotated tag "Release v0.1.0: initial
  tagged release … for submodule pinning", on `main`). The `v0.1.0` tag ALREADY contains the rich
  `R/` factories (`te_utils.R`, `create_combined_dge.R`, `create_te_genesets.R`,
  `validate_te_input.R`) byte-identical to the working tree (`git diff v0.1.0 -- R/<f>` = 0 lines for
  all four), and `git log v0.1.0..HEAD` is docs-only. The editable checkout sits on
  `feat/strand-split-qc`, but the branch a checkout happens to be on is irrelevant to which tag a
  skill doc cites — pin to the tag `v0.1.0` (NOT the in-file `Version: 2.0.0` header, which
  corresponds to no git tag and would re-introduce a phantom pin). All TE pins reconcile to `v0.1.0`.
- **TE-RNAseq-toolkit is NOT a submodule of 14839** and the survey's editable checkout is in a
  separate pipeline tree. TE code is editable there, but wiring it into any analysis project is a
  *project* concern (W1), out of scope for this toolkit doc.

### User-style target (what the refactor must let an agent reproduce)

From `13036/02_Analysis/00_annotate_data.R` and the TE reference semantics:

**Gene path**
- Ensembl `ENSMUSG` is the **primary key**. The matrix is read with versioned rownames, version is
  stripped (`sub("\\..*$","")`) for mapping, and `agg_rows_by_id()` REBUILDS the matrix with the
  **stripped stable IDs** as rownames (reference lines 52/55). So the canonical DGEList rownames,
  `gene_index` key, and `input_gene_name` join key are all the **stripped stable ID** (post-aggregate),
  NOT versioned. Do not tell the recipe to "keep versioned rownames" — it would mis-key the join.
- **SUM-aggregate** duplicate stable IDs (the script's `agg_rows_by_id()` — collapse rows that
  share a stripped Ensembl ID by `colSums`). Output rownames come out in `split()` (alphabetical-by-ID)
  order, NOT input order; the reference survives this only by re-matching annotation to
  `rownames(count_mat)` afterward (lines 219/267 `match()`). Any recipe that aggregates MUST re-order
  annotation to `rownames(mat)` before building the DGEList/writers (see S3/S6 alignment checks).
- `org.Mm.eg.db` `mapIds()` → `SYMBOL` + `ENTREZID`.
- biomaRt is **optional/best-effort** (`try()`, `TRY_BIOMART` toggle) → adds `gene_biotype` +
  `mgi_symbol`; on failure these stay `NA`/fallback.
- `SYMBOL` falls back to the Ensembl ID when missing/empty.
- Integer counts with a `ROUND_NONINT` toggle; `stopifnot` row/col alignment.
- Salmon gene-level ingestion captures `gene_name` as `input_gene_name` (featureCounts path sets it
  `NA`).
- Canonical deliverable = edgeR `DGEList(counts + genes + samples + TMM)` saved as
  `<PROJECT_ID>_DGEList.rds`.

**Signature exports** — the reference writes THREE tables (do not conflate):
- (A) "Wide-by-sample matrix" `<PROJECT_ID>_counts_wide_by_sample.tsv` (reference lines 233-239):
  rows = samples, columns = factor columns + gene columns. **In scope for S6.**
- (B) "Transposed annotated matrix" `<PROJECT_ID>_counts_transposed.tsv` (reference lines 262-299):
  lead columns `Symbol, Ensembl`, then sample columns, with the experiment **factor rows embedded at
  the TOP**. Exact layout (reference line 277): factor name in the FIRST annotation column (`Symbol`),
  the SECOND annotation column (`Ensembl`) BLANK (`""`), and ALL other annotation columns blank for the
  top rows; factor values across the sample columns. `nrow == length(factor_cols) + nrow(counts)`.
- (C) Gene dictionary `<PROJECT_ID>_gene_index.tsv` with columns
  `ensembl_gene_id, SYMBOL, mgi_symbol, ENTREZID, gene_biotype, input_gene_name` (reference line 308+).
  NOTE the column-name rename: `annotate_genes_from_ensembl()` emits `Symbol`/`Ensembl` (capital S /
  E), so the writer maps `Ensembl→ensembl_gene_id` and `Symbol→SYMBOL` (NOT a straight projection).

**TE path**
- featureCounts SAF IDs `Subfamily:Family:Class` → parse into `subfamily/family/class`;
  `feature_type = "te"`, `is_te = TRUE`, gene/symbol slot = `subfamily`, `replication_type`
  (cytoplasmic if class ∈ {LINE, SINE, LTR, Retroposon}; nuclear if ∈ {DNA, RC}).
- `rbind` genes + TEs into **ONE combined `DGEList`** (TMM on the combined library; genes dominate).
- Family-level + class-level **GMT genesets** for TE GSEA.

---

## 2. Repo topology — which repo owns each fix, and the editable checkout to edit

> The skill is a router. Fixes split across three repos. The table maps every gap to its repo +
> the exact editable checkout an implementer must edit. "BLOCKED" = no editable, project-wired
> checkout reachable from this toolkit's contract.

| Gap | Owning repo | Editable checkout to edit | Editable? |
|-----|-------------|---------------------------|-----------|
| G1 gene_index export | RNAseq-toolkit | `/scratch/current/antonz/projects/14839-DM-cGAS/01_modules/RNAseq-toolkit` (branch `dev`, HEAD = `v0.2.0`) | yes |
| G2 project-ID filenames from config | RNAseq-toolkit | same as G1 | yes |
| G3 Salmon gene_name capture | RNAseq-toolkit | same as G1 | yes |
| G4 parameterize metadata factor rows | RNAseq-toolkit | same as G1 | yes |
| G5 `aggregate_duplicate_ids()` add / doc fix | RNAseq-toolkit | same as G1 | yes |
| G6 reproducibility (pin Ensembl, sessionInfo, build) | RNAseq-toolkit | same as G1 | yes |
| T1 repoint TE skill at rich `R/` factories | SciAgent-toolkit (doc) | `/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit` (branch `dev`) | yes |
| T1b TE `R/` `source()` hardcoded path | TE-RNAseq-toolkit | `/data1/users/antonz/pipeline/TE-RNAseq-toolkit` (HEAD `feat/strand-split-qc`; cite tag `v0.1.0`) | yes (but see note) |
| T2 combined gene+TE DGEList via factory | SciAgent-toolkit (doc) | same as T1 | yes |
| T3 family/class GMT handoff to GSEA | SciAgent-toolkit (doc) | same as T1 | yes |
| T4 reconcile phantom version pins | SciAgent-toolkit (doc) | same as T1 | yes |
| N1 new `te-geneset-gsea` skill | SciAgent-toolkit (new skill) | same as T1 | yes |
| W1 add TE-RNAseq-toolkit submodule to 14839 | analysis project `14839-DM-cGAS` | NOT a toolkit path | OUT OF SCOPE |

**Editable-checkout notes**
- **RNAseq-toolkit:** the 14839 submodule checkout is on branch `dev` at the `v0.2.0` tag and is
  writable. Editing it here advances `dev`; the released `v0.2.0` tag will lag and a later
  re-tag + submodule bump is needed (see Risks). All gene fixes (G1–G6) land here.
- **TE-RNAseq-toolkit:** the only fully-populated editable checkout is in the pipeline tree, on
  branch `feat/strand-split-qc` with 2 uncommitted doc edits in its working tree. The released tag
  `v0.1.0` (on `main`) carries the rich `R/` factories byte-identical to that working tree, so the
  skill pins `v0.1.0` regardless of which branch the checkout sits on. It is NOT a submodule of any
  analysis project relevant to this toolkit. Therefore TE *code* edits (T1b/S12) are technically
  editable but **operationally BLOCKED**: changing code in a checkout that no consuming project wires
  in changes nothing reproducible.
- **T1b is NOT neutralized by pre-sourcing — corrected.** R's `source()` opens a file BY PATH; it does
  not skip because the symbols are already in scope. `R/create_te_genesets.R:109` (inside the
  `level == "class"` branch) and `:175` (`create_te_genesets_from_dge`) UNCONDITIONALLY run
  `source("01_modules/TE-RNAseq-toolkit/R/te_utils.R")` (verified). When the recipe's cwd is not a
  root where that exact relative path resolves (which is the case from the 14839 root — TE is not a
  14839 submodule, W1), `create_te_genesets(level = "class")` HARD-ERRORS with "cannot open file …".
  `level = "family"` has NO internal `source()` and works. **Consequence:** the class-level GMT path
  is genuinely blocked unless one of these holds, so S7 carries a class-level mitigation (see S7
  change text): build the class T2N in the recipe from the already-sourced `TE_CYTOPLASMIC_CLASSES` /
  `TE_NUCLEAR_CLASSES` constants WITHOUT calling `create_te_genesets(level = "class")`, OR `setwd()`
  to a root where the relative path resolves, OR defer the class GMT to S12+W1. The upstream
  self-locating `source()` fix is a recorded TE-RNAseq-toolkit follow-up (S12, BLOCKED).

---

## 3. Per-gap solution

### Gene path (RNAseq-toolkit `scripts/General/`)

**G1 — Missing `*_gene_index.tsv` dictionary export.**
Add a helper `write_gene_index(ann_df, outfile)` (home = the NEW `provenance.R`, per S5 — NOT
`io_helpers.R`) that writes exactly
`ensembl_gene_id, SYMBOL, mgi_symbol, ENTREZID, gene_biotype, input_gene_name` (column order from
the reference script). This is NOT a straight projection: `annotate_genes_from_ensembl()` emits
`Symbol`/`Ensembl` (capital S / capital E), so `write_gene_index` must RENAME inside the helper —
`ensembl_gene_id = Ensembl`, `SYMBOL = Symbol` — then carry `mgi_symbol, ENTREZID, gene_biotype,
input_gene_name` through. Cite the helper in `gene-annotation.md`.

**G2 — No project-ID-prefixed filenames driven by config.**
The skill's `gene-annotation.md` already has a "paths (customize per project)" block. Make it read
`project.id` from `02_analysis/config/analysis_config.yaml` (via `yaml::read_yaml`) and build output
names as `file.path(outdir, paste0(project_id, "_DGEList.rds"))`,
`paste0(project_id, "_gene_index.tsv")`, `paste0(project_id, "_counts_transposed.tsv")`. This is a
**reference-doc** change (the recipe), not a toolkit function — filenames are caller policy, the
toolkit just provides the writers. Record the config key as `project.id`.

**G3 — No Salmon gene-level `gene_name` capture (`input_gene_name`).**
Extend `read_counts_matrix()` to detect the Salmon gene-level shape
(`gene_id` + `gene_name` + samples) in addition to the featureCounts shapes it already handles, and
to return the captured `gene_name` vector (e.g. as an attribute `attr(mat, "input_gene_name")`,
keyed by stripped stable ID, `NA` for featureCounts). Then `annotate_genes_from_ensembl()` gains an
`input_gene_name = NULL` argument that is threaded into the returned tibble (column
`input_gene_name`, `NA` when absent). This mirrors the reference script's auto-detect (`fmt`) +
`input_gene_name` plumbing. **Cross-step hazard:** `aggregate_duplicate_ids()` (G5) rebuilds the
matrix and DROPS attributes, so it must re-subset `attr(mat,"input_gene_name")` to the collapsed
unique-id order (first non-NA per group) or the captured names are lost before annotation (see S3).

**G4 — `write_annotated_matrix()` factor rows HARDCODED to 13036; `read_metadata()` stops without
them.**
Generalize both:
- `read_metadata(fp, sample_col_candidates = c("Sample_ID","Sample ID"), required_cols =
  character(0))` — drop the hardcoded 8-column `needed` `stop()`; keep only the Sample-ID resolution.
  Required columns become an *opt-in* argument (empty by default) so non-13036 sheets load. Also
  accept CSV as well as `.xlsx` (the reference reads `13036-DM_Samples_Metadata.csv` via
  `data.table::fread`, line 136; the current toolkit hardcodes `readxl::read_xlsx`) — dispatch on
  extension, or document that the recipe reads the CSV directly.
- `write_annotated_matrix(mat, md, add_cols, outfile, factor_cols = NULL)` — replace the hardcoded
  `top_map` list with a generic loop over `factor_cols` (a character vector of metadata column
  names). Each factor name goes in the FIRST annotation column, the SECOND annotation column blank
  (`""`), all other annotation columns blank for the top rows, values across sample columns (the
  reference `factor_rows` pattern, line 277). When `factor_cols = NULL`, write no top block (plain
  annotated matrix).
This makes the writer reproduce the reference "transposed annotated matrix with factor rows at TOP"
while being design-agnostic. Factor column names come from the caller; factor DERIVATION (e.g.
`LPS_treat`/`IFNg_treat`/`genotype` from `Group`/`Treatment` via `if_else`, reference lines 165-172)
is **recipe/caller responsibility** done before `write_annotated_matrix(factor_cols=)` — the toolkit
function does not derive factors, it only embeds the named columns.

**G5 — Doc bug: `gene-annotation.md` calls `aggregate_duplicate_ids()` which is absent.**
**Resolve by ADDING the function** (the reference script proves the user wants SUM-aggregation, so
the doc's intent is correct). Add `aggregate_duplicate_ids(mat, ids = rownames(mat))` to
`io_helpers.R` implementing the reference `agg_rows_by_id()` SUM-collapse: split row indices by
`ids`, `colSums` each group, rebuild the matrix with the unique IDs as rownames; if no dups, just
set rownames. Two contracts the implementer MUST honor:
(1) output rownames come out in `split()` order = `sort(unique(ids))`, NOT input order, so the recipe
must re-match annotation to `rownames(mat)` downstream (reference lines 219/267);
(2) preserve `attr(mat,"input_gene_name")` (from G3/S1) by re-subsetting it to the collapsed unique-id
order (first non-NA per group) — `colSums` rebuild would otherwise drop it.
Keep the `gene-annotation.md` citation as-is (now valid).

**G6 — Reproducibility: unpinned `useEnsembl()`, no `sessionInfo`, genome build unrecorded.**
- Add a `biomart_version = NULL` (and/or `biomart_host = NULL`) argument to
  `annotate_genes_from_ensembl()`, passed through to `useEnsembl()` so the Ensembl release can be
  pinned. Default `NULL` preserves current floating behavior. **`reference.ensembl_version` does NOT
  currently exist** in the 14839 config (verified: only `project.id`, `project.genome_build`,
  `project.species_db`). S6 must therefore EITHER add `reference.ensembl_version` (+ optional
  `biomart_host`) to the config template so the pinned path is the default, OR state plainly that the
  key is absent and the run stays floating with provenance-only reproducibility (Risk 4).
- Add a tiny `write_session_provenance(outfile, genome_build = NULL, ensembl_version = NULL)` helper
  (new `scripts/General/provenance.R`, per S5) that writes `sessionInfo()` PLUS the genome build /
  Ensembl release. Critically, `sessionInfo()` records the biomaRt PACKAGE version, NOT the live
  remote archive release that actually varies — so the recipe must ALSO capture the resolved release
  (e.g. via `biomaRt::listEnsemblArchives()` / the archive host actually used) into the provenance
  file, not just `sessionInfo()`.
- The recipe records `project.genome_build` (mm10 here) into the provenance file and the
  `gene_index` header comment. Genome-build recording is a *recipe* responsibility tied to config.

### TE path

**T1 — Skill routes to the THIN `scripts/te_utils.R`; the RICH `R/` factories already match the
spec but are uncited. REPOINT the skill (do not rewrite code).**
Rewrite `te-annotation.md` to `source()` the rich `R/` files and call their factories:
- `R/te_utils.R` — `parse_te_id`, `build_te_annotation` (emits `feature_id, is_te, subfamily,
  family, class, replication_type` — NOT `feature_type`), `annotate_dge_features`,
  `detect_feature_type`, `is_te`, and the constants `TE_CYTOPLASMIC_CLASSES`/`TE_NUCLEAR_CLASSES`.
- `R/validate_te_input.R` — `validate_combined_input()` (a HARD dependency of `create_combined_dge`
  with its default `validate = TRUE`; verified the function is defined only here).
- `R/create_combined_dge.R` — `create_combined_dge()` (T2) — this is where `feature_type` ("gene"/"te")
  is added to `$genes`.
- `R/create_te_genesets.R` — `create_te_genesets()`, `filter_te_genesets()`,
  `export_te_genesets_gmt()` (T3).
Remove the thin-path citation (the real `scripts/te_utils.R::build_te_annotation`, which emits
`Symbol/Ensembl/type`) and the hand-rolled `rbind`/`Symbol = Symbol` recipe the current doc shows —
they contradict the rich factories (rich emits `is_te`/`replication_type` from `build_te_annotation`,
and `feature_type` from `create_combined_dge`). Document that the gene/symbol slot for a TE is the
**`subfamily`** (via `create_te_genesets(..., use_subfamily_as_gene = TRUE)`), naming the column
accurately as `gene_symbol`.

**T1b (upstream, recorded) — hardcoded `source()` in `create_te_genesets.R` is NOT neutralized by
pre-sourcing.**
`R/create_te_genesets.R:109` (inside the `level == "class"` branch) & `:175`
(`create_te_genesets_from_dge`) UNCONDITIONALLY run `source("01_modules/TE-RNAseq-toolkit/R/te_utils.R")`,
a cwd-relative path with no `exists()` guard. R's `source()` opens the file by PATH regardless of
whether the symbols are already in scope, so pre-sourcing does NOT make it harmless: when the relative
path does not resolve from the recipe's cwd (the 14839 root case — TE is not wired in, W1),
`create_te_genesets(level = "class")` HARD-ERRORS. `level = "family"` has no internal `source()` and
works. **Therefore the S7 recipe does NOT rely on pre-sourcing to neutralize this.** S7 obtains the
class-level GMT by ONE of: (a) building the class T2N in the recipe from the already-sourced
`TE_CYTOPLASMIC_CLASSES`/`TE_NUCLEAR_CLASSES` constants (no `create_te_genesets(level="class")` call);
(b) `setwd()` to a root where `01_modules/TE-RNAseq-toolkit/R/te_utils.R` resolves; or (c) deferring
the class GMT until W1+S12. The proper upstream fix (make `R/` self-locating via `here`/`system.file`)
is the TE-RNAseq-toolkit follow-up S12 (BLOCKED), owned by the TE maintainer.

**T2 — Combined gene+TE output is a flat TSV; user wants ONE combined `DGEList`.**
Replace the `te-annotation.md` "Combined gene + TE matrix" `rbind`-of-TSVs recipe with a
`create_combined_dge(gene_counts, te_counts, samples, gene_annotation = ann_gene_df, te_annotation =
build_te_annotation(rownames(mat_te)))` call. **`create_combined_dge` defaults `validate = TRUE`**,
which calls `validate_combined_input()` (defined only in `R/validate_te_input.R`) — so the recipe MUST
`source("R/validate_te_input.R")` (per T1) OR pass `validate = FALSE` explicitly, else it errors with
"could not find function validate_combined_input". It rbinds counts, builds the `$genes` slot with
`feature_type/is_te/subfamily/family/class/replication_type`, sets `analysis_mode = "combined"`, and
TMM-normalizes. Save `<PROJECT_ID>_combined_DGEList.rds`. Keep the existing graded joint-DE caveats
block (already present, defers to `te-gene-featurecounts` evidence scale).

**T3 — No family/class GMT geneset handoff to GSEA.**
In `te-annotation.md`, after building the TE annotation, call (EXACT verified signatures —
`te_annotation` is a required positional arg; `export_te_genesets_gmt(te_genesets, output_file)` takes
the genesets LIST, not a level):
- Family GMT: `gs_fam <- create_te_genesets(te_ann, level = "family")` then
  `export_te_genesets_gmt(gs_fam, "<PROJECT_ID>_TE_family.gmt")`. Family has no internal `source()` and
  runs anywhere.
- Class GMT: do NOT call `create_te_genesets(te_ann, level = "class")` from the 14839 root — its
  `R/create_te_genesets.R:109` internal `source()` hard-errors (T1b). Build the class T2G/T2N IN THE
  RECIPE: `T2G <- distinct(select(te_ann, gs_name = class, gene_symbol = subfamily))`, derive
  descriptions from the already-sourced `TE_CYTOPLASMIC_CLASSES`/`TE_NUCLEAR_CLASSES`, then
  `export_te_genesets_gmt(list(T2G = T2G, T2N = T2N), "<PROJECT_ID>_TE_class.gmt")`. (Or `setwd()` to a
  resolving root, or defer to S12+W1.)
These GMTs become the documented input to the new `te-geneset-gsea` skill (N1). Note the small set
count (~dozens) → relaxed `min_size` (do NOT apply the default 5–500 filter blindly).

**T4 — Reconcile phantom version pins to REAL git tags.**
- `annotate-bulk-rnaseq-data/SKILL.md:103` `RNAseq-toolkit v2.0.0` → **`v0.2.0`** (real tag).
- `SKILL.md:104` `TE-RNAseq-toolkit v2.0.0` AND `te-annotation.md` `v2.0.1` (4 places) → reconcile to
  the **real annotated tag `TE-RNAseq-toolkit v0.1.0`** (on `main`; verified to contain the cited rich
  `R/` factories byte-identical to the working tree). Do NOT pin to the in-file `Version: 2.0.0`
  header (no git tag → would re-introduce a phantom pin) and do NOT pin to a feature branch (the
  branch a checkout sits on is independent of which tag a doc cites). Use the single consistent string
  `TE-RNAseq-toolkit v0.1.0` everywhere (SKILL.md Resources line + every te-annotation.md mention).
- `gene-annotation.md` `v2.0.0` mentions → `v0.2.0`.

**N1 — New `te-geneset-gsea` skill.**
`bulk-rnaseq-gsea` is gene/MSigDB/symbol-centric; the TE-GMT → GSEA handoff is homeless. Scaffold a
new `implementation`-scope, `standard`-tier skill `te-geneset-gsea` that:
- Takes the family/class GMTs from T3 + a ranked TE list (DE on the TE/combined DGEList, ranked by
  the TE **subfamily** IDs — the same ID space as the GMT `gene_symbol` column).
- Reuses the `bulk-rnaseq-gsea` T2G/T2N → `clusterProfiler::GSEA()` machinery at the ALGORITHM level
  — so `te-geneset-gsea` is a thin specialization, not a fork. **The column names DIFFER and snippets
  must be adapted, not copied:** TE GMTs use `gs_name`/`gene_symbol` (verified `create_te_genesets.R`),
  whereas the bulk-rnaseq-gsea `@geneSets` snippet uses `T2G$gene`/`T2G$term`. So the TE `@geneSets`
  fix is `gsea_result@geneSets <- split(T2G$gene_symbol, T2G$gs_name)` and the overlap check is
  `sum(db$T2G$gene_symbol %in% names(ranked))`. A verbatim copy of the bulk snippet references
  non-existent `$gene`/`$term` columns and silently yields empty gene sets (the 0-results trap N1
  exists to prevent).
- Documents the TE-specific traps: ID space must be subfamily (not Ensembl symbols — the #1 GSEA
  pitfall is non-symbol rownames → 0 results); the overlap check (`<30%` → wrong ID space); the
  `@geneSets` slot fix for `enrichplot`; and a **relaxed `min_size`** (the default 5 drops most TE
  classes — only ~6 classes exist). Complementary to `bulk-rnaseq-gsea` and
  `annotate-bulk-rnaseq-data`; contraindication points back to `bulk-rnaseq-gsea` for gene/MSigDB
  GSEA.
- Frontmatter follows the LIVE `_TEMPLATE` flat-metadata convention (NOT the unadopted ADR-001
  `metadata.sciagent.*` / ADR-003 orchestrator-atomic proposals): top-level keys only
  `name/description/license/metadata/...`; under `metadata:` use `scope: implementation`,
  `version: 0.1.0`, `category: analysis`, `tier: standard`, `tags: [pathway]` (exists in tags.yaml),
  `requires: []`, `complementary-skills`, `contraindications`. Description must contain no angle
  brackets and quote any raw colon. Validate with
  `python skills/skill-creator/scripts/quick_validate.py skills/te-geneset-gsea/`.

**W1 — 14839 lacks the TE-RNAseq-toolkit submodule (DEPENDENCY, OUT OF SCOPE).**
The 14839 `.gitmodules` declares only `RNAseq-toolkit` and `SciAgent-toolkit`. The TE path cannot run
inside 14839 until TE-RNAseq-toolkit is wired in at `01_modules/TE-RNAseq-toolkit`. That is an
analysis-project change, not a SciAgent-toolkit change — recorded here as a dependency only.
**Testability caveat:** no step S7-S11 executes the TE R/ factories — they edit docs and validate only
via `quick_validate.py` (a frontmatter/markdown linter). So within this plan's own steps the TE
recipe's correctness rests on CODE REVIEW against the verified `v0.1.0` factory signatures, not on an
end-to-end run. An optional, out-of-scope smoke test (record only) would source the `v0.1.0` `R/`
factories from a resolving cwd and run
`build_te_annotation(c("HAL1:L1:LINE","ENSMUSG00000000001"))` + a tiny `create_combined_dge(...)` /
`create_te_genesets(te_ann, level="family")` / `export_te_genesets_gmt(...)` round-trip to demonstrate
the repoint runs. That requires W1 (the TE submodule) and is therefore deferred.

---

## 4. Ordered, non-overlapping, sequential step plan

> One row per step. No two steps edit the same file unless one `dependsOn` the other.
> Order: RNAseq-toolkit code (G-fixes) first, then the SciAgent-toolkit doc repoint (T/G recipe),
> then the new skill scaffold last. BLOCKED steps tagged.
> **Step granularity & overlap, stated honestly:** S1, S2a, S2b, S3 all edit the SAME file
> (`io_helpers.R`); they are serialized into a `dependsOn` chain ONLY to satisfy the non-overlap rule
> (file contention), NOT because of any logical dependency — they COMMUTE, a re-orderer may sequence
> them in any order. S4 (`annotate_genes.R`) and S5 (new `provenance.R`) are the genuinely
> independent files (different files, no shared state). S6 and S7 are NOT byte-sized — they are
> **composite reference-doc rewrites** (each touches one markdown file but rewrites several recipe
> blocks); they are labelled "(composite)" and their verifies pair `quick_validate.py`
> (frontmatter/markdown shape ONLY) with authoritative `grep` content assertions.
> All paths absolute. `RT` = `/scratch/current/antonz/projects/14839-DM-cGAS/01_modules/RNAseq-toolkit`.
> `SA` = `/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit`.
> `TE` = `/data1/users/antonz/pipeline/TE-RNAseq-toolkit`.

| id | repo+abs path | files (abs paths) | change (precise, implementable) | verify (exact command/check) | dependsOn | commit msg | contextForNext (sliding-window handoff) |
|----|---------------|-------------------|---------------------------------|------------------------------|-----------|------------|------------------------------------------|
| S1 | RNAseq-toolkit `RT` (branch `dev`) | `RT/scripts/General/io_helpers.R` | In `read_counts_matrix()` add a Salmon gene-level branch: when `all(c("gene_id","gene_name") %in% names(dt))`, set rownames from `gene_id`, samples = `setdiff(names, c("gene_id","gene_name"))`, and attach `attr(mat,"input_gene_name")` = first `gene_name` per stripped stable id (`NA` for featureCounts/generic). Do NOT touch `read_metadata`/`write_annotated_matrix` yet. | `Rscript -e 'source("RT/scripts/General/io_helpers.R"); cat("ok\n")'` parses without error; manual read confirms new branch + `attr` set. | — | `RNAseq-toolkit: read_counts_matrix detects Salmon gene-level format + captures input_gene_name` | io_helpers.R now exposes input_gene_name via attr. G3 partly done; S5 threads it into annotate_genes. Dependency to S2a/S2b/S3 is FILE-CONTENTION only (same file), NOT semantic — they commute. read_metadata/write_annotated_matrix still 13036-hardcoded (S2a/S2b). |
| S2a | RNAseq-toolkit `RT` | `RT/scripts/General/io_helpers.R` | `read_metadata` only: rename arg to `read_metadata(fp, sample_col_candidates=c("Sample_ID","Sample ID"), required_cols=character(0))`; keep Sample-ID resolution; REMOVE the hardcoded 8-col `needed` `stop()`; only `stop()` if opt-in `required_cols` missing. Accept CSV as well as `.xlsx` (dispatch on file extension: `fread` for `.csv`, `read_xlsx` for `.xlsx`). | `Rscript -e 'source("RT/scripts/General/io_helpers.R"); print(names(formals(read_metadata)))'` lists the 3 args; functional check: a 2-col CSV with only `Sample_ID`+one factor loads WITHOUT error (old 8-col `stop()` gone). grep shows no hardcoded 8-col `needed` vector. | S1 | `RNAseq-toolkit: read_metadata drops 13036 8-col stop, adds required_cols opt-in + CSV/xlsx dispatch` | read_metadata now design-agnostic (G4 part 1). FILE-CONTENTION dep on S1 only. write_annotated_matrix still hardcoded (S2b). |
| S2b | RNAseq-toolkit `RT` | `RT/scripts/General/io_helpers.R` | `write_annotated_matrix` only: `write_annotated_matrix(mat, md, add_cols, outfile, factor_cols=NULL)`: replace hardcoded `top_map` with a loop over `factor_cols` — each factor name in the FIRST annotation column, SECOND annotation col blank (`""`), all other annot cols blank for top rows, values across sample columns (reference line 277); `factor_cols=NULL` → no top block. | `Rscript -e 'source("RT/scripts/General/io_helpers.R"); cat("ok\n")'` parses; functional check: `write_annotated_matrix(..., factor_cols=NULL)` writes header `c(annot_cols, samples)` with NO top rows; `factor_cols="grp"` adds exactly 1 top row whose 2nd annot cell is `""`. grep shows no `"Treatment 1"`/`top_map` remains. | S2a | `RNAseq-toolkit: write_annotated_matrix generic factor_cols loop (drop 13036 top_map hardcoding)` | io_helpers writer now reproduces reference top-block layout (G4 part 2). FILE-CONTENTION dep on S2a only. aggregate_duplicate_ids still missing (S3). |
| S3 | RNAseq-toolkit `RT` | `RT/scripts/General/io_helpers.R` | Add `aggregate_duplicate_ids(mat, ids=rownames(mat))` SUM-collapse (split row idx by `ids`, `colSums` per group, rebuild matrix with unique ids as rownames; no dups → just set rownames). Output rownames = `sort(unique(ids))` (document the reorder). MUST preserve `attr(mat,"input_gene_name")`: re-subset it to the collapsed unique-id order (first non-NA per group). Mirrors reference `agg_rows_by_id()`. | `Rscript -e 'source("RT/scripts/General/io_helpers.R"); m<-matrix(1:8,2,4,dimnames=list(c("A","A"),NULL)); attr(m,"input_gene_name")<-c("g","g"); r<-aggregate_duplicate_ids(m,c("A","A")); print(c(dim(r), identical(rownames(r),"A"), attr(r,"input_gene_name")))'` prints `1 4 TRUE g` (dims `1 4`, rowname preserved, attr survives). | S2b | `RNAseq-toolkit: add aggregate_duplicate_ids() SUM-collapse + preserve input_gene_name attr (fixes gene-annotation.md citation)` | G5 closed; doc citation valid. io_helpers.R done. Output rownames are alphabetical-by-id → S6 recipe MUST re-match annotation to rownames(mat). Gene annotation still lacks mgi_symbol/biomart pin (S4). |
| S4 | RNAseq-toolkit `RT` | `RT/scripts/General/annotate_genes.R` | Add args `input_gene_name=NULL, biomart_version=NULL, biomart_host=NULL`. Thread `version=`/`host=` into `useEnsembl()` (pass only when non-NULL). The returned tibble currently has `Symbol, Ensembl, ENTREZID, gene_biotype` (capital S/E preserved — do NOT rename here). ADD `mgi_symbol` (from biomaRt when it runs, else `mgi_symbol=Symbol`) and `input_gene_name` (joined by stripped id; `NA` when arg NULL); the final `select` must KEEP `mgi_symbol`/`input_gene_name` (current code drops `mgi_symbol`). Keep biomaRt inside `try()`. The 6-output set is exactly `Symbol, Ensembl, mgi_symbol, ENTREZID, gene_biotype, input_gene_name`. | `Rscript -e 'source("RT/scripts/General/annotate_genes.R"); print(names(formals(annotate_genes_from_ensembl)))'` lists the new 3 args; run on a 2-id stub and assert `setequal(names(out), c("Symbol","Ensembl","mgi_symbol","ENTREZID","gene_biotype","input_gene_name"))`. | S3 | `RNAseq-toolkit: annotate_genes pins Ensembl version + emits mgi_symbol/input_gene_name (keeps Symbol/Ensembl casing)` | G3+G6(biomart pin) done. annotate output = the 6 gene_index source cols, named `Symbol`/`Ensembl` (NOT yet renamed to SYMBOL/ensembl_gene_id — S5's writer does the rename). Independent file (annotate_genes.R) — no contention with S5's new provenance.R. |
| S5 | RNAseq-toolkit `RT` | `RT/scripts/General/provenance.R` (NEW) | Create `provenance.R` (G1's `write_gene_index` lives HERE, not io_helpers.R). `write_gene_index(ann_df, outfile)` RENAMES `Ensembl→ensembl_gene_id`, `Symbol→SYMBOL`, then carries `mgi_symbol, ENTREZID, gene_biotype, input_gene_name`, writing exactly `ensembl_gene_id,SYMBOL,mgi_symbol,ENTREZID,gene_biotype,input_gene_name` via `data.table::fwrite`. `write_session_provenance(outfile, genome_build=NULL, ensembl_version=NULL)` writes `sessionInfo()` + build/release lines AND, when a biomaRt session is available, the resolved archive release (not just the package version). New file — no overlap with S1-S4. | `Rscript -e 'source("RT/scripts/General/provenance.R"); df<-data.frame(Symbol="A",Ensembl="ENSMUSG1",mgi_symbol="A",ENTREZID="1",gene_biotype="x",input_gene_name=NA); write_gene_index(df,tf<-tempfile()); print(strsplit(readLines(tf,1),"\t")[[1]])'` prints exactly `ensembl_gene_id SYMBOL mgi_symbol ENTREZID gene_biotype input_gene_name`; `exists("write_session_provenance")` TRUE. | S4 | `RNAseq-toolkit: add provenance.R (write_gene_index G1 rename+project + write_session_provenance G6)` | G1+G6 closed. All gene-path helpers exist. write_gene_index does the Symbol→SYMBOL / Ensembl→ensembl_gene_id rename so the dictionary header matches spec. SciAgent doc repoint can now cite real functions (S6+). |
| S6 (composite) | SciAgent-toolkit `SA` (branch `dev`) | `SA/skills/annotate-bulk-rnaseq-data/references/gene-annotation.md` | Composite recipe rewrite. (a) source `io_helpers.R`/`annotate_genes.R`/`dge_helpers.R`/`provenance.R`; (b) read `project.id`+`project.genome_build` from config; for `reference.ensembl_version`: EITHER add it to the config template here (then pin) OR state it is ABSENT so `biomart_version` stays NULL (floating) with provenance-only repro (Risk 4) — do not silently read a missing key. (c) build `<id>_DGEList.rds`, `<id>_gene_index.tsv`, `<id>_counts_transposed.tsv`, AND `<id>_counts_wide_by_sample.tsv` (table A — both signature matrices, reference lines 233-239 + 262-299). (d) call `aggregate_duplicate_ids` then RE-ORDER annotation to `rownames(mat)` via `match()` with a `stopifnot(identical(ann$ensembl_gene_id, rownames(mat)))` alignment check (reference 219/267) BEFORE DGEList/writers. (e) call `annotate_genes_from_ensembl(..., input_gene_name=attr(mat,'input_gene_name'), biomart_version=)`, `write_annotated_matrix(..., factor_cols=<derived in recipe>)`, `write_gene_index`, `write_session_provenance`. Pin `v2.0.0`→`v0.2.0`. | `python SA/skills/skill-creator/scripts/quick_validate.py SA/skills/annotate-bulk-rnaseq-data/` passes (frontmatter shape ONLY). Authoritative content greps: `_gene_index.tsv`, `_counts_wide_by_sample.tsv`, `aggregate_duplicate_ids`, `match(`, `stopifnot`, `project.id`, `v0.2.0` all present; `v2.0.0` absent. | S5 | `annotate-bulk-rnaseq-data: gene recipe — config names, both signature matrices, alignment re-match, gene_index, provenance, real v0.2.0 (G1-G6)` | gene-annotation.md reproduces reference gene path generically incl. the re-match after aggregate. te-annotation.md still cites thin path + flat-TSV combine + phantom TE pins (S7). |
| S7 (composite) | SciAgent-toolkit `SA` | `SA/skills/annotate-bulk-rnaseq-data/references/te-annotation.md` | Composite recipe rewrite. (T1) Replace the THIN `source(".../scripts/te_utils.R")` line and the hand-rolled `rbind`/`Symbol=Symbol` recipe; instead source the rich `R/` in order: `R/te_utils.R` (constants `TE_CYTOPLASMIC_CLASSES`/`TE_NUCLEAR_CLASSES`, `parse_te_id`, `build_te_annotation`), `R/validate_te_input.R` (REQUIRED — `create_combined_dge` defaults `validate=TRUE` → calls `validate_combined_input`), `R/create_combined_dge.R`, `R/create_te_genesets.R`. (T2) `create_combined_dge(gene_counts, te_counts, samples, gene_annotation=, te_annotation=build_te_annotation(rownames(mat_te)))` → `<id>_combined_DGEList.rds`. (T3) Family GMT: `gs_fam<-create_te_genesets(te_ann, level="family"); export_te_genesets_gmt(gs_fam, "<id>_TE_family.gmt")` (note `te_ann` is a REQUIRED positional; `export_te_genesets_gmt` takes the LIST). Class GMT: do NOT call `create_te_genesets(..., level="class")` (its `R/create_te_genesets.R:109` internal `source()` hard-errors from a non-resolving cwd — T1b); build class T2G/T2N in the recipe from `class`/`subfamily` + the sourced constants, then `export_te_genesets_gmt(list(T2G=,T2N=), "<id>_TE_class.gmt")`. Relaxed `min_size` note. Document gene slot=`subfamily` (col `gene_symbol`), `feature_type` from `create_combined_dge`, `is_te`/`replication_type` from `build_te_annotation`. (T4) Reconcile EVERY TE pin to `TE-RNAseq-toolkit v0.1.0`. Keep graded caveats block. | `python .../quick_validate.py SA/skills/annotate-bulk-rnaseq-data/` passes (shape ONLY). Authoritative greps: `create_combined_dge`, `validate_te_input`, `export_te_genesets_gmt`, `replication_type`, `TE-RNAseq-toolkit v0.1.0` present; `v2.0.1`, `v2.0.0`, and the thin `scripts/te_utils.R` source line all ABSENT. | S6 | `annotate-bulk-rnaseq-data: TE recipe repointed to rich R/ v0.1.0 factories (T1-T3, class-GMT via recipe-built T2N), reconcile pins to v0.1.0 (T4)` | te-annotation.md uses rich v0.1.0 factories; family GMT via factory, class GMT recipe-built (avoids T1b hard-error). Reconciled pin string = `TE-RNAseq-toolkit v0.1.0` (S8 greps for this literal). SKILL.md Resources line still has phantom pins + no te-geneset-gsea edge (S8). |
| S8 | SciAgent-toolkit `SA` | `SA/skills/annotate-bulk-rnaseq-data/SKILL.md` | Update Resources lines 103-104: `RNAseq-toolkit v0.2.0`; `TE-RNAseq-toolkit v0.1.0` (the exact reconciled literal from S7) and point TE helpers at the rich `R/` factory files. Add `te-geneset-gsea` to the Complementary Skills table as the TE-GMT→GSEA next step and to the handoff chain. Bump `metadata.version` 1.1.1→1.2.0 + add CHANGELOG comment lines. | `python .../quick_validate.py SA/skills/annotate-bulk-rnaseq-data/` passes (frontmatter shape ONLY — does NOT check version/tags). Authoritative greps: `RNAseq-toolkit v0.2.0`, `TE-RNAseq-toolkit v0.1.0`, `te-geneset-gsea`, `version: 1.2.0` present; phantom `v2.0.0`/`v2.0.1` ABSENT (the reconciled string is `v0.1.0`, so a blanket "no v2.0.0" grep is safe and does not collide with the chosen pin). | S7 | `annotate-bulk-rnaseq-data: SKILL.md real pins (RNAseq v0.2.0 / TE v0.1.0) + te-geneset-gsea edge, bump 1.2.0` | annotate skill fully consistent on real tags. Next: scaffold the new te-geneset-gsea skill (S9-S11). |
| S9 | SciAgent-toolkit `SA` | `SA/skills/te-geneset-gsea/SKILL.md` (NEW) | Copy `_TEMPLATE/SKILL.md` → new dir (create `references/` dir but do NOT add `.gitkeep` — S10 immediately adds a real reference file, so the dir is never empty-tracked). Set `name: te-geneset-gsea`, description (TE family/class GMT to TE-aware GSEA; "Use when running GSEA on TE subfamily ranks"; "For gene/MSigDB GSEA use bulk-rnaseq-gsea"; no angle brackets, quote any raw colon), `scope: implementation`, `version: 0.1.0`, `category: analysis`, `tier: standard`, `tags:[pathway]`, `requires:[]`, complementary `[annotate-bulk-rnaseq-data, bulk-rnaseq-gsea]`, contraindication → `bulk-rnaseq-gsea`. Body: Overview/Decision Tree/Quick Start+Verify/Progressive Depth/Verification/Pitfalls (ID space=subfamily, overlap<30% trap, @geneSets fix, relaxed min_size)/Complementary/Resources. | `python .../quick_validate.py SA/skills/te-geneset-gsea/` passes (it only gates allowed top-level keys + kebab name + description shape; it does NOT validate tags/category/scope/version — so ALSO grep-assert `name: te-geneset-gsea`, `tags`, `pathway`, `version: 0.1.0` manually, and confirm `pathway` is in `SA/tags.yaml`). | S8 | `te-geneset-gsea: scaffold new skill (TE family/class GMT GSEA) — N1` | New skill SKILL.md exists + validates (shape). `references/` exists, no `.gitkeep`. Next: write the reference doc (S10). |
| S10 | SciAgent-toolkit `SA` | `SA/skills/te-geneset-gsea/references/te-gsea.md` (NEW) | Write the load-on-demand reference (<=300 lines): how the family/class GMT (from annotate S7) + a TE-subfamily-ranked DE list feed `clusterProfiler::GSEA(TERM2GENE,TERM2NAME)`; reuse bulk-rnaseq-gsea machinery at the ALGORITHM level but state the column-name contract explicitly — TE GMTs use `gs_name`/`gene_symbol`, so the @geneSets fix is `gsea_result@geneSets <- split(T2G$gene_symbol, T2G$gs_name)` and the overlap check is `sum(db$T2G$gene_symbol %in% names(ranked))` (do NOT copy the bulk `T2G$gene`/`T2G$term` snippet verbatim — empty-set trap); relaxed `min_size` for ~dozens of sets; SCREAMING_SNAKE prefix naming (`TE_FAMILY_*`/`TE_CLASS_*`). No `.gitkeep` to remove (S9 didn't create one). | `python .../quick_validate.py SA/skills/te-geneset-gsea/` passes; `wc -l te-gsea.md` <=300; grep shows `@geneSets`, `split(T2G$gene_symbol, T2G$gs_name)`, `gene_symbol %in% names(ranked)`, `min_size`; confirm NO `references/.gitkeep` exists. | S9 | `te-geneset-gsea: add references/te-gsea.md (GMT→clusterProfiler GSEA, TE gs_name/gene_symbol contract)` | te-geneset-gsea complete with correct TE column-name contract. Last: cross-link bulk-rnaseq-gsea so the TE handoff is discoverable (S11). |
| S11 | SciAgent-toolkit `SA` | `SA/skills/bulk-rnaseq-gsea/SKILL.md` | Add a one-line contraindication/complementary pointer: "For TE family/class geneset GSEA use te-geneset-gsea" (Complementary Skills table + a contraindication entry). Do NOT change custom-db.md. Bump `metadata.version` 1.0.0→1.0.1 + CHANGELOG comment. | `python SA/skills/skill-creator/scripts/quick_validate.py SA/skills/bulk-rnaseq-gsea/` passes; grep shows `te-geneset-gsea`, version `1.0.1`. | S10 | `bulk-rnaseq-gsea: cross-link te-geneset-gsea for TE-set GSEA handoff (N1)` | All skills cross-consistent. T1b (TE R/ self-locating source) + W1 (14839 submodule) remain as out-of-scope follow-ups (see doc §5/§6). |
| S12 | TE-RNAseq-toolkit `TE` (HEAD `feat/strand-split-qc`) — **BLOCKED / OUT-OF-SCOPE FOLLOW-UP, do NOT run this phase** | `TE/R/create_te_genesets.R` | (FOLLOW-UP) Make lines 109,175 self-locating (e.g. `here::here`/`system.file`) instead of hardcoded `source("01_modules/TE-RNAseq-toolkit/R/te_utils.R")`, so `create_te_genesets(level="class")` runs from any cwd. BLOCKED: no analysis project wires this checkout, and it sits on a feature branch with uncommitted doc edits — editing changes nothing reproducible. NOTE: S7 does NOT depend on this — S7 sidesteps the hard-error by recipe-building the class T2N (T1b), so this is a true follow-up, not a prerequisite. Unblocks with W1. | n/a — not executed this phase. If executed later: `Rscript -e 'source("TE/R/create_te_genesets.R"); create_te_genesets(build_te_annotation("HAL1:L1:LINE"), level="class")'` from an arbitrary cwd succeeds. | W1 (out of scope) | `TE-RNAseq-toolkit: make create_te_genesets source() self-locating (follow-up)` | Recorded only; no dependency on this phase. S7 avoids the class-level `source()` hard-error by recipe-building the class GMT, NOT by pre-sourcing (pre-sourcing does not neutralize a path-based `source()`). |

---

## 5. Risks

1. **RNAseq-toolkit tag drift.** S1–S5 advance `dev` past the `v0.2.0` tag. The 14839 submodule
   gitlink still points at `v0.2.0`; six other projects share `v0.2.0` checkouts and will NOT
   auto-update. Mitigation: after S5, cut a new release tag (e.g. `v0.3.0`) in RNAseq-toolkit and do
   a coordinated submodule bump per the toolkit's `AGENTS.md` pinning procedure. The refactor doc
   does not perform the bump (project concern), but flags it.
2. **TE pin — RESOLVED (was: "no real git tag").** Earlier drafts claimed TE-RNAseq-toolkit had no
   annotated release tag. That is FALSE: tag `v0.1.0` exists on `main` and contains the cited rich `R/`
   factories byte-identical to the working tree (`git diff v0.1.0 -- R/<f>` = 0 lines). T4 therefore
   pins every TE reference to the real tag `v0.1.0`. The in-file `Version: 2.0.0` header is NOT a git
   tag and must NOT be used (it would re-introduce the phantom-pin bug T4 exists to kill). No residual
   risk beyond keeping the doc and SKILL.md citing the same `v0.1.0` literal (S7→S8 hand off it).
3. **TE-RNAseq-toolkit checkout has uncommitted doc edits** on a feature branch. Any code edit there
   (S12, deferred) must branch/coordinate to avoid mixing with in-progress work. We avoid this by
   keeping S12 out of scope. Note: the class-level GMT `source()` hard-error (T1b) is sidestepped in
   the S7 recipe by building the class T2N from sourced constants — NOT by pre-sourcing (which does
   not neutralize a path-based `source()`).
4. **biomaRt floating release + absent config key.** The `biomart_version` arg (G6) defaults to
   floating, and `reference.ensembl_version` does NOT currently exist in the 14839 config (only
   `project.{id,genome_build,species_db}`). So byte-reproducibility holds only if S6 ADDS that key (the
   preferred default path) — otherwise the run floats. Even then, `sessionInfo()` records only the
   biomaRt PACKAGE version, NOT the live remote Ensembl archive release that actually varies.
   Mitigation: `write_session_provenance` must ALSO capture the RESOLVED archive release
   (`biomaRt::listEnsemblArchives()` / the host actually used), so the varying quantity is recorded
   even when not pinned.
5. **TE GSEA ID-space mismatch.** The #1 GSEA failure mode (non-symbol rownames → 0 results) maps
   onto TEs as "ranked list keyed by Ensembl/gene symbols, GMT keyed by TE subfamily". S10 must make
   the overlap check + ID-space requirement load-bearing, and relax `min_size` (default 5 drops most
   of the ~6 TE classes).
6. **Validator scope is NARROW — it does NOT check tags/category/scope/version.** `quick_validate.py`
   only enforces: allowed top-level keys (`name, description, license, allowed-tools, metadata,
   compatibility`), kebab-case `name` (not `SKILL_IDENTIFIER`), and `description` shape (no angle
   brackets, length). It does NOT validate `tags`-in-`tags.yaml`, `category`, `scope`, or `version` at
   all (verified). Therefore the version bumps (S6/S8/S11), pin reconciliation (S7/S8), and tag
   validity are caught ONLY by the accompanying `grep` assertions — those greps are the authoritative
   content gate; `quick_validate` is only a frontmatter-shape gate. `pathway` IS in `tags.yaml`
   (verified) but the implementer must confirm this by grep, not rely on the validator. The
   ADR-001/ADR-003 proposals are NOT adopted — scaffold against the live flat `_TEMPLATE` convention.

## 6. Out of scope

- **W1 — adding `01_modules/TE-RNAseq-toolkit` as a submodule of the `14839-DM-cGAS` analysis
  project.** This is an analysis-project `.gitmodules` change, not a SciAgent-toolkit change.
  Recorded as a hard dependency: the TE path is non-runnable inside 14839 until this lands. The doc
  must not present analysis-project-private paths (`/scratch/.../14839-DM-cGAS/...`) as toolkit
  contracts — they appear here only as the editable RNAseq-toolkit submodule checkout and as the W1
  dependency note.
- **S12 — TE-RNAseq-toolkit `R/` `source()` self-location fix.** Deferred upstream follow-up
  (BLOCKED, see §4); dependsOn W1, not on any step in this phase. The S7 recipe sidesteps the
  class-level `source()` hard-error by recipe-building the class T2N from sourced constants (NOT by
  pre-sourcing, which does not neutralize a path-based `source()`), so S12 is a true follow-up rather
  than a prerequisite of S7.
- **The RNAseq-toolkit release re-tag + submodule bump** (Risk 1 mitigation) — project/release
  operation, not part of the doc-and-helpers refactor.
- **Any change to `bulk-rnaseq-gsea/references/custom-db.md`** — left untouched; only the SKILL.md
  cross-link (S11) is in scope.
