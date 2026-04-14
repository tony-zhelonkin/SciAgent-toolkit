---
name: bulk-rnaseq-pathway-explorer
description: "pathway-explorer — generate standalone interactive HTML dashboards from unified GSEA / TF / PROGENy / TE master tables, embedding pathways on a UMAP-of-gene-sets scatter with Jaccard/Overlap neighbor edges, per-pathway running-sum plots, and database/entity filters. Use when turning master_unified.csv (or legacy master_gsea_table.csv + master_tf_activities.csv + master_progeny_activities.csv) into a shareable .html for pathway exploration, or when a user asks for an interactive pathway explorer, pathway scatter, or cross-database dashboard. A contrast switcher (collapsible toggle, single HTML across all contrasts) is planned — today the tool ships per-contrast HTMLs plus an index.html. For running GSEA itself use bulk-rnaseq-gsea-msigdb or bulk-rnaseq-gsea-custom-db. For building the master tables this skill consumes, use bulk-rnaseq-gsea-master-tables. For static publication figures (dotplot, barplot, running-sum PDFs) use bulk-rnaseq-gsea-visualization."
license: MIT
metadata:
  skill-author: SciAgent-toolkit
  last-reviewed: 2026-04-14
  version: 0.1.0
  upstream-docs: https://github.com/tony-zhelonkin/pathway-explorer
  category: workflow
  tier: standard
  tags:
    - gsea
    - pathway-explorer
    - interactive-dashboard
    - html
    - plotly
    - umap
    - jaccard
    - transcription-factors
    - progeny
    - transposable-elements
    - bulk-rnaseq
  complementary-skills:
    - bulk-rnaseq-gsea-msigdb
    - bulk-rnaseq-gsea-custom-db
    - bulk-rnaseq-gsea-master-tables
    - bulk-rnaseq-gsea-visualization
  contraindications:
    - "Do not use for running GSEA itself. Use bulk-rnaseq-gsea-msigdb or bulk-rnaseq-gsea-custom-db instead."
    - "Do not use for assembling the master tables this skill consumes. Use bulk-rnaseq-gsea-master-tables instead."
    - "Do not use for static publication figures (PDF/PNG dotplots, barplots, running-sum plots). Use bulk-rnaseq-gsea-visualization instead."
    - "Do not use on raw gseaResult RDS checkpoints directly. The tool reads CSVs, not R objects — normalize first via bulk-rnaseq-gsea-master-tables."
---

# Pathway Explorer — Interactive Cross-Entity Pathway Dashboards

## Overview

`pathway-explorer` is a Python CLI (module `pathway_explorer`, v2.0.0, MIT) that turns unified GSEA / TF / PROGENy / TE result tables into self-contained interactive HTML dashboards. Each dashboard shows all pathways as points on a 2D embedding of their leading-edge gene sets, with database and entity-type filters, FDR/NES sliders, gene-set neighbor edges, a gene table, and per-pathway running-sum plots.

The tool sits downstream of `bulk-rnaseq-gsea-master-tables` (which produces the CSV inputs) and alongside `bulk-rnaseq-gsea-visualization` (which produces static R/ggplot PDFs). It is intentionally decoupled from R: the master CSV schema is the one bridge.

> **Planned change — contrast switcher.** Today the tool emits one HTML per contrast plus an `index.html` landing page. We intend to migrate to a **single HTML that embeds all contrasts and exposes a collapsible toggle** to switch between them in-place, so reviewers can compare contrasts without navigating to a new file. When writing or extending this skill, prefer changes that keep the JSON payload contrast-aware (nested by contrast, not pre-filtered) so the future switcher is a UI-only addition.

**When to use this skill:**
- Producing an interactive pathway HTML for a project that already has `master_unified.csv` (or the legacy triplet: `master_gsea_table.csv` + `master_tf_activities.csv` + `master_progeny_activities.csv`).
- Exploring pathway–TF–PROGENy–TE relationships in one view (different shapes per entity type, shared color scale on signed significance).
- Sharing a pathway view with collaborators who do not have R / clusterProfiler installed (output is a standalone `.html`).
- Debugging when a dashboard is empty, mis-colored, or missing a database — nearly always a master-table schema issue.

**When NOT to use this skill:**
- Running GSEA itself → use `bulk-rnaseq-gsea-msigdb` or `bulk-rnaseq-gsea-custom-db`.
- Building or appending to `master_*.csv` → use `bulk-rnaseq-gsea-master-tables`.
- Static publication figures (dotplot, barplot, running-sum PDFs) → use `bulk-rnaseq-gsea-visualization`.
- Metabolic network (atom-transition) visualization → use `gatom-metabolomic-predictions`.

---

## Decision Tree

```
Need an interactive pathway view?
│
├─ Do you already have master_unified.csv (or master_gsea_table.csv)?
│   ├─ No  → run bulk-rnaseq-gsea-master-tables first, then come back
│   └─ Yes → continue
│
├─ One contrast or many?
│   ├─ One        → `pathway-explorer --contrast <name>`
│   ├─ All        → `pathway-explorer --all` (emits per-contrast HTMLs + index.html)
│   └─ Future     → one HTML with a collapsible contrast toggle (see "Planned Work")
│
├─ Want to restrict entity types or TE level?
│   ├─ Only pathways + TFs  → `--entity-types Pathway TF`
│   └─ TE class vs family   → `--te-level class|family`
│
└─ Custom data path or output location?
    → `--data <csv> --output <html>`
```

---

## Quick Start

Run from a project root that has `03_results/tables/master_unified.csv`:

```bash
# Single-contrast dashboard
pathway-explorer --contrast Sema_WL_vs_Base

# All contrasts (writes per-contrast HTMLs + index.html)
pathway-explorer --all

# Explicit I/O
python -m pathway_explorer \
    --data 03_results/tables/master_unified.csv \
    --output 03_results/interactive/my_dashboard.html \
    --contrast Sema_WL_vs_Base
```

Outputs land in `03_results/interactive/pathway_explorer_<contrast>.html`.

**Verify it worked:**

- The file is non-empty (`> 100 KB` typically — it embeds Plotly and all pathway JSON).
- Opening in a browser shows a scatter plot, left-side database/FDR filters, and non-empty entity-type counts in the sidebar.
- Clicking a point opens a gene table and a running-sum plot (the running sum requires `master_de_table.csv` next to the input; without it the plot is empty but the rest of the dashboard still works).

---

## Progressive Depth

### Basic Usage — What Artifacts to Ingest, and in What Format

The tool reads CSVs only; no RDS, no Parquet, no pickled R objects. Put them in `03_results/tables/` of the project (the path is discovered by walking up from the CLI's CWD until a `03_results/tables` directory is found — see `config.py::get_project_root`).

**Primary input (preferred):** `03_results/tables/master_unified.csv`

Columns read by the loader (`data_loader.py::load_gsea_data`):

| Column | Type | Required? | Purpose |
|---|---|---|---|
| `pathway_id` | str | yes | Stable ID (MSigDB / CollecTRI / PROGENy / TE identifier). Used as neighbor key. |
| `pathway_name` | str | yes | Display name (cleaned to Title Case, truncated to 60 chars). |
| `database` | str | yes | Source (`Hallmark`, `KEGG`, `Reactome`, `GO_BP`, `MitoPathways`, `MitoXplorer`, `CollecTRI`, `PROGENy`, `TE_Class`, `TE_Family`, `TransportDB`, `GATOM`, …). Drives `DB_COLORS`. |
| `entity_type` | str | yes (or derivable) | One of `Pathway`, `TF`, `PROGENy`, `TE`. If absent, inferred: `CollecTRI → TF`, `PROGENy → PROGENy`, `TE_* → TE`, else `Pathway`. Drives `ENTITY_SHAPES` (circle / diamond / square / triangle-up). |
| `nes` / `NES` | float | yes | Normalized enrichment score. Either casing is accepted; the loader normalizes to lowercase `nes`. |
| `padj` / `adj.P.Val` | float | yes | FDR-corrected p-value. Either naming is accepted; normalized to `padj`. Clipped to `1e-50` before `-log10`. |
| `pvalue` | float | yes | Raw p-value (shown in tooltip). |
| `set_size` | int | yes | Total gene-set size (for the size encoding). |
| `leading_edge_size` | int | no (derived from `genes`) | Count of leading-edge genes. |
| `core_enrichment` | str | yes | Slash-separated list of leading-edge gene symbols (`GENE1/GENE2/…`). This is what seeds the similarity matrix. |
| `contrast` | str | required for `--all` / contrast filtering | Contrast label; `pathway_id` values may repeat across contrasts. |
| `direction` | str | no | `Up` / `Down`; else computed from `sign(nes)`. |
| `signed_sig` | float | no (derived) | Internal score `-log10(padj) * sign(nes)`, capped at ±50. Set by `standardize_scores` if missing. |

**Legacy three-file mode (still supported):** if `master_unified.csv` is missing, the loader falls back to:

- `master_gsea_table.csv` (MSigDB / custom pathway GSEA; schema above minus `entity_type`).
- `master_tf_activities.csv` (CollecTRI TF activities in the GSEA column schema).
- `master_progeny_activities.csv` (PROGENy pathways in the GSEA column schema).

The loader concatenates all three, reclassifies a generic `Mitochondria` database into `MitoPathways` / `MitoXplorer` based on `pathway_id` prefix, and adds `entity_type` from the database name. Prefer unified — the legacy branch exists for backward compatibility only.

**Optional input for running-sum plots:** `03_results/tables/master_de_table.csv` with columns `gene_symbol, t, logFC, adj.P.Val`. Used only to render the per-pathway running sum when a point is clicked. If absent, the dashboard still works; the running-sum panel is empty.

**Output:** one standalone HTML per contrast at `03_results/interactive/pathway_explorer_<contrast>.html`, plus `index.html` when `--all`. Plotly is loaded from CDN (`PLOTLY_JS_URL` in `config.py`); all pathway / neighbor / gene-ranking data is inlined as JSON.

### Intermediate Usage — Pipeline and Tunables

End-to-end flow (`main.py::generate_dashboard`):

1. Load master table (unified or legacy triplet).
2. Filter by `contrast` (if `--contrast` / `--all`).
3. Filter by `entity_types` (if `--entity-types`) and TE level (`--te-level family|class`).
4. `standardize_scores` → `signed_sig = -log10(padj) * sign(nes)`, clipped to ±50.
5. `compute_hybrid_similarity` over leading-edge gene sets:
   - Same-type pairs (Pathway–Pathway, TF–TF, PROGENy–PROGENy) → **Jaccard**.
   - TF ↔ Pathway / TF ↔ PROGENy → **Overlap coefficient** (`|A∩B| / min(|A|, |B|)`) because TF regulons are much larger than Hallmark-style gene sets and Jaccard collapses.
   - PROGENy ↔ Pathway → Jaccard (both pathway-scale).
6. `extract_top_neighbors(k=5)` — keep the top-5 neighbors per pathway above `MIN_JACCARD_EDGE = 0.15` for the edge overlay.
7. `compute_embedding` — prefers UMAP, falls back to PCA, then random projection. Pick is auto (`embedding.get_best_method`).
8. `prepare_pathway_data` → list of dicts with `{id, name, database, entity_type, nes, padj, x, y, genes, neighbors, …}` serialized to JSON.
9. `generate_html` — assembles CSS + JS + JSON into a single file.

Thresholds worth knowing (in `config.py`):

- `MIN_JACCARD_EDGE = 0.15` — minimum similarity to draw an edge; raise to de-clutter, lower to show more.
- `NES_MAX = 3.5` — color-scale cap on `NES`.
- `DEFAULT_FDR_SLIDER = 0.05` — initial FDR filter in the UI slider.
- `DB_COLORS`, `ENTITY_SHAPES` — change these to re-theme the dashboard.

Install with the optional UMAP backend to avoid the PCA fallback:

```bash
pip install -e "01_modules/pathway-explorer[full]"
```

### Advanced Usage — Extending and Theming

- **New entity type.** Add to `ENTITY_TYPES` and `ENTITY_SHAPES` in `config.py`, update `_add_entity_types` if the source is the legacy triplet, and make sure the upstream master table sets `entity_type` for the new rows.
- **New database.** Append to `DB_COLORS`. The loader treats anything it sees as valid; missing colors fall through to a default.
- **Changing the similarity metric.** `compute_hybrid_similarity` picks metric per pair based on `entity_types`. If a new cross-type comparison is added, extend the dispatch there rather than forcing one global metric — Jaccard and Overlap answer different questions and the rest of the pipeline assumes the hybrid matrix.
- **Custom output path or landing page.** `main.py::generate_index_page` builds a minimal landing page; swap in your own template string if the project needs branded indexes.

---

## Planned Work — Contrast Switcher (single HTML)

The current design is one HTML per contrast plus an `index.html` with per-contrast links. The **next iteration** will emit a single HTML that:

- Embeds JSON for **all** contrasts (nested under a `contrasts` key, not flattened).
- Exposes a **collapsible sidebar toggle** (or a dropdown) to switch the active contrast in-place — no page reload.
- Keeps filters, selections, and the running-sum panel state consistent when switching, so a reviewer can compare "Sema_WL_vs_Base" vs "Ctrl_WL_vs_Base" on the same pathway without losing context.

Design hints for contributors:

- The JSON payload today is keyed flat (`pathways: [...]`). Change the producer (`html_generator.prepare_pathway_data` / `generate_html`) to emit `{contrasts: {"<name>": {pathways, metadata, neighbors}, ...}, default_contrast: "<name>"}`.
- The UI change is JS-only: on toggle, swap `pathways` and re-render the Plotly trace; the existing filter/neighbor code should keep working.
- Keep the `--contrast` and `--all` CLI flags; add a third mode (e.g. `--all-in-one`) that emits the multi-contrast HTML.
- Gene-ranking running sum is contrast-specific — make sure `master_de_table.csv` is loaded per contrast (or replace with a contrast-keyed `master_de_tables.csv`).

Until then, treat `--all` as the working mode for multi-contrast projects.

---

## Verification Checklist

After running this skill, confirm:

- [ ] **Output file exists and is non-trivial.** `ls -la 03_results/interactive/pathway_explorer_*.html` — should be hundreds of KB, not a few KB.
- [ ] **Sidebar entity-type counts are sane.** All four types present for a full project; if `TF: 0` or `PROGENy: 0`, the legacy triplet was likely loaded without the corresponding file.
- [ ] **Database legend covers what you expect.** Missing `MitoPathways` / `MitoXplorer` means the reclassifier didn't see the expected `pathway_id` prefix — check the master table.
- [ ] **Edges draw between related pathways.** Zero edges → `MIN_JACCARD_EDGE` is too strict for your leading-edge sizes, or `core_enrichment` is blank.
- [ ] **Clicking a point populates the gene table** (and the running-sum if `master_de_table.csv` is present).

---

## Common Pitfalls

### Pitfall: Empty dashboard after filtering

- **Symptom:** `ValueError: No pathways found after filtering!` or an HTML with zero points.
- **Cause:** `--contrast` value doesn't match any row in `master_unified.csv`, or the TE filter removed everything because the source only had one TE level.
- **Fix:** Run `pathway-explorer --data <csv>` with no filters once, read the "Found N contrasts" log line, and re-run with the exact contrast string. For TEs, verify `database` values (`TE_Class`, `TE_Family`) before setting `--te-level`.

### Pitfall: `NES column not found. Available: [...]`

- **Symptom:** Crash in `standardize_scores`.
- **Cause:** Master table uses neither `nes` nor `NES` — usually a toolkit-vs-project naming split. The R side of the RNAseq-toolkit writes `NES`; some project-side normalizers rename to `nes`.
- **Fix:** Rename in the master-table assembler (see `bulk-rnaseq-gsea-master-tables` pitfalls). Don't monkey-patch `data_loader.py` — the loader accepts either casing, so the real bug is usually a third column name.

### Pitfall: All points cluster in one blob

- **Symptom:** UMAP scatter is a single tight cluster; no spatial separation by database.
- **Cause 1:** `core_enrichment` has empty strings → every pathway is an empty set → similarity matrix is all zeros → embedding is meaningless.
- **Cause 2:** Fell back to `random` embedding because neither UMAP nor sklearn is installed.
- **Fix 1:** Inspect `master_unified.csv` — `core_enrichment` must be a `/`-joined list of gene symbols. If your upstream step drops this column, see `bulk-rnaseq-gsea-master-tables`.
- **Fix 2:** `pip install umap-learn scikit-learn` or install the `[full]` extra.

### Pitfall: Mitochondria database shows as one color instead of two

- **Symptom:** Dashboard lumps MitoCarta and mitoXplorer under a single `Mitochondria` color.
- **Cause:** `_reclassify_mito_databases` reads `pathway_id` prefixes (`MITOPATHWAYS_`, `MITOXPLORER_`). If the master-table assembler stripped or lower-cased those prefixes, reclassification silently fails.
- **Fix:** Preserve the prefix in the master table, or pre-split the `database` column upstream so `MitoPathways` / `MitoXplorer` arrive already-separated.

### Pitfall: Running-sum panel is always empty

- **Symptom:** Click a point, gene table loads, running-sum plot stays blank.
- **Cause:** `master_de_table.csv` is missing (the loader warns but continues). Gene rankings default to an empty DataFrame.
- **Fix:** Place `master_de_table.csv` with columns `gene_symbol, t, logFC, adj.P.Val` next to the master table. One file covers all contrasts today; for the planned contrast switcher this will need to become contrast-keyed.

### Pitfall: Dashboard size balloons past 10 MB

- **Symptom:** HTML is slow to open, GitHub won't render a preview.
- **Cause:** No upstream FDR filter and very large leading-edge gene lists get inlined as JSON.
- **Fix:** Pre-filter the master table (e.g. `padj < 0.25` and/or top-N per database) before invoking the tool, or set `MAX_PATHWAYS` in `config.py`. The UI slider then filters within the embedded set.

---

## Complementary Skills

| When you need... | Use skill | Relationship |
|---|---|---|
| Running GSEA (MSigDB H/C2/C3/C5) | `bulk-rnaseq-gsea-msigdb` | Prerequisite |
| Running GSEA with custom gene sets (MitoCarta, TransportDB, GATOM) | `bulk-rnaseq-gsea-custom-db` | Prerequisite |
| Assembling / appending `master_*.csv` tables | `bulk-rnaseq-gsea-master-tables` | Prerequisite — owns the inputs this skill reads |
| Static publication figures (dotplot, barplot, running-sum PDFs) | `bulk-rnaseq-gsea-visualization` | Alternative — static track of the same data |
| Metabolic network (atom-transition) visualization | `gatom-metabolomic-predictions` | Alternative — different entity model |

---

## Resources

- **Module source:** `01_modules/pathway-explorer/pathway_explorer/` (this project).
- **Upstream repo:** https://github.com/tony-zhelonkin/pathway-explorer
- **RNAseq-toolkit workflow docs:** `01_modules/RNAseq-toolkit/docs/WORKFLOWS.md`, `docs/GSEA-workflow/04-output-artifacts-and-visualization.md`.
- **Sibling visualization skill (R/ggplot2 track):** `bulk-rnaseq-gsea-visualization`.
- **Master-table schema:** `bulk-rnaseq-gsea-master-tables`.
