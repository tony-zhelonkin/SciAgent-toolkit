# Tracer-bullet check — comp-bio, DC_hum_verse scenario

**Date:** 2026-06-02 · **Persona:** Anton-as-working-comp-bio, mid-project on DC_hum_verse
**Toolkit:** `/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit`
**Method:** Read the real project to ground the scenario; ran all mutating sciagent ops
in a throwaway `mktemp -d` scratch project. The real project was read-only.

---

## (a) Scientific scenario adopted

Grounded in the real project at `/scratch/current/antonz/projects/DC-nexus/DC_hum_verse`
(`context.md`, `README.md`, `02_analysis/scripts/`):

> **DC_hum_verse** is a *human dendritic-cell atlas* project. Two scRNA-seq atlases
> (DC-VERSE Seurat RDS; Shi et al. PanCancer h5ad), integrated with **scVI**, to ask
> whether **iron/transferrin (TFRC/CD71) metabolism** defines functionally distinct DC
> states across normal vs tumor tissue. Real scripts show the live workflow:
> `S01_subset_cdc1.py → S02_harmonize → S03_integrate (scVI) → S04_compute_iron`.

So my concrete need as the scientist: **scRNA QC → batch integration (scVI/harmony) →
cell-type/subset annotation → iron/transferrin signature scoring → deliver Seurat/Loupe
objects to the PI.** Pure transcriptomics + pathway scoring. No ATAC, no chromatin.

The question I brought to sciagent cold: *"Which role do I activate for this, and which
skills am I missing?"*

## (b) Self-directed path (actual commands + output)

1. **`sciagent --help`** — clean, complete usage. No hint on role→science mapping (by design).
2. **`sciagent list roles`** — the single most useful screen. Five roles, each with a
   one-line scientific blurb + skill/agent/cmd counts:
   ```
   base               skills:36 ... Default bioinformatics analysis role with full agent suite
   pathway-signature  skills:8  ... Pathway/TF/signature functional interpretation — GSEA + decoupleR + CoReSh
   scatac-regulatory  skills:13 ... scATAC-seq regulatory analysis — CREscendo, ChromVAR, ...
   architect          skills:2  agents:13 cmds:19  Architecture-first dev harness
   software-tool      skills:1  ... Standalone software library / CLI development
   ```
   From this alone I correctly concluded: my project is `base` (it lists scvi/scanpy/qc),
   and `pathway-signature` *sounds* like my iron-scoring step.
3. **`sciagent list role base` / `list role pathway-signature`** — drilled in. Discovered
   **pathway-signature's 8 skills are a strict SUBSET of base's 36** (anndata, scanpy,
   conversion, gsea, activity-inference, pathway-explorer, coresh, gatom — all already in base).
4. **`sciagent new project <scratch> --type analysis --title DC_hum_verse --species human`**
   — scaffolded a clean `00_data/01_modules/02_analysis/03_results/docs/_internal` tree with
   `.gitkeep`s, `analysis_config.yaml`, `scientific-context.md` (species pre-filled), and a
   well-routed `docs/_internal/{reasoning,sessions,research}` with naming conventions. Good
   "Next steps" hint pointing at scbio-docker's `init-container.sh`.
5. **`sciagent activate base`** → `skills:39 agents:7 commands:1`. (36 catalog + 3 pulled
   `via requires:` — tobias/hint-atac/signac-footprint behind tf-footprint.)
6. **`sciagent activate base pathway-signature`** → **still 39 skills.** Overlay added zero.
7. **`sciagent status`** — excellent transparency: shows provenance per skill, `(shadows base)`,
   and `inherited via requires:`. Confirmed the overlay is pure redundancy (everything is
   "pathway-signature (shadows base)").
8. **Perks-left-out view** — diffed catalog vs active:
   ```
   comm -23 all_skills active_skills  →
   cellxgene-census-annotation, mllmcelltype-consensus-annotation,
   harmonypy-batch-integration, louper-seurat-conversion, scired,
   seurat-* , single-cell-vector-search, mofa-*, ...
   ```
   The annotation + harmony + Loupe-export skills my workflow needs are **NOT** in base
   or pathway-signature.
9. **`sciagent inject mllmcelltype-consensus-annotation` / `inject harmonypy-batch-integration`**
   — both injected cleanly into the overlay.
10. **`sciagent inject --tag annotation`** — worked, idempotent (`already injected ... nothing
    to inject`), but it *also* pulled ATAC annotation skills (`scembed-atac-annotation`) that
    don't fit a scRNA project.
11. **`sciagent list tags`** → `unknown category 'tags'`. The tag vocabulary (rich, with
    descriptions in `tags.yaml`) is **not discoverable from the CLI**.
12. **`sciagent validate`** → passed (two benign name-overlap warnings for architect/treemap).
13. **`sciagent deactivate`** → clean teardown: `.claude/skills` gone, AGENTS.md managed block
    stripped (`grep -c SCIAGENT:ROLES` → 0). Lifecycle is sound.

**Source-code fallbacks I was forced into** (the thing to hunt): I had to read
`tags.yaml` to learn the tag vocabulary, and I had to read role YAML / run `list role X`
twice to discover that pathway-signature ⊂ base. The CLI surface alone did not tell me either.

## (c) FINDINGS

1. **[CONFUSING] `pathway-signature` is a strict subset of `base`, so stacking it is a no-op.**
   `activate base` = 39 skills; `activate base pathway-signature` = 39 skills. Every overlay
   skill reports `(shadows base)`. For a scientist who reads "I'm doing pathway scoring → pick
   pathway-signature," the natural move (overlay it on base) buys *nothing* and the status
   screen is wall-to-wall "shadows base" noise. Either pathway-signature should be a
   *standalone* lean role (use it INSTEAD of base, not on top) and the docs should say so, or
   base shouldn't swallow all of pathway-signature's skills. Right now the taxonomy implies
   composition that doesn't compose.

2. **[BLOCKER-for-discoverability] No CLI path from "my science" → "the right skill" without
   reading source.** `list skills` prints **names only** — no descriptions, no tags. The
   annotation skills I most needed (`mllmcelltype-consensus-annotation`,
   `cellxgene-census-annotation`) are guessable from names, but `coresh-signature-search`,
   `scired`, `gatom` are not. A busy scientist cannot answer "which skill scores a gene
   signature?" or "which does cell annotation?" from the CLI. `tags.yaml` has exactly these
   descriptions but is invisible to the tool.

3. **[BUG] `sciagent list tags` does not exist**, yet `inject --tag <name>` is a first-class
   command. You can inject by tag but cannot list the tags you're allowed to inject by — you
   must open `tags.yaml`. The tag system is half-surfaced.

4. **[CONFUSING] No scRNA-analysis role exists** that matches the single most common comp-bio
   workflow (QC → integration → annotation → DE/scoring). `base` is "everything" (36 skills,
   half of them ATAC/footprinting irrelevant to my project); `pathway-signature` is interpretation-
   only and subset-of-base; `scatac-regulatory` is chromatin. The taxonomy carves
   *base / chromatin / pathway-interp / architecture / library* — but the bread-and-butter
   **scRNA integration+annotation** lane has no dedicated role. I land on `base` and inherit a
   lot of ATAC noise (chromvar, pycistopic, crescendo, tf-footprint, scembed) I'll never use.

5. **[ORPHAN] The REAL project's live manifest is drifted from the catalog.**
   `/scratch/.../DC_hum_verse/.sciagent/manifest.json` has `"stack": ["architect","planning"]`
   and an injected skill `single-cell-vector-search`. **`planning` is not in `list roles`**
   (only architect/base/pathway-signature/scatac-regulatory/software-tool), and at the time the
   manifest's skill set predates the current catalog. A returning scientist re-running `status`
   in that real project would see a role that no longer exists. There's no migration/`doctor`
   path for stale manifests.

6. **[CONFUSING] `inject --tag annotation` pulls cross-modality skills into a scRNA project.**
   The `annotation` tag spans `scembed-atac-annotation` and ATAC reference-mapping. For my
   transcriptomics project that's over-injection. Tags need a modality facet (rna/atac) or the
   inject should warn which it's pulling.

7. **[NIT] `--title DC_hum_verse` doesn't reach the scientific-context.md H1.** The H1 reads
   "Scientific Context — DC_hum_verse" (good) but the title flag's effect is otherwise invisible;
   the body is a blank template. Minor, but the flag feels under-wired.

8. **[NIT] `sciagent new project --help` errors** (`unknown flag '--help'`) — usage only prints
   from `sciagent new --help`. Inconsistent help routing for a sub-sub-command.

9. **[POSITIVE] Status/validate/deactivate are genuinely excellent.** Provenance + shadowing +
   `via requires:` transitive resolution in `status`; idempotent tag inject; clean symlink +
   AGENTS.md-block teardown on deactivate. The *mechanics* are trustworthy; the *navigation* is
   what's thin.

10. **[POSITIVE] docs/_internal routing is clear and would orient an agent.** The README's
    path-table (scientific-context / reasoning / sessions / research), naming conventions, and the
    one-way reference rule (`_internal` may cite public, never reverse) are unambiguous. A scientist
    or an AI agent would know where a handoff vs a decision vs a lit-note goes.

## (d) Intuitiveness scorecard

| Dimension | Score | One-line justification |
|---|---:|---|
| **Role selection** (pick right role w/o reading YAML) | **6/10** | `list roles` blurbs got me to `base` correctly, but the base⊃pathway-signature subset trap and the missing scRNA-annotation lane mean I picked an over-broad role and an overlay that did nothing — and only *source-reading* revealed why. |
| **Discoverability** (find the skills I'm missing) | **3/10** | `list skills` is names-only; `list tags` doesn't exist; the rich `tags.yaml` descriptions are CLI-invisible. I found my annotation skills by eyeballing names + diffing, not by the tool guiding me. |
| **Docs-routing clarity** | **9/10** | `docs/_internal` table + naming conventions + one-way rule are crisp; clearly orients both human and agent. |
| **Zero-to-ready convenience** | **8/10** | `new project` → `activate base` → working harness in two commands, clean symlinks, AGENTS.md managed block, honest "Next steps." Teardown is clean. Only friction is over-injected ATAC skills and the redundant overlay. |

---

### Bottom line
The **machinery** (scaffold, activate, status provenance, requires-resolution, tag-inject,
deactivate) is solid and trustworthy. The **map from science to harness** is where a busy
comp-bio stalls: there's no dedicated scRNA-integration-and-annotation role, `pathway-signature`
silently duplicates `base`, and the only way to discover "which skill does X" is to leave the CLI
and read `tags.yaml`. Surfacing tags/descriptions in `list` and adding a real scRNA-annotation
role (or making pathway-signature a lean standalone) would move role-selection and discoverability
from "read the source" to "ask the tool."
