# TE Reference Build Record — `<GENOME_BUILD>`

> **Recipe:** `te-reference-saf-build` v`<X.Y.Z>`; build script:
> `scripts/build_te_saf.sh`. Build-once, shared across every dataset on this genome build.
>
> One genome build = one record. Fill a copy each time the TE reference is (re)built and keep
> it beside the artifacts (e.g. `/data1/shared/ref/<species>/<provider>/<build>/`). This closes
> the acquisition + build provenance gaps flagged in
> `docs/_internal/research/08-shared-artifact-provenance-and-skill-decision.md` §B.

---

## 1. Identity

| Field | Value |
|---|---|
| Genome build | `<e.g. mm39 / GRCm39>` |
| Species / provider | `<e.g. Mus musculus, Ensembl>` |
| Gene annotation | `<e.g. GENCODE vM37>` |
| Build date | `<YYYY-MM-DD>` |
| Operator | `<name>` |
| Home directory | `<e.g. /data1/shared/ref/mouse/Ensembl/mm39/>` |

## 2. TE GTF source (DOWNLOADED, not built)

- Tool / provider: TEtranscripts (https://www.mghlab.org/software/tetranscripts), pre-generated.
- File: `<e.g. GRCm39_Ensembl_rmsk_TE.gtf.gz>`
- Download URL: `<exact Dropbox/source URL>`
- Download date: `<YYYY-MM-DD>`
- md5: `<md5sum of the downloaded .gz>`

## 3. Gene GTF used for the exon BED

- File: `<e.g. gencode.vM37.primary_assembly.annotation.gtf.gz>`
- md5: `<md5sum>`
- Provenance: `<derived freshly from the canonical gene GTF — RECOMMENDED. Do NOT borrow a
  filtered GTF from a sibling project's outdir, see artifact 08 §A.3>`

## 4. Class composition choice

- `<keep ALL classes (canonical) | retro-only via --keep-classes '^(LINE|SINE|LTR|RC)$'>`
- Rationale: `<canonical = keep all, filter downstream>`

## 5. Exact build command (VERBATIM)

```bash
scripts/build_te_saf.sh \
  --te-gtf   <.../GRCm39_Ensembl_rmsk_TE.gtf.gz> \
  --gene-gtf <.../gencode.vM37.primary_assembly.annotation.gtf.gz> \
  --out-dir  <.../mm39> \
  --prefix   GRCm39_rmsk_TE
  # [--keep-classes '<regex>'] [--strip-chr yes|no]
```

## 6. Outputs (paths + md5s)

| Artifact | Path | md5 |
|---|---|---|
| Grouped SAF (intermediate) | `<.../GRCm39_rmsk_TE_GROUPED_all.saf>` | `<md5>` |
| Grouped no-exon SAF (shared input) | `<.../GRCm39_rmsk_TE_GROUPED_all_noExon.saf>` | `<md5>` |

## 7. Verification

| Check | Value |
|---|---|
| Unique groups (grouped SAF) | `<e.g. 1243>` |
| TE loci overlapping exons (pre-subtract) | `<n>` |
| Residual exon overlap (post-subtract) | `<MUST be 0>` |
| Unique groups (no-exon SAF) | `<≈ same as grouped>` |
| Software versions | gawk `<v>`, bedtools `<v>` |

## 8. Datasets consuming this build

> List the datasets that point at this no-exon SAF (so a rebuild's blast radius is known).

- `<e.g. 13036-DM, AdaW_eWAT_WL, 14839-DM>`
