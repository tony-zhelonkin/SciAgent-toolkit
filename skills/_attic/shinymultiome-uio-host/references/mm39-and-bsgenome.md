# mm39 / GRCm39 and BSgenome — what changes vs the upstream hg38 default

## Upstream's hard-coded assumption

`EskelandLab/ShinyMultiomeUiO/global.R` line 18 (verbatim):

```r
library(BSgenome.Hsapiens.UCSC.hg38)
```

This is the only line that assumes a specific genome. The rest of the app uses the loaded BSgenome object indirectly via Signac:

- `Signac::CoveragePlot` reads sequence from the loaded BSgenome (matched by `genome(annotation)`).
- `Signac::AddMotifs(obj, genome = bsg, ...)` (footprint tab) reads sequence to scan motifs.
- `Signac::Footprint(obj, genome = bsg, ...)` reads sequence around motif hits.

For mm39 / GRCm39, three things change:

1. The `library()` call in `global.R` (the skill overlays this).
2. The Bioconductor package install in the Dockerfile (parameterised via `BSGENOME_PKG`).
3. The matching `EnsDb` package for `Annotation` rebuild, if Phase A.2 needs to construct one (`EnsDb.Mmusculus.v107` or later).

## Available mouse BSgenome packages

```
BSgenome.Mmusculus.UCSC.mm10   — GRCm38, ArchR-supported, older analyses
BSgenome.Mmusculus.UCSC.mm39   — GRCm39, our project's assembly
BSgenome.Mmusculus.UCSC.mm9    — legacy
```

`BSgenome.Mmusculus.UCSC.mm39` was added to Bioconductor in release 3.16 (Oct 2022). The skill's Dockerfile uses `bioconductor/bioconductor_docker:RELEASE_3_19`, which is comfortably above that floor.

## Available mouse EnsDb packages (mm39 / GRCm39)

```
EnsDb.Mmusculus.v79     — Ensembl 79, GRCm38 — DO NOT USE for mm39
EnsDb.Mmusculus.v107    — Ensembl 107, GRCm39 — first GRCm39 release
EnsDb.Mmusculus.v108+   — newer Ensembl on GRCm39
```

Pick the EnsDb release **that was used during the analysis** — gene models change between Ensembl releases (new transcripts, retired genes), and a mismatch makes track lookups by gene symbol return inconsistent regions vs the analysis. Default in `.env.template` is `EnsDb.Mmusculus.v107` (the first GRCm39 release); bump if your analysis used a later one.

## Why ArchR is not a path here

ArchR ships pre-built genome annotations only for `hg19, hg38, mm9, mm10`. mm39 requires `createGenomeAnnotation(genome = "mm39")` + `createGeneAnnotation(TxDb = ..., OrgDb = ...)` — a 30-line custom build that adds ~600 MB to the container image.

The skill's deployment is Signac-based, not ArchR-based. Signac's `CoveragePlot` works on any genome that has a BSgenome package, with no ArchR custom-build dance. This is why we picked Signac over ArchR for the H2 deploy (see `docs/multiome-deploy/PHASE3-ALT-SHINYMULTIOME.md` §2 for the full reasoning).

## Building a custom BSgenome (only if needed)

For non-standard assemblies (newer than what Bioconductor ships, or in-house genomes), build a custom BSgenome:

```r
# Outside the container, on a build machine:
library(BSgenome)
seed <- "<custom_seed_file.dcf>"        # describes provider, assembly, FASTA path, etc.
forgeBSgenomeDataPkg(seed)
# Produces BSgenome.<Org>.<Provider>.<Assembly>_<version>.tar.gz
```

Then in the skill's Dockerfile, install from the local tarball:

```dockerfile
COPY bsgenome-custom.tar.gz /tmp/
RUN R -e "install.packages('/tmp/bsgenome-custom.tar.gz', repos = NULL)"
```

The `.env::BSGENOME_PKG` should match the package name produced by `forgeBSgenomeDataPkg`. This adds a dependency on a build artifact that is not in any registry — document the seed file and the rebuild procedure in the deploy repo.

## Image build time on first mm39 build

Approximate timings on a typical CI node (4 cores, 16 GB RAM):

| Stage | Time | Image size delta |
|-------|------|-------------------|
| `bioconductor/bioconductor_docker:RELEASE_3_19` pull | 2–4 min | 1.5 GB |
| CRAN packages (Seurat, shiny, ggplot2, ...) | 8–15 min | +400 MB |
| Bioc packages (Signac, GenomicFeatures, ...) | 6–10 min | +500 MB |
| `BSgenome.Mmusculus.UCSC.mm39` install | 2–4 min | +650 MB |
| `EnsDb.Mmusculus.v107` install | 1–2 min | +120 MB |
| Upstream ShinyMultiomeUiO clone + overlay | <30 s | +15 MB |
| **Total** | **30–45 min** | **~3.0–3.2 GB** |

Subsequent rebuilds (without cache invalidation) finish in <2 min. To preserve the cache across BSgenome/EnsDb bumps, see `Dockerfile.shinymultiome`'s explicit ARG layer split.
