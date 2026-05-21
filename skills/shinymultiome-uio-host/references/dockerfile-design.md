# Dockerfile design — layering and ARG / ENV strategy

## Base image choice: `bioconductor/bioconductor_docker:RELEASE_3_19`

Three alternatives evaluated, then rejected:

| Base | Why rejected |
|------|--------------|
| `python:3.10-slim` + manual R install | R + Bioc compile from source = 3–4× build time |
| `rocker/r-ver:4.4` | No Bioc deps cached; same compile cost as above |
| `rocker/shiny-verse:latest` | Has Shiny + tidyverse but lacks Bioc — Signac install pulls 30+ Bioc deps from source |

`bioconductor/bioconductor_docker:RELEASE_3_19` ships with R 4.4, BiocManager, and most Bioc dep system libraries pre-installed. Image starts at ~1.5 GB but Signac install is a few minutes (vs ~30 min on rocker).

`RELEASE_3_19` is the Bioconductor 3.19 release (April 2024). Bumps:

- Bioc 3.20 (Oct 2024) — Signac 1.14, Seurat 5.1, mm39 EnsDb v110.
- Bioc 3.21 (May 2025) — newer EnsDb releases.

Bumping `RELEASE_3_xx` is a one-line change in the Dockerfile + a smoke-test cycle. Pinning to `RELEASE_3_19` keeps the deploy reproducible across builds; bump deliberately, not accidentally.

## Layer ordering

The Dockerfile is structured so that the most-likely-to-change values invalidate the smallest cache:

```
FROM bioconductor/bioconductor_docker:RELEASE_3_19      # base layer
RUN apt-get install ...                                  # rarely changes
RUN install.packages('Seurat', 'Signac', ...)            # changes on Seurat/Signac bump
RUN BiocManager::install(BSGENOME_PKG)                   # changes on genome change ← MOST LIKELY
COPY overlay/global.R /srv/shiny-server/app/global.R     # changes on overlay edit
```

Genome change → only the BSgenome layer rebuilds. Overlay edit → only the COPY layer rebuilds.

`BSGENOME_PKG` is passed via `ARG` and re-exported as `ENV`. The `ARG` controls the build; the `ENV` ensures the package name is also visible at runtime so the overlay's `library(BSGENOME_PKG)` can pick it up dynamically (it does `library(get(BSGENOME_PKG))` style — actually `library(BSGENOME_PKG, character.only = TRUE)`).

## Two-stage builds — considered, rejected

A two-stage build (Stage 1: install R + Bioc, Stage 2: copy site-library to a runtime image) was evaluated and rejected. Reasons:

1. The runtime image still needs R + system libs; you cannot drop to a `slim` runtime base.
2. The site-library copy across stages is fragile — Bioc packages frequently include compiled artifacts whose paths embed the original `R_HOME`.
3. The savings are small (~200 MB) for the complexity added.

The current single-stage build is simpler and cache-friendly.

## Image tagging

The compose template tags the image as `local/shinymultiome-uio:latest`. For deploy reproducibility, retag after build:

```bash
docker compose -f docker-compose.shinymultiome.yml build
docker tag local/shinymultiome-uio:latest \
           local/shinymultiome-uio:mm39-bioc3.19-pin-$(git rev-parse --short HEAD)
```

The skill does not automate retagging — handled by the deploy operator.

## Build-time smoke test

The Dockerfile ends with:

```dockerfile
RUN R --no-save -e " \
        suppressPackageStartupMessages({ \
          library(shiny); library(Seurat); library(Signac); \
          library(get(Sys.getenv('BSGENOME_PKG')), character.only = TRUE); \
        }); \
        cat('shinymultiome image OK; ...\n')"
```

This runs *during* `docker build` — fails the build, not the runtime, on missing genome packages or version drift. Without this, a misnamed `BSGENOME_PKG` produces a runtime crash on first user session, which is much harder to debug.

## What the Dockerfile does NOT do

- **No SSL termination** — out of scope; nginx-shinymultiome runs HTTP only. Wrap with Caddy/Traefik upstream if SSL needed.
- **No PAM auth** — basic-auth via htpasswd only. Out of scope.
- **No hot-reload** — overlay changes require `docker compose down && up`. Bind-mounting the overlay for dev (instead of `COPY`) is documented in `overlay/README.md` but not the default.
- **No multi-arch** — built for `linux/amd64` only. ARM hosts (M1/M2 Macs) will work via emulation but build time triples.
- **No tini / signal forwarding** — Shiny handles SIGTERM correctly; the `tini` init wrapper is not required for the simple `R -e shiny::runApp()` command.

If any of these become hard requirements, document the change in the project's deploy notes — do not modify the skill's Dockerfile silently.

## Rebuild triggers

| Trigger | Layer rebuild | Time |
|---------|---------------|------|
| Bump `bioconductor_docker:RELEASE_3_xx` | All layers above | 30–45 min |
| Bump `Seurat` / `Signac` versions | CRAN/Bioc layers + below | 10–15 min |
| Change `BSGENOME_PKG` / `ENSDB_PKG` | Genome layer + below | 5–10 min |
| Bump `SHINYMULTIOME_PIN` (upstream commit) | Clone + overlay layers | <1 min |
| Edit `overlay/global.R` | Overlay COPY layer | <30 s |

Rebuild during deploy windows; never push a fresh image to a running deploy without redeploying.
