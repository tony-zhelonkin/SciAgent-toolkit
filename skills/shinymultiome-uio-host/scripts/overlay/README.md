# scripts/overlay/ — Signac-compatible global.R overlay

## What this directory is

A minimal patch on top of upstream `EskelandLab/ShinyMultiomeUiO`. The Dockerfile clones upstream verbatim, then **replaces** `global.R` with `overlay/global.R`. Upstream's `server.R`, `ui.R`, and `app.R` are untouched.

## What the overlay changes vs upstream `global.R`

| Upstream behaviour | Overlay behaviour |
|---|---|
| `library(BSgenome.Hsapiens.UCSC.hg38)` (hard-coded) | `library(${BSGENOME_PKG})` (env var) |
| `seuratObject = "/Users/akshay/.../*.RDS"` | `RDS_PATH = Sys.getenv("RDS_PATH")` |
| `fragFilePath = "/Users/akshay/.../*.tsv.gz"` | `FRAGMENTS_DIR = Sys.getenv("FRAGMENTS_DIR")` (rewrites every Fragment object's path by basename match inside that dir) |
| One hard-coded fragment file | All Fragment objects on the loaded `pbmc` get `UpdatePath` against `FRAGMENTS_DIR` (multi-sample multiome works) |
| No assay-name validation | Fails fast at startup if `ATAC_ASSAY` is absent; warns if `RNA_ASSAY` is absent |
| No annotation seqlevel check | Warns at startup if seqlevels look Ensembl-style |

## Why we do not fork the whole repo

Upstream is small (~6 R files) and rarely changes. Forking forces us to track upstream commits. Overlaying only `global.R` keeps the surface to bump small: bump `SHINYMULTIOME_PIN` in `Dockerfile.shinymultiome`, rebuild, done.

## Pinning

The Dockerfile clones upstream at `--branch ${SHINYMULTIOME_PIN}` (default `main`). For reproducibility on a deployed server, pin to a commit SHA:

```bash
docker compose build \
  --build-arg SHINYMULTIOME_PIN=<commit-sha-or-tag> \
  shinymultiome
```

Currently verified pin: `main` as of 2026-05 (last upstream push 2023-07-21 — repo is stable).

## When the overlay needs updating

If upstream introduces breaking changes to `server.R` or `ui.R` that depend on values previously read from the hard-coded `global.R`, the overlay must mirror those changes. Trigger conditions:

1. Upstream renames the loaded variable from `pbmc` to something else.
2. Upstream introduces new `Sys.getenv` calls of its own (collision risk).
3. Upstream's `server.R` switches to multi-sample fragment lists in a new shape.

Watch upstream commits; rebuild + smoke-test if any of the above land.

## Smoke test

After bumping `SHINYMULTIOME_PIN`:

```bash
docker compose build shinymultiome
docker compose run --rm shinymultiome \
  R -e "source('/srv/shiny-server/app/ui.R'); cat('UI sourced OK\n')" \
  -e "source('/srv/shiny-server/app/server.R'); cat('server sourced OK\n')"
```

Both lines must print before the container exits. A `server.R` source-error in this run reveals upstream API drift.
