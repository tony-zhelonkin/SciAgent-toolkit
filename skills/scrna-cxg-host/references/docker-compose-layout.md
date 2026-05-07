# Docker Compose layout

The reference deployment runs N pairs of (cellxgene container + nginx container). Each pair is independent — separate ports, separate annotations directory, separate dataset file. They share the cellxgene image, the htpasswd file (or per-dataset variants), and the host data path (read-only).

## Topology

```
                     ┌──────────────┐                ┌──────────────┐
host browser ──HTTP──┤ nginx-full   ├──reverse-proxy─┤ cellxgene-full│ -- /data:ro
   :8080             │ :8080        │  + basic auth  │ :5005          │ -- /annotations:rw
                     └──────────────┘                └──────────────┘
                     ┌──────────────┐                ┌──────────────┐
host browser ──HTTP──┤ nginx-t-cells├──reverse-proxy─┤ cellxgene-t-cells│ -- /data:ro
   :8081             │ :8081        │  + basic auth  │ :5006          │ -- /annotations:rw
                     └──────────────┘                └──────────────┘
                     ... one pair per dataset ...
```

cellxgene listens *only* on the docker-network interface (`5005+`); nginx listens on the host (`8080+`). The browser never sees `5005`. This is what `--host 0.0.0.0 --port <internal>` accomplishes — within-network, not host-network.

## Per-cellxgene container

```yaml
cellxgene-<dataset>:
  image: local/cellxgene:latest          # built from Dockerfile (Python 3.10-slim + pip install cellxgene)
  command: >
    launch
    --backed
    --title "<<PROJECT_NAME>> • <dataset display name>"
    --annotations-dir /annotations
    --max-category-items 5000
    --host 0.0.0.0
    --port <<CXG_PORT_<dataset>>>
    /data/<dataset>.h5ad
  user: "<<UID>>:<<GID>>"                 # match host owner; prevents root-owned annotation files
  volumes:
    - <<DATA_HOST_PATH>>:/data:ro
    - <<ANNOTATIONS_HOST_PATH>>/<dataset>:/annotations
  deploy:
    resources:
      limits:
        memory: 30G                       # tune for project; larger datasets need more
      reservations:
        memory: 8G
  restart: unless-stopped
```

Key flags:

- **`--backed`** — keeps `X` on disk; required for AnnData larger than memory
- **`--annotations-dir /annotations`** — autosave path; must be a host-mounted volume, not container ephemeral
- **`--max-category-items 5000`** — UI handles up to 5000 unique values per category; below this, the discrete-color picker works; above, CXG falls back to a free-text label that confuses wet-lab users
- **`user: "${UID}:${GID}"`** — host-aligned UID; without this, autosave files end up `root:root` and the wet lab cannot delete them from the host

## Per-nginx container

```yaml
nginx-<dataset>:
  image: nginx:alpine
  depends_on:
    - cellxgene-<dataset>
  ports:
    - "<<NGINX_PORT_<dataset>>>:<<NGINX_PORT_<dataset>>>"
  volumes:
    - ./nginx/nginx-<dataset>.conf:/etc/nginx/nginx.conf:ro
    - ./nginx/htpasswd:/etc/nginx/htpasswd:ro
  restart: unless-stopped
```

`depends_on` does not wait for the cellxgene server to *be ready* — it just orders container starts. Nginx will return 502s until cellxgene's server loop is up; this typically takes 5–30 seconds for `--backed` mode on large files.

## How port pairs scale

`docker-compose.yml.template` uses `<<CXG_PORT_<dataset>>>` and `<<NGINX_PORT_<dataset>>>` placeholders, one pair per dataset. The `render_compose.py` helper takes `--datasets` mappings and emits one pair per. Defaults: `CXG_PORT_FULL=5005`, `NGINX_PORT_FULL=8080`, then `+1` per additional dataset.

The Decision Pause for port collisions intercepts when defaults are taken. After remap, the resolved ports go to `.env` and `analysis_config.yaml::decisions::scrna-cxg-host::ports`.

## Resource sizing heuristics

| Dataset size | `limits.memory` | `reservations.memory` |
|--------------|-----------------|----------------------|
| ≤ 50k cells  | 8G              | 2G                   |
| 50k–200k     | 30G (default)   | 8G                   |
| 200k–500k    | 64G             | 16G                  |
| ≥ 500k       | benchmark first; consider per-celltype subsets to keep each instance under 64G |

`--backed` mode dramatically reduces resident memory but the inferred-dtype + sparse-encoding overhead is still ~2–3× nominal. Watch `docker stats` during a typical wet-lab session.

## Image build

The skill ships [`scripts/Dockerfile.cellxgene`](../scripts/Dockerfile.cellxgene) — copy it verbatim into your deploy directory. Key elements:

```dockerfile
FROM python:3.10-slim
ENV PIP_NO_CACHE_DIR=1 PYTHONUNBUFFERED=1
RUN apt-get update && apt-get install -y --no-install-recommends build-essential curl \
 && rm -rf /var/lib/apt/lists/*

# Pin numpy + pandas alongside cellxgene so wheels resolve to a single ABI
# generation. cellxgene 1.2.0's own pins (numpy>1.22, pandas<2.0.0) under-
# constrain the resolver and let it pair numpy 2.x with pandas 1.5.x — those
# wheels disagree on PyArray_Descr (88 vs 96 bytes) and pandas crashes at
# import. See SKILL.md "numpy/pandas binary incompatibility" pitfall.
RUN pip install --no-cache-dir \
      "numpy==1.23.5" \
      "pandas==1.5.3" \
      "cellxgene==1.2.0" \
 && mkdir -p /data /annotations

# Build-time smoke test — fail the build, not runtime, if the ABI is wrong.
RUN python -c "import pandas, numpy; from server.cli.cli import cli; \
print('cellxgene env OK; numpy', numpy.__version__, 'pandas', pandas.__version__)"

EXPOSE 5005
ENTRYPOINT ["cellxgene"]
CMD ["launch", "--host", "0.0.0.0", "--port", "5005"]
```

`build-essential` covers the wheel-less builds among `cellxgene`'s deps; `curl` is convenient for in-container debugging. The image stays under 500MB.

`compose.yml` references the image as `image: local/cellxgene:latest` (or a project-specific tag) and is built once via `docker compose build cellxgene-full`. Bumping `cellxgene` past 1.2.0 (e.g. to 1.3.0) lets you drop the explicit numpy/pandas pins — 1.3.0's upstream `requirements.txt` already pins `numpy==2.0.1` + `pandas>=2.2.2`.

## Tear down

```bash
docker compose down                  # stop and remove containers, keep volumes
docker compose down -v               # also remove named volumes (annotations are bind-mounted, so they survive)
```

Bind-mounted host paths (`./annotations`, the data dir) are *not* deleted by `down -v`; they live on the host fs and are durable across deploys.
