# Shiny memory profile and timeouts

## What lives in container RAM

ShinyMultiome.UiO's `global.R` does `pbmc <- readRDS(RDS_PATH)` at app startup, **once per container**. The Seurat object stays in RAM for the lifetime of the container; Shiny sessions share the same global. Per-session reactivity adds only small overhead (per-user input state).

The fragment files are **not** loaded into RAM — `Signac::CoveragePlot` reads them via tabix random-access. This is what keeps the memory ceiling tractable on multi-sample multiomes.

## What's heavy

For a typical 10x Multiome object:

| Slot | Size |
|------|------|
| RNA assay counts (sparse) | ~1–3 GB per 50k cells × 30k genes |
| ATAC peaks counts (sparse) | ~3–8 GB per 50k cells × 200k peaks |
| Reductions (UMAP, LSI, PCA) | <100 MB |
| Fragments (object metadata only — paths) | <1 MB |
| `Annotation` (GRanges, ~1M ranges for mouse/human) | ~30–80 MB |
| `Links` (GRanges, ~10k–100k links) | <100 MB |
| `@motifs` (Motif object) | ~5–20 MB per 600 PFMs |
| `@positionEnrichment` (footprint SE) | ~50–200 MB per motif × group.by combo |

A 50k-cell mouse multiome with no Links, no Motifs typically uses **6–12 GB RAM** at app idle. Add Links → +0.5 GB. Add 50 motif footprints across 20 cell types → +2 GB. Coverage track render adds 1–3 GB transient (released between requests).

## Sizing guidance

| Cell count | Memory limit | Memory reservation | Notes |
|------------|--------------|---------------------|-------|
| <30k       | 8G           | 2G                  | Plenty of headroom |
| 30k–80k    | 16G          | 4G                  | Default for `.env.template` is `32G` to be safe |
| 80k–200k   | 32G          | 8G                  | Skill's default |
| 200k–500k  | 64G          | 16G                 | Strip `@motifs` + `@positionEnrichment` if not used |
| >500k      | benchmark first | — | Consider downsampling for the Shiny instance, or fall back to per-celltype subsets |

Always verify with `docker stats --no-stream` during typical wet-lab usage. The "first session loads in" pattern is misleading — long-running track renders peak memory transiently and the limit must accommodate the peak.

## Memory pitfalls

### Cold-start OOM

A container with the limit set just at idle-RAM can OOM on the first `CoveragePlot` because the transient memory exceeds the limit. **Always set the limit ~30% above measured idle**.

### Multiple Shiny sessions

Shiny in single-process mode shares one R process across all concurrent sessions. The `pbmc` global is shared. Per-session reactives are small, but per-session `CoveragePlot` renders happen **in the same R process** — they queue, not parallel. Two simultaneous Coverage redraws across cell-type-rich `group.by` columns can serialize for ~2 minutes.

For heavy concurrent usage, consider running multiple shinymultiome containers behind a load-balancer config. Out of scope for this skill.

### Garbage that doesn't get collected

Signac's plot internals create large `data.table`s during `CoveragePlot`. R's GC is conservative; with a 32 GB limit, you may see RAM creep upward over a long-lived session before the GC reclaims. Symptoms: container memory at 28 GB after 10 redraws, even when each plot renders in a few seconds. Fix: `gcinfo(TRUE)` in dev, or schedule `gc()` calls on a Shiny `observe` timer (out of scope of the skill but documented as a knob).

## nginx timeouts

Default nginx proxy timeouts (60s) are too short for ShinyMultiome's worst-case redraws:

| Scenario | Typical redraw time | Notes |
|----------|---------------------|-------|
| Coverage on 1 gene × 5 cell types × 10kb | 3–8 s | Fits within default |
| Coverage on 1 gene × 30 cell types × 10kb | 30–90 s | **Exceeds 60s** |
| Coverage on 1 gene × 30 cell types × 200kb | 90–180 s | Exceeds default |
| Footprint on 1 TF × 10 cell types | 60–180 s | Exceeds default |
| Feature plot (UMAP) on 1 gene | 1–3 s | Fits |

The shipped `nginx-shinymultiome.conf` sets:

```
proxy_read_timeout    300s;
proxy_send_timeout    300s;
proxy_connect_timeout 60s;
```

5-minute read/send timeout covers all the above scenarios. The connect timeout stays at 60s because that's a fast-fail for unreachable backends (no need to bump).

## WebSocket disconnect

Shiny's reactivity uses WebSocket frames to push UI updates. nginx must:

1. Set `proxy_http_version 1.1;`
2. Pass `Upgrade` and `Connection: upgrade` headers.
3. Disable buffering (`proxy_buffering off`) so frames arrive promptly.

Without (1)+(2): the browser reports "Disconnected from server" within ~30s.
Without (3): the UI feels laggy — frames arrive in chunks of 4–8 KB.

The shipped config sets all three. Verify with the WebSocket smoke test in `validate_shinymultiome_env.sh`.

## Long-running session keepalive

If a single Shiny session sits idle for >2 minutes (default `Shiny:::application$timeoutInterval`), the connection drops. Wet lab users coming back from coffee see "Disconnected". Two mitigations:

1. Bump Shiny's session timeout in the overlay (`global.R`):
   ```r
   options(shiny.idleTimeout = 60 * 30)   # 30 min
   ```
2. The browser's heartbeat ping keeps the WebSocket alive automatically when the tab is foregrounded; the issue is only on backgrounded tabs.

Not currently set in the skill's overlay (avoids divergence from upstream behaviour). Add it if the wet lab reports drops.
