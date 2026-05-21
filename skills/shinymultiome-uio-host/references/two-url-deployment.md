# Two-URL deployment — CellxGene + ShinyMultiome alongside

## The pattern

```
Wet lab opens two browser tabs:
   ┌────────────────────────────────┐    ┌────────────────────────────────┐
   │  http://<HOST>:8080/           │    │  http://<HOST>:8090/           │
   │  CellxGene (scrna-cxg-host)    │    │  ShinyMultiome.UiO (this skill)│
   │  - RNA UMAP + autosave annot   │    │  - ATAC tracks                 │
   │  - violin/dotplot              │    │  - Peak2Gene links             │
   │  - feature search              │    │  - Footprints (optional)       │
   └────────────────────────────────┘    └────────────────────────────────┘
                  │                                    │
                  │       same htpasswd file           │
                  ↓                                    ↓
              writes to                           reads from
              annotation/<dataset>/                  annotation/htpasswd
```

Same host, same auth, two URLs. The wet lab labels in CellxGene (autosave persists) and views chromatin in ShinyMultiome (read-only, ephemeral).

## Why two URLs and not one

Three forces push us this way:

1. **CellxGene has no ATAC concept.** Schema does not represent ChromatinAssay; the H1 .h5ad is RNA-only.
2. **ShinyMultiome.UiO is read-only by design** — no autosave, no annotation tooling. Pretending to label there is a footgun.
3. **Trying to combine forces VIP** (cellxgene-VIP), which is a fragile patch on top of cellxgene's frontend (version-pinned to 1.1.1, breaks on upstream bumps). The two-URL path is the lower-risk default.

The trade-off is wet-lab UX: two URLs to remember, two browser tabs to context-switch between. Mitigation: ship a one-page wet-lab onboarding doc with both URLs side-by-side and one paragraph each describing what they're for.

## Why share htpasswd

The cohort that needs CellxGene access is the same cohort that needs ShinyMultiome access. Two separate htpasswd files would mean rotating credentials in two places — invariably one drifts and a wet-lab user is locked out of one tab without knowing why.

The skill's Decision Pause S4 default is **share** — symlink `shinymultiome/annotations/htpasswd` to `scrna-cxg-host/annotations/htpasswd`. To compartmentalise (e.g., chromatin viewer only for ATAC team), use Option B.

## Resource sharing on the same host

Both stacks run on the same Docker host. The cellxgene container has its own memory limit (typically 30 GB), the shinymultiome container has its own (typically 32 GB). Together: 62 GB of `limits.memory` allocated. Most servers have headroom; on a 64 GB host these two stacks alone leave nothing for the OS — bump to 128 GB host for comfort.

The cellxgene cohort is "many small queries" (each request is a few MB of data per panel). ShinyMultiome is "fewer, heavier queries" (each Coverage redraw is 1–3 GB transient). They contend for I/O on the disk holding the data dir but not on Docker network.

## Deployment ordering

The H1 (scrna-cxg-host) deploy goes first — it's the annotation tool. ShinyMultiome is the **chromatin companion**, not the standalone deliverable. If the wet lab cohort cannot label, the chromatin viewer is academic.

```
Phase 1 (week 1): scrna-cxg-host up, wet lab labelling RNA UMAP
Phase 2 (week 2): ShinyMultiome up, wet lab now also viewing chromatin tracks
                  for the labels they just defined
```

Both stacks read from the same `03_results/` directory. Annotations made in CellxGene end up in `03_results/annotation/<dataset>/`; ShinyMultiome ignores this directory entirely (it just shows chromatin for whatever labels exist on the .rds at deploy time).

## Refreshing ShinyMultiome with new labels

ShinyMultiome is **read-only**. New CellxGene annotations do not flow into it automatically. To push new labels:

1. Pull new labels from CellxGene's annotation CSV: `03_results/annotation/<dataset>/<user>.csv`
2. Re-run Phase A.4 (group.by validation) on the .rds, ensuring the new label column is present in `obj@meta.data`
3. Re-save the .rds
4. `docker compose -f docker-compose.shinymultiome.yml down`
5. `docker compose -f docker-compose.shinymultiome.yml up -d`

The container reload cost is the .rds load time at startup (30–120s depending on object size). Schedule these refreshes during low-usage windows (e.g., overnight).

A future enhancement (see `docs/multiome-deploy/LONG-TERM-VISION.md`): a "label sync" daemon that watches the CellxGene annotation dir and triggers ShinyMultiome refreshes automatically. Out of scope for this skill.

## What if the wet lab wants annotation in the chromatin viewer

This is Decision Pause S4-equivalent in the H2 architecture, escalated to LONG-TERM-VISION.md territory:

- **Option A (recommended):** keep the read-only model. Wet lab labels in CellxGene, views in ShinyMultiome. Two URLs.
- **Option B (heavier):** patch ShinyMultiome to write annotations back to a dir; build a sync layer that pushes them into CellxGene's annotation format. Out of scope for this skill — it's a 1–2 month project.
- **Option C (alternate):** drop ShinyMultiome and use cellxgene-VIP, which embeds chromatin tracks in the same UI as the labelling tool. Trades VIP's maintenance fragility for UX consolidation. See `docs/multiome-deploy/PHASE3-VIP-DEPLOY.md`.

The skill defaults to A. Document the trade-off explicitly in onboarding so wet lab does not expect autosave in the chromatin viewer.
