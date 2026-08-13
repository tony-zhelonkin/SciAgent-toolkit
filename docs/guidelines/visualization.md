# Visualization Standards — REDIRECTED

> **This document is superseded.** The `base_size = 12` / `font.size: 12` /
> `theme_publication(base_size=12)` guidance defined here is no longer authoritative and
> must not be used. The `save_publication_plot` helper defined here is retired.

---

## Single figure SSOT

All figure styling, saving, and captioning is governed by the **figure-style contract**.
Consult these three sources — they are the only authoritative references:

| What | Where |
|------|-------|
| Contract skill (usage, anti-patterns, done-when) | `skills/figure-style/SKILL.md` |
| Contract implementation (R + Python helpers) | `lib/figure-style/figure_helpers.{R,py}` |
| Per-project geometry + font floors | `02_analysis/config/analysis_config.yaml` under `figures:` (`base_size` floor 14 pt — one legible tier) |

The per-project shim (`02_analysis/helpers/figure_style.{R,py}`) is the import entry point
for every viz script. Bind the toolkit catalog with `sciagent link`.

This file is kept as a stub so that old links resolve; its previous content has been removed.
