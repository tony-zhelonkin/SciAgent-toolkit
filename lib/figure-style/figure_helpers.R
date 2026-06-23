## figure_helpers.R — the SciAgent-toolkit cross-language FIGURE-STYLE CONTRACT (R side).
## =====================================================================================
## ONE place that owns the project's figure format so phase viz scripts never reinvent it. This
## is the R half of a two-language contract; the Python half (`figure_helpers.py`, same dir) has
## FUNCTION PARITY — identical public names, equivalent semantics — and both read the SAME
## `analysis_config.yaml:figures` block. Centralizing styling here is the load-bearing capability
## behind the owner's #1 recurring pain point (figure legibility) and #2 (results placement).
##
## Config keys read (from `analysis_config.yaml:figures`):
##   base_size, title_size, axis_title_size, axis_text_size, strip_size, legend_text_size,
##   label_size, line_width, point_size, base_size_column, width, height, width_column,
##   height_column, dpi, top_n, volcano_label_top, z_clamp, nes_cap, caption_wrap_column,
##   variants, by_contrast_dir, overview_dir
## Plus, from elsewhere in the config: `paths.results` (results root), `paths.master`
## (master-table root), `paths.stage_tables_subdir` / `paths.stage_figures_subdir`.
##
## LAZY HEAVY-DEP DESIGN (important — read before editing):
##   This file MUST be `source()`-able and its path / caption / table / config functions callable
##   WITHOUT ggplot2 / cairo on a bare box. Plotting deps (ggplot2, cairo) are loaded LAZILY,
##   inside the function that needs them via requireNamespace() — never at top level — so a parity
##   check / code review on a toolchain-less box can source this file. Plotting functions degrade
##   to a clear stop() when ggplot2 is absent. There are NO top-level side effects.
##
## Reuse (one source, a few calls):
##   source("01_modules/SciAgent-toolkit/lib/figure-style/figure_helpers.R")
##   cfg <- load_figure_config("02_analysis/config/analysis_config.yaml")
##   p <- ggplot(...) + ... + project_theme(config = cfg)        # the SINGLE theme entry point
##   save_overview(p, "04_gsea", "gsea_hallmark_heatmap", table = df,
##                 finding = "Hallmark IFN-alpha/gamma dominate the ISD90 response.",
##                 script = "02_analysis/scripts/11_gsea_viz.R", fn = "save_overview",
##                 config_kv = "figures.nes_cap = 3.5", input = "03_results/objects/gsea.rds",
##                 how_to_read = "Rows = pathways; color = NES (orange up / blue down).",
##                 config = cfg)            # figure + sibling table + caption, atomic

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

## Per-variant fallback FLOORS, used only when a key is absent from the project config; the
## project config is authoritative. Keep in sync with figure_helpers.py:_FIG_DEFAULTS + template.
.FIG_DEFAULTS <- list(
  base_size = 16, title_size = 18, axis_title_size = 15, axis_text_size = 13, strip_size = 14,
  legend_text_size = 13, label_size = 5, line_width = 0.8, point_size = 2.0,
  base_size_column = 9, width = 10, height = 8, width_column = 3.5, height_column = 3.0,
  dpi = 300, top_n = 20, volcano_label_top = 10, z_clamp = 2.5, nes_cap = 3.5,
  caption_wrap_column = 70, variants = c("print", "screen"),
  by_contrast_dir = "by_contrast", overview_dir = "_overview")

.DEFAULT_CONFIG_PATH <- "02_analysis/config/analysis_config.yaml"

## =====================================================================================
## 0. CONFIG — read the figures block once; everything else takes `config=` (no naked literals)
## =====================================================================================
load_figure_config <- function(path = NULL) {
  ## Load the project analysis_config.yaml and return the full parsed list. Pass it to every
  ## other function as `config=`. Reading once per script keeps these helpers side-effect-free.
  if (!requireNamespace("yaml", quietly = TRUE))
    stop("load_figure_config() needs the 'yaml' package to read analysis_config.yaml.")
  yaml::read_yaml(path %||% .DEFAULT_CONFIG_PATH)
}

.figures <- function(config) {
  ## Merge: project figures values over .FIG_DEFAULTS floors.
  f <- .FIG_DEFAULTS
  pf <- (config %||% list())$figures %||% list()
  modifyList(f, pf)
}

.fig_get <- function(config, key) {
  ## One figures-key lookup with default fallback (mirrors %||% on CONFIG$figures).
  .figures(config)[[key]] %||% .FIG_DEFAULTS[[key]]
}

.results_root <- function(config) {
  (config %||% list())$paths$results %||% "03_results/"
}

.master_root <- function(config) {
  (config %||% list())$paths$master %||% file.path(.results_root(config), "master")
}

.stage_dir <- function(config, stage, kind) {
  ## 03_results/<stage>/<figures|tables>/ — the one place the {figures,tables} subdir is named.
  stopifnot(kind %in% c("figures", "tables"))
  paths <- (config %||% list())$paths %||% list()
  subdir <- if (kind == "figures") (paths$stage_figures_subdir %||% "figures")
            else                   (paths$stage_tables_subdir %||% "tables")
  file.path(.results_root(config), stage, subdir)
}

## =====================================================================================
## 1. PATHS — per-contrast and cross-contrast; the ONLY sanctioned way to build these dirs
## =====================================================================================
contrast_path <- function(stage, contrast, kind = "figures", config = NULL) {
  ## Build + mkdir 03_results/<stage>/<kind>/by_contrast/<contrast>/; return the dir.
  ## kind in {"figures","tables"}. Contrast dir name MUST be the exact config contrast name.
  d <- file.path(.stage_dir(config, stage, kind), .fig_get(config, "by_contrast_dir"), contrast)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

overview_path <- function(stage, kind = "figures", config = NULL) {
  ## Build + mkdir 03_results/<stage>/<kind>/_overview/ (cross-contrast); return the dir.
  d <- file.path(.stage_dir(config, stage, kind), .fig_get(config, "overview_dir"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

.resolve_fig_dir <- function(stage, contrast, overview, config) {
  ## Pick the figures dir: by_contrast/<c>/ if `contrast` given, else _overview/ if `overview`.
  if (!is.null(contrast)) return(contrast_path(stage, contrast, "figures", config))
  if (isTRUE(overview))   return(overview_path(stage, "figures", config))
  d <- .stage_dir(config, stage, "figures")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

## =====================================================================================
## 2. THEME — the SINGLE style entry point. project_theme() is the R canonical name;
##    set_paper_style() is the thin same-named alias so a cross-language parity grep finds BOTH
##    contract names in this file (and the Python file likewise defines both).
## =====================================================================================
project_theme <- function(base_size = NULL, legend = TRUE, variant = "screen", config = NULL) {
  ## Return a ggplot2 theme built from the `figures:` config (LAZY ggplot2 load).
  ## Bold axis titles, no top/right spines, decluttered minor grid. Sizes come from the
  ## per-variant font tier (base_size for screen, base_size_column for print) and are enforced as
  ## FLOORS (clamped up + warned if a passed base_size is below). Use cairo on PDF export (see
  ## save_figure) so Unicode direction glyphs render. stop()s with context if ggplot2 is absent.
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("project_theme() needs ggplot2 (the plotting backend). Install ggplot2, or call only ",
         "the path/caption/table helpers (which need no backend).")
  f <- .figures(config)
  floor <- .variant_base_floor(variant, config)
  bs <- .enforce_floor(base_size %||% floor, floor, "base_size", variant)
  base <- ggplot2::theme_minimal(base_size = bs)
  base + ggplot2::theme(
    text             = ggplot2::element_text(size = bs),
    plot.title       = ggplot2::element_text(size = f$title_size, face = "bold"),
    axis.title       = ggplot2::element_text(size = f$axis_title_size, face = "bold"),  # bold axes
    axis.text        = ggplot2::element_text(size = f$axis_text_size),
    legend.text      = ggplot2::element_text(size = f$legend_text_size),
    legend.title     = ggplot2::element_text(size = f$legend_text_size),
    strip.text       = ggplot2::element_text(size = f$strip_size),
    legend.position  = if (isTRUE(legend)) "right" else "none",
    panel.grid.minor = ggplot2::element_blank(),               # declutter
    axis.line        = ggplot2::element_line(),                # keep bottom/left
    panel.border     = ggplot2::element_blank(),               # no top/right box (spine removal)
    plot.title.position = "plot")
}

set_paper_style <- function(base_size = NULL, legend = TRUE, variant = "screen", config = NULL) {
  ## Cross-language alias of project_theme() (the Python-side canonical name). Present so a parity
  ## grep finds `set_paper_style` in BOTH files. Returns the same ggplot2 theme object.
  project_theme(base_size = base_size, legend = legend, variant = variant, config = config)
}

## =====================================================================================
## 2b. Per-variant font floors — the legibility contract (clamp up, warn; never silently shrink)
## =====================================================================================
.variant_base_floor <- function(variant, config) {
  if (identical(variant, "print")) as.numeric(.fig_get(config, "base_size_column"))
  else                             as.numeric(.fig_get(config, "base_size"))
}

.enforce_floor <- function(value, floor, what, variant) {
  ## Clamp `value` up to `floor` if below it, warning. Legibility is never silently lost.
  value <- as.numeric(value)
  if (value < floor) {
    warning(sprintf("[figure-style] %s=%g below %s floor %g; clamping up to %g (legibility contract).",
                    what, value, variant, floor, floor), call. = FALSE)
    return(floor)
  }
  value
}

.variant_geometry <- function(variant, config) {
  ## c(width, height) inches for a variant: width/height (screen) / *_column (print).
  if (identical(variant, "print"))
    c(as.numeric(.fig_get(config, "width_column")), as.numeric(.fig_get(config, "height_column")))
  else
    c(as.numeric(.fig_get(config, "width")), as.numeric(.fig_get(config, "height")))
}

.normalize_variants <- function(variant, config) {
  vs <- if (identical(variant, "both")) (.fig_get(config, "variants") %||% c("print", "screen"))
        else variant
  vs <- as.character(unlist(vs))
  bad <- setdiff(vs, c("print", "screen"))
  if (length(bad)) stop(sprintf("variant must be print/screen/both, got %s", paste(bad, collapse = ",")))
  vs
}

## =====================================================================================
## 3. EXPORT — one plot object, two variant artifacts, from ONE config-driven code path
## =====================================================================================
save_figure <- function(plot, stage, name, variant = "both", contrast = NULL,
                        overview = FALSE, config = NULL) {
  ## Render ONE ggplot object to dual variants from one call (LAZY ggplot2/cairo load).
  ##   print  -> <stem>.print.pdf  : vector PDF via cairo_pdf, column geometry
  ##             (width_column x height_column), print font tier (base_size_column) so Unicode
  ##             glyphs render and text stays editable.
  ##   screen -> <stem>.screen.png : raster PNG @ dpi, screen geometry (width x height),
  ##             screen font tier (base_size).
  ## variant in {print, screen, both}. Output dir resolved via contrast_path()/overview_path()/
  ## the plain stage figures dir. Stale <name>*.{png,pdf} are purged first so the run owns its
  ## namespace. The SAME plot object is re-themed + re-sized per variant -> both files from one plot.
  ## Returns a named list (variant -> filepath).
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("save_figure() needs ggplot2 to render. On a backend-less box call the path/caption/",
         "table helpers instead, which need no plotting backend.")
  variants <- .normalize_variants(variant, config)
  out_dir  <- .resolve_fig_dir(stage, contrast, overview, config)
  purge_figures(stage, name, contrast = contrast, overview = overview, config = config)

  f <- .figures(config)
  written <- list()
  for (v in variants) {
    ext   <- if (identical(v, "print")) "pdf" else "png"
    out   <- file.path(out_dir, sprintf("%s.%s.%s", name, v, ext))
    geom  <- .variant_geometry(v, config)
    ## Re-theme the SAME plot object with this variant's font tier.
    styled <- tryCatch(plot + project_theme(variant = v, config = config),
                       error = function(e) plot)
    ## cairo_pdf for the print/PDF variant so Unicode direction glyphs (arrows, Delta) render;
    ## the default pdf() device substitutes "." for them. device = NULL lets ggsave infer for png.
    dev <- if (identical(v, "print") && isTRUE(capabilities()[["cairo"]]))
             grDevices::cairo_pdf else NULL
    ggplot2::ggsave(out, styled, width = geom[1], height = geom[2],
                    dpi = as.numeric(f$dpi), device = dev)
    written[[v]] <- out
    message(sprintf("  [figure-style] save_figure: %s (%gx%gin, base>= %gpt)",
                    basename(out), geom[1], geom[2], .variant_base_floor(v, config)))
  }
  invisible(written)
}

## =====================================================================================
## 3b. SERIES POST-STYLER — fix axis/legend for cross-panel comparability (was style_running_sum)
## =====================================================================================
style_series <- function(plot, ylim = NULL, config = NULL) {
  ## Pin a shared y-range + a fixed inside legend so a SERIES of figures stays comparable.
  ## Ported from 14839's style_running_sum: across a family of figures (e.g. one running-sum per
  ## database) a name-length-dependent outside legend silently resizes the plotting panel, so the
  ## curves stop being comparable. This clamps the y-axis to a fixed range via coord_cartesian
  ## (zoom, never drops data) and pins the legend INSIDE (zero layout width) so every figure in
  ## the series has identical panel proportions. `ylim` defaults to a symmetric clamp from config
  ## (z_clamp) when not given. LAZY ggplot2 load. Returns the styled plot.
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("style_series() needs ggplot2.")
  if (is.null(ylim)) {
    z <- .fig_get(config, "z_clamp")
    if (!is.null(z)) ylim <- c(-as.numeric(z), as.numeric(z))
  }
  styled <- tryCatch(plot + project_theme(config = config), error = function(e) plot)
  if (!is.null(ylim)) {
    ylim <- as.numeric(unlist(ylim))
    stopifnot(length(ylim) == 2, all(is.finite(ylim)))
    styled <- styled + ggplot2::coord_cartesian(ylim = ylim)
  }
  styled <- styled + ggplot2::theme(
    legend.position        = "inside",
    legend.position.inside = c(0.98, 0.98),
    legend.justification   = c(1, 1),
    legend.background      = ggplot2::element_rect(fill = "white", colour = "grey90"))
  styled
}

## =====================================================================================
## 4. PURGE — delete stale figures before a fresh write so a run OWNS its figure namespace
## =====================================================================================
purge_figures <- function(stage, prefix, contrast = NULL, overview = FALSE, config = NULL) {
  ## Delete <prefix>*.{png,pdf} under the resolved stage figures dir (DC semantics). Removes
  ## orphaned stems a fresh run no longer produces. Scoped by `prefix` so scripts sharing a
  ## figures/ dir don't clobber each other. Needs NO plotting backend (base R only). Returns count.
  d <- if (!is.null(contrast))
         file.path(.stage_dir(config, stage, "figures"), .fig_get(config, "by_contrast_dir"), contrast)
       else if (isTRUE(overview))
         file.path(.stage_dir(config, stage, "figures"), .fig_get(config, "overview_dir"))
       else
         .stage_dir(config, stage, "figures")
  n <- 0L
  if (dir.exists(d)) {
    ## Match by literal prefix + .png/.pdf suffix (no regex metacharacter hazards from `prefix`).
    all_files <- list.files(d, full.names = FALSE)
    hit <- startsWith(all_files, prefix) &
           (endsWith(all_files, ".png") | endsWith(all_files, ".pdf"))
    stale <- file.path(d, all_files[hit])
    if (length(stale)) { file.remove(stale); n <- length(stale) }
  }
  if (n > 0L)
    message(sprintf("  [figure-style] purge_figures('%s*'): removed %d stale file(s)", prefix, n))
  invisible(n)
}

## =====================================================================================
## 5. CAPTION — idempotent create/UPDATE of the sibling stage README.md (no plotting backend)
## =====================================================================================
write_caption <- function(stage, filename, finding, script, fn, config_kv, input,
                          how_to_read, config = NULL) {
  ## Idempotently write/replace ONE artifact's caption section in 03_results/<stage>/README.md.
  ## The section is a path-qualified heading (## figures/_overview/<file> etc.), a one-sentence
  ## FINDING, a mandatory **How to read** subsection (glyphs / sign convention / claim tier), and
  ## a `Script | Function | Config | Input` table. Re-running with the SAME `filename` REPLACES
  ## that file's section in place (idempotent — never duplicates). Prose wrapped to
  ## caption_wrap_column. Needs NO plotting backend (base R only).
  readme <- file.path(.results_root(config), stage, "README.md")
  dir.create(dirname(readme), recursive = TRUE, showWarnings = FALSE)
  wrap    <- as.integer(.fig_get(config, "caption_wrap_column"))
  heading <- sprintf("## %s", filename)
  section <- .render_caption_section(heading, finding, script, fn, config_kv, input, how_to_read, wrap)

  if (file.exists(readme)) {
    existing <- readLines(readme, warn = FALSE)
    new_text <- .replace_section(existing, heading, section)
  } else {
    new_text <- c(sprintf("# %s — artifact captions", stage), "", section)
  }
  writeLines(new_text, readme)
  invisible(readme)
}

.wrap <- function(text, width) {
  if (is.null(text) || !nzchar(as.character(text))) return("")
  paste(strwrap(as.character(text), width = width), collapse = "\n")
}

.render_caption_section <- function(heading, finding, script, fn, config_kv, input, how_to_read, wrap) {
  ## Render the canonical caption section block for one artifact (a character vector of lines).
  c(heading,
    "",
    .wrap(finding, wrap),
    "",
    paste0("**How to read:** ", .wrap(how_to_read, wrap)),
    "",
    "| Script | Function | Config | Input |",
    "|---|---|---|---|",
    sprintf("| `%s` | `%s` | `%s` | `%s` |", script, fn, config_kv, input),
    "")
}

.replace_section <- function(lines, heading, new_section) {
  ## Replace the `## <heading>` section in `lines` with `new_section`; append if absent. A section
  ## runs from its `## ` heading to the next `## `/`# ` heading or EOF. Exact heading match.
  out <- character(0)
  i <- 1L; n <- length(lines); replaced <- FALSE
  while (i <= n) {
    if (identical(trimws(lines[i]), trimws(heading))) {
      i <- i + 1L
      while (i <= n && !grepl("^## ", lines[i]) && !grepl("^# ", lines[i])) i <- i + 1L
      out <- c(out, new_section)
      replaced <- TRUE
      next
    }
    out <- c(out, lines[i]); i <- i + 1L
  }
  if (!replaced) out <- c(out, "", new_section)
  ## Trim trailing blank lines to a single terminator.
  while (length(out) && !nzchar(out[length(out)])) out <- out[-length(out)]
  c(out, "")
}

## =====================================================================================
## 6. OVERVIEW — the ATOMIC adjacency mechanism: figure + sibling table + caption in ONE call
## =====================================================================================
save_overview <- function(plot, stage, name, table, finding, script, fn, config_kv, input,
                          how_to_read, contrast = NULL, config = NULL) {
  ## Write a figure AND its same-stem source table AND its README caption in one call. The only
  ## sanctioned path for an overview/by-contrast figure: you cannot make the figure without its
  ## neighbor table + caption (source-table-adjacency + README-adjacency contracts, enforced
  ## mechanically). Writes:
  ##   figures/_overview/<name>.<variant>.<ext>   (via save_figure; or by_contrast/<c>/ if contrast)
  ##   tables/_overview/<name>.csv                 (data behind the figure; round_numeric_cols)
  ##   03_results/<stage>/README.md caption        (via write_caption, path-qualified, idempotent)
  ## Returns list(figures = <named list>, table = <path>, readme = <path>).
  overview <- is.null(contrast)
  figs <- save_figure(plot, stage, name, variant = "both", contrast = contrast,
                      overview = overview, config = config)

  if (!is.null(contrast)) {
    tdir    <- contrast_path(stage, contrast, "tables", config)
    rel_sub <- file.path("tables", .fig_get(config, "by_contrast_dir"), contrast)
    fig_rel <- file.path("figures", .fig_get(config, "by_contrast_dir"), contrast,
                         sprintf("%s.screen.png", name))
  } else {
    tdir    <- overview_path(stage, "tables", config)
    rel_sub <- file.path("tables", .fig_get(config, "overview_dir"))
    fig_rel <- file.path("figures", .fig_get(config, "overview_dir"),
                         sprintf("%s.screen.png", name))
  }
  table_path <- file.path(tdir, sprintf("%s.csv", name))
  if (!is.null(table)) utils::write.csv(round_numeric_cols(table), table_path, row.names = FALSE)

  readme <- write_caption(stage, fig_rel, finding = finding, script = script, fn = fn,
                          config_kv = config_kv, input = input, how_to_read = how_to_read,
                          config = config)
  message(sprintf("  [figure-style] save_overview: figure + %s/%s.csv + README caption",
                  rel_sub, name))
  invisible(list(figures = figs, table = table_path, readme = readme))
}

## =====================================================================================
## 7. MASTER TABLE — idempotent, byte-stable cross-stage accumulator append (no plotting backend)
## =====================================================================================
append_master_table <- function(df_or_rows, database, stage, name, config = NULL) {
  ## Idempotently append rows to 03_results/master/<name>.csv, deduped on the `database` column.
  ## Re-running for the SAME `database` REPLACES those rows (filter-out then append), so the master
  ## table is a stable accumulator across stages. round_numeric_cols(sig=9) is applied for
  ## byte-stability (re-runs produce identical files). Accepts a data.frame OR a list-of-lists
  ## (rows). Needs NO plotting backend. `stage` is recorded in a `stage` column for provenance.
  df <- .to_df(df_or_rows)
  if (is.null(df[["database"]]) || !"database" %in% names(df)) df[["database"]] <- database
  df[["database"]][is.na(df[["database"]]) | df[["database"]] == ""] <- database
  if (!"stage" %in% names(df)) df[["stage"]] <- stage
  df <- round_numeric_cols(df)

  out <- file.path(.master_root(config), sprintf("%s.csv", name))
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)

  if (file.exists(out)) {
    existing <- utils::read.csv(out, stringsAsFactors = FALSE, check.names = FALSE)
    kept <- existing[as.character(existing[["database"]]) != as.character(database), , drop = FALSE]
    ## Union the columns so old + new schemas coexist (stable order: existing first).
    all_cols <- union(names(kept), names(df))
    for (c0 in setdiff(all_cols, names(kept))) kept[[c0]]   <- NA
    for (c0 in setdiff(all_cols, names(df)))   df[[c0]]     <- NA
    combined <- rbind(kept[, all_cols, drop = FALSE], df[, all_cols, drop = FALSE])
  } else {
    combined <- df
  }
  utils::write.csv(combined, out, row.names = FALSE)
  invisible(out)
}

## =====================================================================================
## 8. ROUNDING — round every numeric column to `sig` significant digits (byte-stable outputs)
## =====================================================================================
round_numeric_cols <- function(df_or_rows, sig = 9) {
  ## Round all numeric columns to `sig` significant digits for byte-stable re-runs. Accepts a
  ## data.frame or a list-of-lists (rows); returns a data.frame. Non-numeric columns pass through.
  ## Needs NO plotting backend.
  df <- .to_df(df_or_rows)
  for (col in names(df)) {
    if (is.numeric(df[[col]])) df[[col]] <- signif(df[[col]], sig)
  }
  df
}

## =====================================================================================
## 9. DIRECTION CUE — map a sign to an unambiguous glyph/label (avoid a bare `*`); port of 14839
## =====================================================================================
direction_cue <- function(value) {
  ## Map a signed value to an unambiguous directional glyph/label (never a bare `*`). Positive ->
  ## "up" cue; negative -> "down" cue; zero / non-finite / non-numeric -> neutral. Glyphs are
  ## arrows + words so the cue is unambiguous in both color-blind and grayscale views.
  v <- suppressWarnings(as.numeric(value[1]))
  if (is.na(v)) return("· n/a")
  if (!is.finite(v) || v == 0) return("· n.s.")
  if (v > 0) "↑ up" else "↓ down"
}

## =====================================================================================
## Internal converter — keep the public API data.frame-or-rows agnostic
## =====================================================================================
.to_df <- function(df_or_rows) {
  ## Normalize input to a data.frame. Accepts a data.frame or a list of named row lists.
  if (is.null(df_or_rows)) return(data.frame())
  if (is.data.frame(df_or_rows)) return(df_or_rows)
  if (is.list(df_or_rows) && length(df_or_rows) && is.list(df_or_rows[[1]])) {
    ## list-of-rows -> data.frame
    cols <- unique(unlist(lapply(df_or_rows, names)))
    cells <- lapply(cols, function(k) sapply(df_or_rows, function(r) {
      v <- r[[k]]; if (is.null(v)) NA else v
    }))
    names(cells) <- cols
    return(as.data.frame(cells, stringsAsFactors = FALSE, check.names = FALSE))
  }
  as.data.frame(df_or_rows, stringsAsFactors = FALSE, check.names = FALSE)
}
