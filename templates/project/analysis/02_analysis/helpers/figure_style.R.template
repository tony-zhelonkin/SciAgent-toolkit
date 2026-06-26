## 02_analysis/helpers/figure_style.R — per-project figure-style shim.
## =====================================================================
## ONE import per viz script. Delegates to the SciAgent-toolkit contract
## lib (02_analysis/helpers/figure-style/figure_helpers.R, symlinked by
## `sciagent activate`) and loads the project analysis_config.yaml.
##
## Usage in any viz script:
##   source("02_analysis/helpers/figure_style.R")
##   p <- ggplot(...) + project_theme(config = FIG_CFG)
##   save_overview(p, "04_gsea", "name", table = df, ..., config = FIG_CFG)

## ---------------------------------------------------------------------------
## 1. Resolve the symlinked contract lib and source it (graceful fallback).
## ---------------------------------------------------------------------------
.HELPERS_LIB <- file.path(
    dirname(sys.frame(1L)$ofile %||% "02_analysis/helpers"),
    "figure-style",
    "figure_helpers.R"
)

# Portable %||% for use before the lib is sourced.
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

.FIGURE_STYLE_LOADED <- FALSE
if (file.exists(.HELPERS_LIB)) {
    source(.HELPERS_LIB)
    .FIGURE_STYLE_LOADED <- TRUE
} else {
    warning(
        "[figure_style] toolkit lib not found at: ", .HELPERS_LIB,
        "\nRun `sciagent activate` in this repo to link lib/figure-style/. ",
        "Falling back to minimal stubs.", call. = FALSE
    )
    # --- MINIMAL FALLBACK stubs (keep scripts from hard-crashing) -----------
    load_figure_config <- function(path = "02_analysis/config/analysis_config.yaml") {
        if (requireNamespace("yaml", quietly = TRUE)) yaml::read_yaml(path)
        else { warning("[figure_style] yaml unavailable; using empty config."); list() }
    }
    # `variant` accepted-but-ignored (drop-in compat with old call sites).
    project_theme <- function(base_size = NULL, legend = TRUE, variant = NULL,
                               config = NULL, ...) {
        if (!requireNamespace("ggplot2", quietly = TRUE))
            stop("[figure_style] ggplot2 required for project_theme().")
        ggplot2::theme_minimal(base_size = base_size %||% 14)
    }
    set_paper_style <- function(...) project_theme(...)
}

## ---------------------------------------------------------------------------
## 2. Load the project config once (FIG_CFG is the stable project-wide handle).
## ---------------------------------------------------------------------------
FIG_CFG <- tryCatch(
    load_figure_config("02_analysis/config/analysis_config.yaml"),
    error = function(e) {
        warning("[figure_style] Could not load analysis_config.yaml: ", conditionMessage(e),
                call. = FALSE)
        list()
    }
)
