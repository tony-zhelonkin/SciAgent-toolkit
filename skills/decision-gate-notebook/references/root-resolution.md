# Root resolution — why the setup chunk is written the way it is

A review notebook has to find the compartment root (where `02_analysis/`, `03_results/`,
`02_analysis/config/analysis_config.yaml` live) so `source()` and `readRDS()` resolve. It must
do this identically across three execution contexts that disagree about the working directory:

| Context | How it runs | What the wd is |
|---|---|---|
| `rmarkdown::render` / `quarto render` | knitr knits the doc | knitr sets it per-chunk from `opts_knit$root.dir` |
| VS Code "Run Cell" / R-extension send-to-terminal | lines sent to a radian session | the R terminal's launch dir (often the workspace root) |
| `radian` line-by-line | you paste chunks | wherever you launched `radian` |

## The bug the naive version has

```r
knitr::opts_knit$set(root.dir = .find_root())    # affects the RENDER pass only
here <- function(...) file.path(getwd(), ...)     # <-- resolves against the SESSION wd
```

`opts_knit$set(root.dir = ...)` changes the working directory **only for the knit/render
pass**. It does nothing to an interactive session's `getwd()`. So when you launch `radian`
in `02_analysis/notebooks/[nb]/` and paste the chunks, `here("02_analysis/helpers/figure_style.R")`
expands to `.../notebooks/[nb]/02_analysis/helpers/figure_style.R` — which doesn't exist, and
`source()` fails. You have to `setwd()` by hand to recover. That is the reported footgun.

## The fix

```r
.find_root <- function(sentinel = "02_analysis/config/analysis_config.yaml") {
  start <- tryCatch(dirname(knitr::current_input(dir = TRUE)), error = function(e) NA_character_)
  if (is.na(start) || !nzchar(start)) start <- getwd()
  d <- normalizePath(start, mustWork = FALSE)
  while (!file.exists(file.path(d, sentinel))) {
    parent <- dirname(d)
    if (identical(parent, d))
      stop("compartment root (", sentinel, ") not found above ", start,
           "\n  -> launch radian from inside the compartment, or render via notebooks/render.R")
    d <- parent
  }
  d
}
.root <- .find_root()
knitr::opts_knit$set(root.dir = .root)
if (interactive()) setwd(.root)
here <- function(...) file.path(.root, ...)
```

Three deliberate choices:

1. **Seed the search from the notebook's own location during render.**
   `knitr::current_input(dir = TRUE)` returns the absolute path of the file being knitted, so
   the sentinel walk starts *inside* the compartment regardless of the render process's wd. In
   an interactive session `current_input()` errors → the `tryCatch` falls back to `getwd()`.

2. **Bind `here()` to `.root`, not `getwd()`.** This is what actually fixes the footgun. Once
   `.root` is correct, `here()` is correct even if the session wd is somewhere else. This alone
   makes `source(here(...))` work line-by-line — the `setwd()` below is a convenience on top.

3. **`if (interactive()) setwd(.root)`.** So *bare* relative paths (`read.csv("03_results/...")`
   without `here()`) also work in a line-by-line session, matching what render does. Guarded by
   `interactive()` so it never fires during render (render already runs from `root.dir`), and it
   runs only *after* a successful `.find_root()`, so it never sets the wd to a wrong place.

## Why it is safe in umbrella containers

The sentinel is **compartment-relative** (`02_analysis/config/analysis_config.yaml`), never an
absolute `/workspaces/...` literal. So `.root` resolves to the same logical location whether the
compartment is mounted standalone (`/workspaces/[compartment]`) or nested under the umbrella
(`/workspaces/[umbrella]/[subproject]/[compartment]`). This is the same container-invariant,
resolve-from-a-marker discipline the umbrella `AGENTS.md` path note prescribes for scripts
(env var → `git rev-parse` toplevel → sentinel), applied to notebooks.

## The one case it cannot fix (by design)

If you launch `radian` **above** the compartment (e.g. at the repo/workspace root that contains
several compartments), the upward walk cannot know *which* compartment you meant — the sentinel
is a descendant, not an ancestor. `.find_root()` then `stop()`s with a clear message telling you
to launch from inside the compartment. This is correct: fail loud rather than guess. For
line-by-line work, start `radian` from the compartment root (or from the notebook's folder).
