#!/usr/bin/env Rscript
# render.R — render a decision-gate review notebook (.qmd) to two targets:
#   gfm  → GitHub-Flavored Markdown ([nb].md + [nb]_files/figure-gfm/*.png, committed):
#          renders inline in the GitHub file browser AND VS Code's Markdown preview.
#   html → self-contained HTML ([nb].html, figures embedded): quick browser view.
#
# Uses rmarkdown + system pandoc — no Quarto CLI required. Once the Quarto CLI is
# installed (scdock-r-dev >= v0.5.6), you may instead:  quarto render [nb].qmd --to gfm,html
#
# NB: GitHub strips base64 `data:` image URIs in markdown, so gfm MUST use committed PNG
# files (rmarkdown::github_document default) — do NOT set knitr upload.fun/image_uri for gfm.
#
# Each notebook lives in its OWN folder (02_analysis/notebooks/[name]/[name].qmd) so its
# rendered outputs (.md, .html, _files/) stage cleanly as one unit.
#
# Usage (from the compartment root):
#   Rscript 02_analysis/notebooks/render.R [notebook.qmd] [gfm|html|both]

# Headless raster: R's default bitmapType is Xlib (needs X11), so on a container png() fails
# and knitr silently falls back to svg. Force cairo here; the notebook setup chunk additionally
# sets dev="ragg_png". Both are headless and produce the committed PNGs the .md references.
options(bitmapType = "cairo")

args       <- commandArgs(trailingOnly = TRUE)
default_nb <- "02_analysis/notebooks/17_signature_review/17_signature_review.qmd"
nb   <- if (length(args) >= 1 && nzchar(args[1])) args[1] else default_nb
what <- if (length(args) >= 2 && nzchar(args[2])) args[2] else "both"
stopifnot("notebook not found" = file.exists(nb),
          "target must be gfm|html|both" = what %in% c("gfm", "html", "both"))

if (what %in% c("gfm", "both"))
  rmarkdown::render(nb, output_format = rmarkdown::github_document(html_preview = FALSE), quiet = TRUE)
if (what %in% c("html", "both"))
  rmarkdown::render(nb, output_format = rmarkdown::html_document(self_contained = TRUE,
                                                                toc = TRUE, toc_float = TRUE), quiet = TRUE)
message("rendered ", nb, "  ->  ", what)
