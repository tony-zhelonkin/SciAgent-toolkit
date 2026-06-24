# Known Issues — bulk-rnaseq-activity-inference

## `get_collectri()` / `get_progeny()` fail — OmniPath static-table fallback bug

**Symptom.** `decoupleR::get_collectri(organism=...)` and `decoupleR::get_progeny(organism=...)` error, typically with:
- `Error: argument is of length zero`
- an OmnipathR warning: *"Accessing `collectri` as a static table: ... only a backup plan for situations when our server or your computer is experiencing issues"*
- `OmnipathR::collectri()` failing on an internal `ncbi_tax_id` join.

**Root cause.** A version-compatibility bug between decoupleR and OmnipathR (observed with **decoupleR 2.16.0 + OmnipathR 3.18.4**, the latest Bioconductor-*release* versions as of 2026-06). OmnipathR decides the web service is unavailable and falls back to a bundled static-table loader that then fails to parse. This is **not necessarily a network outage** — confirm with:

```bash
curl -s -o /dev/null -w '%{http_code}\n' https://omnipathdb.org/queries/enzsub   # 200 => server reachable; it's the loader/version bug
```

Updating within the current Bioconductor release does **not** fix it (these are already the newest versions; there is no newer release to move into).

**Fix — build the networks locally from authoritative primary sources, cache as RDS.** None of these touch OmniPath:

```r
# PROGENy (14 pathways) — progeny Bioconductor data package
prog <- progeny::getModel("Mouse", top = 500)                 # genes x 14
net_progeny <- as.data.frame(prog) |>
  tibble::rownames_to_column("target") |>
  tidyr::pivot_longer(-target, names_to = "source", values_to = "weight") |>
  dplyr::filter(weight != 0)                                  # -> run_mlm(.mor = "weight")

# DoRothEA (mouse regulons) — dorothea Bioconductor data package
data("dorothea_mm", package = "dorothea")
net_dorothea <- dorothea_mm |>
  dplyr::filter(confidence %in% c("A", "B", "C")) |>
  dplyr::transmute(source = tf, target, mor)                  # -> run_ulm(.mor = "mor")

# CollecTRI — Zenodo (human) mapped to mouse via babelgene
ctri <- read.csv("https://zenodo.org/records/8192729/files/CollecTRI_regulons.csv")  # source,target,weight(±1)
orth <- babelgene::orthologs(unique(c(ctri$source, ctri$target)), species = "mouse", human = TRUE)
map  <- setNames(orth$symbol, orth$human_symbol)
tc   <- function(x){ m <- unname(map[x]); ifelse(is.na(m),
          paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x)))), m) }   # title-case fallback
net_collectri <- ctri |>
  dplyr::transmute(source = tc(source), target = tc(target), mor = weight) |>
  dplyr::distinct(source, target, .keep_all = TRUE)           # -> run_ulm(.mor = "mor")
```

Then call `run_ulm(mat, net, .source="source", .target="target", .mor="mor", minsize=5)` /
`run_mlm(..., .mor="weight")` exactly as in **Quick Start**, just with the local `net_*` objects in place of the `get_*` calls.

**Caveats.**
- The Zenodo→mouse CollecTRI is the published human network mapped by babelgene orthology (+ title-case fallback for unmapped symbols); it may differ slightly from canonical `get_collectri(organism="mouse")`. Note this where rankings are reported.
- `saveRDS()` the built networks so the analysis is reproducible and offline.
- Verify target overlap (`sum(rownames(mat) %in% net$target)` should be in the thousands — see the under-50-TFs pitfall in SKILL.md).

**Reference implementation.** `02_analysis/scripts/00c_prepare_networks.R` (STING-cGAS-GSE329522) builds and caches all three to `03_results/objects/net_{collectri_mouse, dorothea_mouse_ABC, progeny_mouse}.rds`.
