# Fragment paths and tabix — host vs container

## The problem

`Fragments(obj)` stores **absolute paths** from the analysis machine. Inside a container, those paths do not exist. The Coverage tab errors with `file not found` (or worse, returns an empty panel without erroring, depending on the Signac version).

Two layers must agree:

1. The `.rds` must hold paths that are valid **inside the container**.
2. The container's bind-mount must surface the actual files at exactly those paths.

## Strategy: rewrite to a stable container path

The skill standardises on `/data/fragments/` as the container path. Phase A.3 calls `UpdatePath(frag, new.path = "/data/fragments/<basename>")` for every Fragment object. The compose template bind-mounts the host fragments dir at `/data/fragments:ro`.

This means: the `.rds` becomes coupled to the deployment convention. To use it elsewhere, rewrite paths back. This is intentional — the alternative (a config layer that rewrites at app start) requires patching upstream `server.R`, which we are avoiding.

## Where do the fragment files come from?

`fragments.tsv.gz` files are produced by Cell Ranger ARC (or `arc-aggr`) — typically next to `outs/atac_fragments.tsv.gz` per sample. The companion `.tbi` file (tabix index) must be in the same directory and have the same basename + `.tbi`.

**Layout the bind-mount expects** (host side):

```
<FRAGMENTS_HOST_PATH>/
├── sample1_atac_fragments.tsv.gz
├── sample1_atac_fragments.tsv.gz.tbi
├── sample2_atac_fragments.tsv.gz
├── sample2_atac_fragments.tsv.gz.tbi
└── ...
```

Single-sample objects need exactly one pair. Multi-sample objects need one pair per sample, named to match what the `Fragment` object's `@path` (after Phase A.3 rewrite) refers to.

## How Phase A.3 verifies before save

```r
rewrite_fragment_paths(obj, container_dir = "/data/fragments",
                            host_dir      = "/scratch/.../fragments")
# Pseudocode:
# for each Fragment:
#   base <- basename(@path)
#   host_path <- file.path(host_dir, base)
#   if (!file.exists(host_path)) stop(...)
#   if (!file.exists(host_path + ".tbi")) stop(...)
#   if (@path != container_dir/base) UpdatePath(...)
```

The verification uses `host_dir` because Phase A runs on the host (or in dev), not inside the container. The path written into the `.rds` uses `container_dir`. If `host_dir == container_dir` (e.g., running Phase A on the deploy server with the bind-mount already in place), one path serves both.

## tabix index integrity

Signac uses `Rsamtools::scanTabix` for fragment random-access. The `.tbi` is **not** validated at object load — only at first `CoveragePlot` call. A broken or stale `.tbi` produces:

```
Error in Rsamtools::scanTabix(...): index file is older than data file
```

To rebuild from the host:

```bash
tabix -p bed sample1_atac_fragments.tsv.gz   # produces sample1_atac_fragments.tsv.gz.tbi
```

Or in R:

```r
Rsamtools::indexTabix("sample1_atac_fragments.tsv.gz", format = "bed")
```

The `validate_signac_rds.R` check tests both the `.tsv.gz` and the `.tbi` exist; it does not currently test that the index is newer than the data file. Add this if you regenerate fragment files frequently.

## Read-only mounting and barcode case

The compose template mounts the fragments dir `:ro`. Signac only reads — never writes. The `:ro` flag prevents accidental overwrites by a misbehaving Shiny session and helps document intent.

Note: `Fragments` carries a per-cell barcode list (`@cells`) that is matched against the file's `cellname` column (column 4 of the BED). Case-sensitive. If the file has `AAACGCT-1` but the object's `@cells` has `aaacgct-1`, the join silently produces zero coverage. This is upstream of the skill — fix at object construction time, not here.

## Why we don't move fragments inside the .rds

Some teams embed fragment files in the `.rds` via Signac's `CreateFragmentObject` with an in-memory representation. The skill does **not** support this for two reasons:

1. The `.rds` becomes huge (tens of GB) and slow to load at app startup.
2. Tabix random-access on disk is faster than scanning an in-memory representation.

If you have such an .rds, use Signac's `Fragments(obj) <- ` to convert to an external-file backing before Phase A.3.

## Multi-platform path normalisation

Windows → Linux container path conversion is not supported by the skill. Run Phase A on a POSIX host (Linux / macOS) so paths are stored as `/foo/bar/...`. If your analysis was on Windows, copy the `.rds` to a Linux box and re-run Phase A there.
