# peak-atlas-framework — test scaffolds

These are **non-brittle test scaffolds** for the framework's shared R primitives.
They are deliberately *pseudocode-complete but not yet runnable*: every test is
guarded by `skip_scaffold()` so the suite is green today, and becomes a real
regression test once you implement the few TODO fixtures inside an R +
Bioconductor environment.

## Why "non-brittle"

The tests assert **invariants / contracts / properties** on **synthetic** inputs
— never byte-for-byte "golden" outputs on a frozen dataset. A golden-output test
breaks the moment a peak shifts by 1bp or a package reorders ties; an invariant
test ("the output is non-overlapping", "every survivor is an input peak", "fewer
peaks out than in", "support is computed before the merge") keeps passing across
versions and only fails when the *behaviour* is actually wrong. That is what you
want from a peak-calling pipeline where exact coordinates are data-dependent.

Principles baked in here:

- **Properties over snapshots** — `width == 501`, `no overlaps`, `out ⊆ in`,
  `monotone in threshold`, `partition of the input`, set-algebra identities.
- **Synthetic, deterministic fixtures** — small GRanges built with a fixed seed
  in `testthat/helper-fixtures.R`; no real data, no network, no large files.
- **Mock the externals** — MACS3, fragment files and BSgenome are stubbed so the
  pure logic is testable without heavyweight tools.
- **Encode the known footguns as regressions** — e.g. *support must be computed
  BEFORE the final merge* (the "bridging" artifact) and *iterative-overlap must
  be summit-faithful, not a midpoint*.

## Layout

```
tests/
  run_skill_tests.sh        # CI entry; skips (exit 0) if R/testthat absent
  testthat/
    helper-fixtures.R       # synthetic builders + skip_scaffold() + sourcing helpers
    test-iterative_overlap.R
    test-normalize_width.R
    test-support_voting.R
    test-frip.R             # needs a toy Signac object (TODO fixture)
    test-blacklist.R
    test-checks.R           # runs the CLI checks/ scripts, asserts exit codes
```

## Running

```bash
bash tests/run_skill_tests.sh          # from the skill dir
# or, directly:
Rscript -e "testthat::test_dir('tests/testthat')"
```

Until the scaffolds are implemented they report as **skipped**, not passed — that
is intentional.

## Developing a scaffold into a real test (in an R env)

1. Install deps: `testthat`, `GenomicRanges`, `IRanges`, `Matrix`, `rtracklayer`,
   and (for FRiP) `Signac`/`Seurat`.
2. Implement the TODO fixtures in `helper-fixtures.R` (e.g. `make_toy_atac_seurat`).
3. Remove the `skip_scaffold(...)` line from the test you are enabling.
4. Run `bash tests/run_skill_tests.sh` and confirm the invariant holds; if it
   fails, that is a real finding about the script, not the test.

## Invariant catalogue (script → what is guaranteed)

| Script | Key invariants tested |
|---|---|
| `iterative_overlap.R` | non-overlapping output; summit-faithful (`out ⊆ in`); one max-score winner per cluster; `|out| ≤ |in|`; idempotent; deterministic |
| `normalize_width.R` | all widths `== 2*extend+1`; centered on the summit; metadata preserved; summit-source precedence |
| `support_voting.R` | `n_strategies ∈ [0, #strategies]`; `adjusted = score*(1+boost*(n-1))`, monotone; **support computed before merge** (bridging regression) |
| `frip.R` | per-cell FRiP `∈ [0,1]`; monotone in the peak set; empty cells defined |
| `blacklist.R` | empty blacklist is a no-op; zero post-filter overlaps; `out ⊆ in`; missing BED degrades to empty |
| `checks/*.R` | valid atlas → exit 0; each injected violation → exit 1; FRiP retention gate exit codes |
