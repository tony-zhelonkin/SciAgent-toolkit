# Strand-split TE counting — QC, the `-R CORE` regime witness, and the regression

The QC reference for the **three-channel** (unstranded / sense / antisense) TE counting this
skill performs. The canonical write-up — the strand-split invariant, the directional meter, the
working principles, the warning-flag taxonomy, and the GREEN/RED gate thresholds — lives in the
toolkit **`docs/QC.md`** (a sibling of `METHODOLOGY.md` / `LIMITATIONS.md`). This page is the
skill-local pointer to it plus the regression that re-validates the behaviour inside the locked
`te-fc:2.0.2` container. Flags are cited **by name only**; their canonical definitions and
severities are owned by `docs/QC.md §5`.

**Convention (reverse-stranded dUTP library).** The featureCounts strand flags map to:
sense = `-s2` assigned, antisense = `-s1` assigned, unstranded = `-s0` assigned. The `-s2` pass
is the correct library strand. After this block the prose uses the words
**sense / antisense / unstranded**.

---

## The `-R CORE` regime-witness procedure (canonical: `docs/QC.md §4`)

Run featureCounts **3x on the SAME BAM** (`-s0` / `-s2` / `-s1`), **identical kernel**
(`-M -F SAF -p --countReadPairs -B -C`, no `-O`, no `--fraction`), adding `-R CORE`. Then:

- **WARNING (load-bearing foot-gun): join the three per-fragment `CORE` tables by a hash-join
  on the fragment ID, NOT a positional `paste`.** Read-ID order differs across threaded runs;
  a `paste` silently misaligns the three passes and fabricates the triple.
- Classify each fragment's `(unstranded, sense, antisense)` assigned-count triple into a regime.
  The three fates of a strand-ambiguous fragment: **AP — antiparallel** `(1,1)`, kept by both
  passes under *different* subfamilies (double-presence — what summing the channels would
  double-count) → excess +2; **M — asymmetric** `(2,1)/(1,2)`, where strand breaks the tie and
  the fragment is genuinely reclaimed → excess +1; **P — parallel / same-strand** `(n,0)`,
  separable by no strand flag and so silently lost in all three passes → excess 0. Plus
  conserved `(1,0)/(0,1)` → 0; bilateral `(>=2,>=2)` → 0 (silent loss); and the **A/A/A leak**
  (unstranded assigned, both stranded passes assign a *different* feature) → +1 each.
- **Mechanism gate:** reconstruct excess and assert **per-fragment excess ∈ {0, +1, +2}**.
  Violation threshold `< 0.1%` (observed well under that on the audited sample — grade A). A
  value outside `{0,+1,+2}` is an instrument bug, not biology.

Of the ambiguous pile, same-strand silent loss dominates — roughly two-thirds; antiparallel
double-presence about a quarter; genuine reclaim the small remainder. At library scale that is
silent loss ≈ **9–10% of assigned TE signal**; bilateral negligible. The
`excess/ambiguous ≈ 0.8` ratio is a **directional meter, never an estimator** of the split
(grade C; grade A on the negative). The exact per-sample tallies live in the regression fixture
(see below).

---

## Warning flags (named in `docs/QC.md §5`; cited here by name only)

Severities: 🔴 **hard** (violating it produces a wrong number / fabricated biology) ·
🟡 **watch** (conditional — must be checked).

| Flag name | Sev | What it guards (one line) |
|---|:--:|---|
| `FLAG-SUM-CHANNELS` | 🔴 | Never use `(sense + antisense)` as a quantity — double-counts the antiparallel + A/A/A piles. |
| `FLAG-S0-DENOM-ERV` | 🔴 | Never use unstranded as denominator for the LTR/ERV/Satellite/DNA tail — denominate on the **SENSE** channel. |
| `FLAG-KERNEL-MISMATCH` | 🔴 | Gene matrix is unique-only; every TE pass is `-M` — not like-for-like, no gene-vs-TE magnitude. |
| `FLAG-RESIDUAL-NOT-LOSS` | 🟡 | A zero/positive `unstranded − sense − antisense` residual does NOT certify no loss (silent loss is residual-invisible). |
| `FLAG-METER-NOT-ESTIMATOR` | 🟡 | `excess/ambiguous ≈ 0.8` is a directional meter, never an estimator of the AP/M/P split. |
| `FLAG-SILENTLOSS-DESIGN` | 🟡 | Silent loss is sample-stable but its share co-varies with the net strand-offset — verify against the design matrix. |
| `FLAG-ALIGNMENT-WEIGHTED` | 🔴 | The SENSE matrix is fragment-weighted under Random-One; if BAMs emit >1 locus/fragment, young high-NH families inflate up to ~meanNH-fold (definition in `docs/QC.md`). |
| `FLAG-DENOM-SCALE` | 🟡 | Per-sample vs library-summed denominators must never be mixed; label the scale on every ratio (definition in `docs/QC.md`). |

---

## The packaged regression

`tests/strand_qc/run_regression.sh` makes the two QC properties above executable. It skips
gracefully (exit 0) if docker or `te-fc:2.0.2` is absent, and otherwise asserts:

1. **Synthetic truth-table re-run.** For each frozen case in
   `fixtures/synthetic/sam/*.sam` (hand-crafted SAF + per-case SAMs), it runs featureCounts 3x
   (`-s0/-s2/-s1`) in `te-fc:2.0.2` with the recorded kernel (single-end `-M -F SAF`, or
   paired-end `-M -F SAF -p --countReadPairs -B -C`), parses the `.summary` + per-feature
   counts, and asserts value-equality to the frozen `fixtures/synthetic/expected_truth_table.tsv`.
   It then asserts every derived per-fragment excess ∈ {0,+1,+2} (the mechanism gate,
   end-to-end on the container).
2. **Regime-classifier validation.** `classify_triples.py` recomputes the regime tallies from
   `fixtures/rcore_fixture_0019.json["joint_counts"]` (a stripped, dependency-free
   re-implementation of the `-R CORE` triple → regime logic) and asserts it reproduces the
   locked tallies recorded in the fixture — the executable test contract: `AP=169032`,
   `M_asymmetric=40362`, `P_silent=448490`, `bilateral=2478`, `excess_observed=506656`,
   `pred_excess` residual `−25.31%`, `resid_s0Amb_pct=0.0`, `truthtable_violations=192`. (These
   are the regression's frozen assertions, not narrative magnitudes — leave them verbatim.)

The multi-GB `-R CORE` BAMs that produced the fixture are **NOT vendored**; the fixture's
`joint_counts` crosstab is the frozen audit they produced.

**Run `tests/strand_qc/run_regression.sh` to revalidate after any container or kernel change.**

---

## Grade tags

The single source of truth for the grade scale is this skill's `## Evidence & open questions`
block in `SKILL.md` (A / B / C / D / GAP); it is not re-defined here. In short: the mechanism
gate and the AP/M/P ledger are grade A (per-read audit); the directional meter is grade C
(directional-only); the silent-loss vs design-matrix confound is the open **YELLOW** item
(grade B), pick-up-ready in `featurecounts_results/DE_precheck/`.

See `docs/QC.md` for the full working principles, the GREEN/RED gates, the minOverlap
sensitivity, and the two open GAPs.
