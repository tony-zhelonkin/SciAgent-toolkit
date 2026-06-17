# Strand-split TE counting — QC, the `-R CORE` witness, and the regression

The QC reference for the **three-channel** (s0 / sense / anti) TE counting this skill
performs. The canonical write-up — the strand-split invariant, the directional meter, the
axioms A1–A6, the warning-flag taxonomy, and the GREEN/RED gate thresholds — lives in the
toolkit **`docs/QC.md`** (a sibling of `METHODOLOGY.md` / `LIMITATIONS.md`). This page is the
skill-local pointer to it plus the regression that re-validates the behaviour inside the locked
`te-fc:2.0.2` container. Flags are cited **by name only**; their canonical definitions and
severities are owned by `docs/QC.md §5`.

Convention (reverse-stranded library, as in 14839-DM / 13036-DM): `sense = -s2 Assigned`,
`anti = -s1 Assigned`, `s0 = -s0`. The `-s2` pass is the correct library strand.

---

## The `-R CORE` regime-witness procedure (canonical: `docs/QC.md §4`)

Run featureCounts **3x on the SAME BAM** (`-s0` / `-s2` / `-s1`), **identical kernel**
(`-M -F SAF -p --countReadPairs -B -C`, no `-O`, no `--fraction`), adding `-R CORE`. Then:

- **WARNING (load-bearing foot-gun): join the three per-fragment `CORE` tables by a hash-join
  on the fragment ID, NOT a positional `paste`.** Read-ID order differs across threaded runs;
  a `paste` silently misaligns the three passes and fabricates the triple.
- Classify each fragment's `(s0, s2, s1)` status triple into a regime: conserved `(1,0)/(0,1)`
  → excess 0; antiparallel `(1,1)` → +2 (**AP**; double-presence); asymmetric `(2,1)/(1,2)` → +1
  (**M**; single reclaim = recovery); parallel-silent `(n,0)` → 0 (**P**; SILENT LOSS); bilateral
  `(>=2,>=2)` → 0 (silent loss); and the **A/A/A leak** (s0 Assigned but both stranded passes
  Assign a *different* feature) → +1 each.
- **Mechanism gate:** reconstruct excess and assert **per-fragment excess ∈ {0, +1, +2}**.
  Violation threshold `< 0.1%` (observed 0.00042% on 14839-DM-0019 — grade A). A value outside
  `{0,+1,+2}` is an instrument bug, not biology.

The witnessed ledger of the s0-ambiguous pile (locked, grade A): **AP 25.7% : M 6.1% :
P-silent 68.2%**; library-scale silent loss ≈ **9–10% of assigned TE signal**; bilateral
negligible (0.38%). The `excess/s0_Amb ≈ 0.78` ratio is a **directional meter, never an
estimator** of the split (grade C → A on the negative).

---

## Warning flags (named in `docs/QC.md §5`; cited here by name only)

Severities: 🔴 **hard** (violating it produces a wrong number / fabricated biology) ·
🟡 **watch** (conditional — must be checked).

| Flag name | Sev | Axiom | What it guards (one line) |
|---|:--:|:--:|---|
| `FLAG-SUM-CHANNELS` | 🔴 | A2 | Never use `(sense + anti)` as a quantity — double-counts the antiparallel + A/A/A piles. |
| `FLAG-S0-DENOM-ERV` | 🔴 | A3 | Never s0-denominate the LTR/ERV/Satellite/DNA tail — denominate on the **SENSE** channel. |
| `FLAG-KERNEL-MISMATCH` | 🔴 | A5 | Gene matrix is unique-only; every TE pass is `-M` — not like-for-like, no gene-vs-TE magnitude. |
| `FLAG-RESIDUAL-NOT-LOSS` | 🟡 | A1 | A zero/positive `s0 − sense − anti` residual does NOT certify no loss (silent loss is residual-invisible). |
| `FLAG-METER-NOT-ESTIMATOR` | 🟡 | A4 | `excess/s0_Amb ≈ 0.78` is a directional meter, never an estimator of the AP/M/P split. |
| `FLAG-SILENTLOSS-DESIGN` | 🟡 | A6 | Silent loss is sample-stable but `P_frac` co-varies with net strand-offset (r=−0.50) — verify against the design matrix. |
| `FLAG-ALIGNMENT-WEIGHTED` | 🔴 | A8 | The SENSE matrix is fragment-weighted under Random-One; if BAMs emit >1 locus/fragment, young high-NH families inflate up to ~meanNH-fold (definition in `docs/QC.md`). |
| `FLAG-DENOM-SCALE` | 🟡 | A9 | Per-sample vs library-summed denominators must never be mixed; label the scale on every ratio (definition in `docs/QC.md`). |

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
   locked ledger: `AP=169032`, `M_asymmetric=40362`, `P_silent=448490`, `bilateral=2478`,
   `excess_observed=506656`, `pred_excess` residual `−25.31%`, `resid_s0Amb_pct=0.0`,
   `truthtable_violations=192`.

The multi-GB `-R CORE` BAMs that produced the fixture are **NOT vendored**; the fixture's
`joint_counts` crosstab is the frozen witness they produced.

**Run `tests/strand_qc/run_regression.sh` to revalidate after any container or kernel change.**

---

## Grade tags

Consistent with this skill's `## Evidence & open questions` block: **A** peer-reviewed /
empirically witnessed · **B** tool default or single strong pipeline's recommendation ·
**C** sound mechanistic inference, not in a TE primary source · **D** folklore / mis-imported ·
**GAP** no adequate primary source — open. The mechanism gate and the AP/M/P ledger are grade
A (witnessed); the directional meter is graded C (directional-only); the silent-loss vs
design-matrix confound is the open **YELLOW** item (A6, grade B), pick-up-ready in
`featurecounts_results/DE_precheck/`.

See `docs/QC.md` for the full axioms, the GREEN/RED gates, the minOverlap sensitivity, and the
two open GAPs.
