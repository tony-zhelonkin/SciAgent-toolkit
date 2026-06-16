# How to QC a new TE dataset end-to-end

This is the **parameterized, operator-invoked** strand-split QC suite for the
`te-gene-featurecounts` skill. Point it at a NEW dataset's BAM dir + SAF + strand and it
**produces the witness** (measurement tables + GREEN/RED verdicts) — it does not assert frozen
values. The frozen regression that guards this suite's own logic lives separately in
`tests/strand_qc/` (it imports tool 05, the shared classifier).

Canonical doctrine — the strand-split invariant, the directional meter, axioms A1–A9, the
GREEN/RED gate thresholds, and every flag's definition — lives in the toolkit `docs/QC.md`. This
README is the runbook; flags are cited **by name**, deferring to `docs/QC.md`.

---

## One command

```bash
qc/run_qc.sh BAM_DIR SAF STRAND OUTDIR [--sample ID] [--design design.tsv] \
    [--matrix te_sense_s2.counts.txt] \
    [--s0-summary F --s2-summary F --s1-summary F --gene-summary F] \
    [--s0-matrix F --s2-matrix F --s1-matrix F] \
    [--young-regex '^(L1MdT|L1MdGf|L1MdA|IAPEz)'] [--image te-fc:2.0.2] [--threads 12]
```

- `BAM_DIR` — dir of star_salmon BAMs; one `--sample` is selected for the per-read witness (default:
  first BAM matching `--bam-glob '*.markdup.sorted.bam'`).
- `SAF` — the grouped TE SAF (`GeneID  Chr  Start  End  Strand`).
- `STRAND` — the dataset's library strand `0|1|2` (recorded; the witness always runs all three passes).
- `OUTDIR` — outputs (per-tool logs, tables, JSON, a consolidated verdict block at the end).

`run_qc.sh` chains the tools in dependency order (**04 feeds 02 / 05 / 06 / 07 / 08**), captures
each tool's verdict line, and prints a consolidated GREEN/RED block. **Container/bedtools tools skip
gracefully** (exit 0, clear message) when their dependency is absent — a SKIP is an absent
dependency, never a failure.

---

## Inputs by tool

| Need | Provide |
|---|---|
| per-read `-R CORE` witness (04, 01, 07, 08) | a BAM in `BAM_DIR` + the SAF + `te-fc:2.0.2` (docker) |
| fragment-weighting check (02) | `--matrix` (the s2 count matrix) + the s2 CORE (from 04) + the BAM + samtools |
| per-sample meter + residual (03) | `--s0-summary --s2-summary --s1-summary --gene-summary` (+ optional `--s*-matrix`) |
| geometry concordance (09) | the SAF + bedtools |
| DE_precheck (de_precheck) | `--design design.tsv` (+ the 03 outputs + `--gene-summary`) |

---

## Order, and what each tool does (with its GREEN/RED verdict + flags)

Run order in `run_qc.sh` (04 first because it is the engine):

1. **`tools/04_core_regime_witness.sh`** — runs featureCounts **3× on the same BAM** (`-s0/-s2/-s1`),
   identical kernel (`-M -F SAF -p --countReadPairs -B -C`, no `-O`, no `--fraction`) + `-R CORE`,
   in `te-fc:2.0.2`. **Hash-joins the three per-fragment CORE tables BY FRAGMENT ID — never a
   positional `paste`** (read-ID order differs across threaded runs) → `joint_counts_<sample>.json`
   (14-cell crosstab) + the `s0/s2/s1` CORE files. GREEN when the 3-pass join completes.
2. **`tools/05_classify_triples.py`** (the SHARED module — also imported by `tests/strand_qc/`) —
   re-derives the regime ledger from `joint_counts`: AP (antiparallel `+2`), M (asymmetric `+1`),
   P (parallel/bilateral silent `0`), the A/A/A leak. GREEN when the s0_Amb identity closes and
   per-fragment excess violations < 0.1%. The `2·AP+M` solve is a meter — `FLAG-METER-NOT-ESTIMATOR`,
   never report it as the split.
3. **`tools/06_closure_audit.py`** — proves the ledger **closes to the read**: identity (1)
   `excess = 2·AP + M + AAA − violations`, identity (2) `s0_Amb = AP + M + P + bilateral`, sense/anti
   sanity, and the full 14-cell contingency **sum vs N_fragments**. GREEN iff every identity residual
   < 0.5% and the table sums exactly to N — else a **hidden fifth bucket** exists, investigate.
   Prints the per-sample vs library-summed denominator table — `FLAG-DENOM-SCALE` (never mix scales;
   A/A/A is s0-Assigned, so `A/A/A ÷ s0_Amb` is not a real fraction).
4. **`tools/02_weighting_check.sh`** — is the SENSE integer matrix fragment- or alignment-weighted?
   Joins each s2-Assigned fragment to its true `NH` from the BAM; GREEN iff `matrix == nfrag` for
   every subfamily (fragment-weighted). Reports the young-set `meanNH` = **the counterfactual fold an
   alignment-weighted run WOULD have inflated young families by** (e.g. IAPEz ~16.6×). RED →
   `FLAG-ALIGNMENT-WEIGHTED`.
5. **`tools/01_random_one_check.sh`** — BAM-witness Random-One: restricted to `NH>1` fragments,
   GREEN iff **0** emit >1 locus. RED → `FLAG-ALIGNMENT-WEIGHTED`.
6. **`tools/07_silent_attribution.sh` + `.py`** — the `-O` target-revealer. `.sh` re-runs the s0 pass
   **kernel-identical + `-O`** (reheader BAM `chr1→1` with `--reheader` if the SAF is unprefixed) to
   emit candidate GeneID lists; `.py` identifies the silent set (fixed by the default-kernel
   statuses), looks up each silent fragment's target list, and reports the young share
   (contains-young = upper bound, exclusively-young = floor). **`-O` is a target-revealer ONLY — it
   does NOT redefine the silent set and is never used on the production matrix** (A7). GREEN at ≥99%
   `-O` coverage of the silent set.
7. **`tools/08_young_gate.py`** — young-autonomous **assignable-evidence** gate: conserved-fraction
   of the young set (`--young-regex`, default `^(L1MdT|L1MdGf|L1MdA|IAPEz)`). **GREEN ≥90% conserved**
   (non-conserved < 5%); **hard-RED at non-conserved ≥10%**; YELLOW between. Folds in 07's silent-loss
   young band. Silent loss is residual-invisible — `FLAG-RESIDUAL-NOT-LOSS`.
8. **`tools/09_geometry_overlap.py`** — read-INDEPENDENT SAF self-intersect (bedtools): antiparallel
   (double-presence) vs parallel (silent-loss) overlap geometry, per subfamily + young bin. Builds
   `geometry_vs_witness.tsv`: **concordance corroborates** the gate; **divergence flags** — resolve
   as read-density (expression × geometry) or a bin-definition artifact (the documented broad-L1Md
   2.5× → OLD-L1Md resolution) before trusting. Skips gracefully if bedtools is absent.
9. **`tools/03_strand_invariant.py`** — per-sample scalars + the `excess/s0_Amb` directional meter
   (mean ≈ 0.78, `FLAG-METER-NOT-ESTIMATOR`) + (optionally) the per-subfamily residual table. The
   residual must be one-sided: a **positive** subfamily residual signals a bug, not loss
   (`FLAG-RESIDUAL-NOT-LOSS`). Writes `per_sample_strand_qc.tsv` for the DE_precheck handoff.
10. **`de_precheck/`** — the 3-metric (silent-loss / multimapper-rate / strand-capture) ×
    design-permutation pre-check. Per-metric GREEN (flat tax, cancels in DE) or FLAG (covariate it /
    re-run gene `-M`-consistently). `FLAG-SILENTLOSS-DESIGN`. Needs `--design`.

---

## Reading the verdict block

`run_qc.sh` ends with a consolidated block listing each tool's last `VERDICT`/`SKIP` line and a
**SUITE VERDICT**: GREEN when no tool is RED (SKIPs = absent dependency, not failure), RED otherwise
with the count of flagged tools. The flags it carries forward — `FLAG-METER-NOT-ESTIMATOR`,
`FLAG-RESIDUAL-NOT-LOSS`, `FLAG-ALIGNMENT-WEIGHTED`, `FLAG-DENOM-SCALE`, `FLAG-SILENTLOSS-DESIGN` —
are defined in `docs/QC.md §5` (and named in `references/strand-split-qc.md`).

---

## Thresholds (summary)

| Check | GREEN | RED |
|---|---|---|
| Random-One (01) | 0 fragments emit >1 locus | any >1-locus fragment |
| Fragment-weighting (02) | `matrix == nfrag` every subfamily | `matrix ≠ nfrag` |
| Mechanism gate (05) | per-fragment excess violations < 0.1% | ≥0.1% (instrument bug) |
| Closure (06) | every identity residual < 0.5% AND sum == N | residual > 0.5% → hidden bucket |
| `-O` coverage (07) | ≥99% of silent set found a target list | < 99% |
| Young gate (08) | conserved ≥90%, non-conserved < 5% | non-conserved ≥10% |
| Geometry concordance (09) | witness within band of geometry | divergent → diagnose |
| DE_precheck (de_precheck) | no metric × design association | a fraction tilts with the design |

The mechanism gate and the AP/M/P ledger are grade A (witnessed); the directional meter is grade C
(directional-only). See `docs/QC.md` for the full grade scale and the open GAPs.
