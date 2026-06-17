# DE_precheck — technical-fraction condition-linkage pre-check (reusable tool)

This is the **dataset-agnostic, relocated copy** of the DE_precheck tool, packaged inside the
`te-gene-featurecounts` skill (`qc/de_precheck/`). It runs the open verification check the
strand-split TE audit defers to DE time: **are the three per-sample technical fractions associated
with the experimental design?** Nothing here changes the count matrices; it gates how they may be
interpreted. The dataset-side instance (`featurecounts_results/DE_precheck/`) stays in place for
the run that produced it; this skill copy is the reusable tool.

**Parameterization.** `build_per_sample_qc.py` takes `--qc-tsv` (the `per_sample_strand_qc.tsv`
written by `qc/tools/03_strand_invariant.py`) and `--gene-summary` (the gene featureCounts
`.summary`) — no hard-coded dataset paths. `check_silentloss_vs_design.py` is already parameterized
(`design_matrix.tsv [--qc ...]`). From the suite, `qc/run_qc.sh ... --design design.tsv` wires both
automatically. Canonical doctrine + flag names live in the toolkit `docs/QC.md`.

---

## OPEN ITEM 1 — technical-fraction condition-linkage check (3 metrics)  🟡 YELLOW (must run before calling derepression)

**One design-matrix pass clears three technical fractions.** The audit found three per-sample
technical fractions that are each sample-STABLE (so each CANCELS in cross-sample DE *under the
usual conditions*) but each becomes a **confound** the moment a DE condition lines up with it.
All three gate on the same test — *is the fraction associated with the design?* — so a **single**
`design_matrix.tsv` and a **single** script run answers all three.

**The three fractions.**

1. **Silent-loss.** Same-strand overlap "silent loss" — reads dropped identically by the
   unstranded, sense and antisense passes and therefore invisible to the
   `unstranded − sense − antisense` residual — is the dominant fate of the strand-ambiguous pile:
   **roughly two-thirds of it, ≈ 9–10% of assigned TE signal** (per-read audit, sample-stable).
   Flat-tax on the variance criterion (its fraction-of-signal varies less, sample to sample, than
   the net strand-offset does), but its per-sample burden **co-varies with the net strand-offset**
   — so if a condition rides that axis it stops being flat. Tracked by `net_offset_frac` (PRIMARY),
   `s0_amb_rate`, `P_pct_s0A`.

2. **Multimapper-rate.** The **gene pass ran without `-M`** (unique-only; a large multimapper pile
   dropped library-wide) while the **TE pass ran with `-M`** (multimappers counted, Random-One).
   That kernel mismatch means a per-sample shift in the multimapper rate tilts gene vs TE signal
   differently; if it lines up with a condition it confounds any gene-vs-TE comparison or shared
   size-factor. Tracked by `multimapper_rate` (gene `Unassigned_MultiMapping` / gene library total,
   from `fc_genes/raw_fc_output/counts_matrix.txt.summary`; the gene library total per sample
   equals TE `N_fragments`, i.e. same library).

3. **Strand-capture.** The net strand-offset itself (sense vs antisense capture) — the primary axis
   silent loss rides on, and a standalone QC for strand-specificity / library-prep drift. Tracked by
   `net_offset_frac` (signed) and `strand_capture` = `sense_Assigned/(sense_Assigned+anti_Assigned)`
   = `(1 − net_offset_frac)/2`.

**Why it is still open.** None of the three could be tested during the audit because the **design
matrix (sample → condition) was not loaded.**

**How to clear all three (≈1 minute once you have the design).**
```bash
# per_sample_strand_qc.tsv comes from qc/tools/03_strand_invariant.py (writes it into the outdir).
# augment it with multimapper_rate + strand_capture (idempotent, stdlib only):
python3 build_per_sample_qc.py --qc-tsv per_sample_strand_qc.tsv --gene-summary gene.summary
cp design_matrix.TEMPLATE.tsv design_matrix.tsv     # fill in real condition/batch per sample
python3 check_silentloss_vs_design.py design_matrix.tsv --qc per_sample_strand_qc.tsv
```
- `per_sample_strand_qc.tsv` — the per-sample data the check uses (all 45 samples). Key columns:
  `net_offset_frac` (read-independent, the axis silent loss co-varies with — PRIMARY; also the
  strand-capture axis), `s0_amb_rate` (derived in-script), `P_pct_s0A` (scalar-solve silent-loss
  proxy; relative tracker only — it **understates** absolute silent loss, which the per-read audit
  puts higher), `multimapper_rate` (gene-no-`-M` vs TE-`-M` mismatch size), `strand_capture`
  (sense capture fraction).
- The script runs a dependency-free, deterministic permutation test of **each** QC metric against
  **each** design column (ANOVA for categorical, Spearman for continuous) and prints, **per technical
  fraction**, **GREEN** (no association → that fraction is a flat tax → safe) or **FLAG** with a
  fraction-specific action.
- **Decision rules.**
  - *Silent-loss / strand-capture:* the call rests on `net_offset_frac × condition`. If FLAGGED,
    treat silent loss as a confound — add `net_offset_frac` (or a surrogate) as a DE covariate and/or
    down-weight the ERV/SINE/satellite overlap tail, then re-examine those subfamilies' calls.
  - *Multimapper-rate:* the call rests on `multimapper_rate × condition`. If FLAGGED, do **not** apply
    gene size factors to the TE matrix unqualified — model `multimapper_rate` (or re-run the gene pass
    `-M`-consistently) for any gene-vs-TE comparison.
  - If all three GREEN, the technical fractions are flat taxes and cancel.

**Scope (not alarmist).** The young full-length autonomous-candidate subfamilies (L1Md_T/Gf/A,
IAPEz) are the LEAST affected — well past the conservation gate at read level (IAPEz-int essentially
immune) — the overlap burden sits on old SINE/MaLR-LTR (silent loss) and ERVK-LTR/satellite
(double-presence). This check is about whether these residual technical taxes tilt with your
specific design, not about the headline young-element signal.

---

## OPEN ITEM 2 — definitive young-silent-loss attribution  (GAP, optional / Rung-2)

Silently-lost reads carry **no GeneID** in `-R CORE` output, so the young-element GREEN above rests
on assignable evidence plus read-independent annotation geometry (only a small fraction of a percent
of silent loss localizes to the young L1Md / IAP-ETn gate set). A *definitive* per-subfamily
silent-loss attribution would need
`-R SAM`/BAM (or coordinate intersection) on a **genic-context-stratified SAF** (loci tagged
intron / intergenic). That is the Rung-2 step (see `docs/LIMITATIONS.md`); only pursue it if a young
autonomous subfamily becomes a load-bearing claim.

---

## Foot-guns these matrices inherit (full list in the per-matrix sidecars and `docs/QC.md`)

- **Never sum** `te_sense` + `te_antisense` — the antiparallel pile double-counts (100% cross-family).
- **Never use unstranded as the denominator** for the LTR/ERV/satellite/DNA tail — the unstranded
  pass drops their antiparallel reads to ambiguity (pathologically low denominator); denominate on
  the **sense** channel.
- **Gene vs TE are different kernels** — gene = unique-only (no `-M`, the multimapper pile dropped),
  TE = `-M` Random-One. Do not apply gene size factors to the TE matrix unqualified; no gene-vs-TE
  magnitude / "fraction of library."
