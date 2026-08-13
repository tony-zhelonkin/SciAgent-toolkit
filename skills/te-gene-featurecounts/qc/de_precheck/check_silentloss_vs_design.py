#!/usr/bin/env python3
"""
check_silentloss_vs_design.py — technical-fraction condition-linkage pre-check (3 metrics).

WHY THIS EXISTS
  The strand-split TE audit (see ../../_strand_invariant/SYNTHESIS_verdict.md) found three
  technical fractions that are each sample-STABLE under the usual conditions and so CANCEL in
  cross-sample DE — but each becomes a CONFOUND if a DE condition lines up with it. All three
  gate on the SAME thing: whether the fraction is associated with the design. So ONE design-matrix
  pass clears all three. This could not be run during the audit because the experimental design
  matrix (which sample is which condition) was not available. THIS SCRIPT runs it once you have it.

  The three technical fractions:
    (i)   silent-loss      — same-strand overlap "silent loss" (~9-10% of assigned TE signal,
                             residual-invisible). Sample-stable (P-fraction SD 0.0028 <
                             net-offset SD 0.0081) but its per-sample burden co-varies with the
                             net strand-offset (r = -0.50): if a condition rides that axis it
                             becomes a confound. Tracked by net_offset_frac / s0_amb_rate / P_pct_s0A.
    (ii)  multimapper-rate — the GENE pass ran WITHOUT -M (unique-only; ~112M multimappers dropped
                             library-wide) while the TE pass ran WITH -M (multimappers counted).
                             That kernel mismatch means a per-sample shift in the multimapper rate
                             tilts gene vs TE signal differently; if it lines up with a condition it
                             confounds any gene-vs-TE comparison. Tracked by multimapper_rate.
    (iii) strand-capture   — the net strand-offset itself (sense vs antisense capture). The primary
                             axis silent loss rides on; also a standalone QC (library prep / strand
                             specificity drift). Tracked by net_offset_frac (signed) and
                             strand_capture = sense/(sense+anti).

WHAT IT TESTS
  For each per-sample QC metric, whether it is ASSOCIATED with each design variable:
    - net_offset_frac  : (anti_Assigned - sense_Assigned)/(anti+sense). WITNESS-INDEPENDENT.
                         PRIMARY silent-loss axis (the r=-0.50 driver) AND the strand-capture axis.
    - s0_amb_rate      : s0_Amb / (s0_Assigned + s0_Amb). Overlap-ambiguity rate, witness-independent.
    - P_pct_s0A        : silent loss as % of assigned signal, from the SCALAR solve.
                         NOTE: the scalar solve UNDER-states absolute silent loss (witnessed ~9-10%
                         vs scalar ~7-8%); use this column only as a RELATIVE per-sample tracker.
    - multimapper_rate : gene Unassigned_MultiMapping / gene library total (per sample). The size
                         of the gene-no-`-M` vs TE-`-M` mismatch. WITNESS-INDEPENDENT.
    - strand_capture   : sense_Assigned/(sense_Assigned+anti_Assigned) = (1-net_offset_frac)/2.
                         Strand-capture as a direct fraction (the net strand-offset, sign-flipped).
  Categorical design var  -> one-way permutation ANOVA (F statistic, label-permutation p).
  Continuous design var   -> Spearman rho with permutation p.

VERDICT
  For each (metric x design_var): FLAG if permutation p < alpha (default 0.05) -> the metric
  differs by that design variable, i.e. that technical fraction may be confounded with it; MODEL IT
  (include the relevant axis as a covariate, or down-weight the affected tail). GREEN if no
  association. The PRIMARY decision for silent loss / strand-capture rests on net_offset_frac vs
  your CONDITION column; the multimapper check rests on multimapper_rate vs CONDITION.

USAGE
  python3 check_silentloss_vs_design.py design_matrix.tsv [--qc per_sample_strand_qc.tsv]
                                         [--alpha 0.05] [--nperm 20000] [--out results.tsv]
  design_matrix.tsv: TSV with a 'sample' column (matching 14839-DM-00NN) + >=1 design columns
  (e.g. condition, genotype, batch). See design_matrix.TEMPLATE.tsv. Deterministic (fixed seed).
"""
import sys, argparse, os
import numpy as np

SEED = 12345
# Three technical-fraction families, all cleared by one design-matrix pass:
#   silent-loss    -> net_offset_frac, s0_amb_rate, P_pct_s0A
#   multimapper    -> multimapper_rate
#   strand-capture -> net_offset_frac (primary), strand_capture
QC_METRICS = ["net_offset_frac", "s0_amb_rate", "P_pct_s0A", "multimapper_rate", "strand_capture"]
PRIMARY = "net_offset_frac"
# which technical fraction each metric belongs to (for the per-fraction roll-up at the end)
METRIC_FAMILY = {
    "net_offset_frac": "silent-loss/strand-capture",
    "s0_amb_rate": "silent-loss",
    "P_pct_s0A": "silent-loss",
    "multimapper_rate": "multimapper-rate",
    "strand_capture": "strand-capture",
}


def load_tsv(path):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = [ln.rstrip("\n").split("\t") for ln in f if ln.strip()]
    cols = {h: [r[i] for r in rows] for i, h in enumerate(header)}
    return header, cols


def to_float(v):
    try:
        return float(v)
    except Exception:
        return np.nan


def is_categorical(values):
    # treat as categorical if any value is non-numeric, or <=6 distinct numeric levels
    nums = []
    for v in values:
        try:
            nums.append(float(v))
        except Exception:
            return True
    return len(set(nums)) <= 6


def spearman_rho(x, y):
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    rx -= rx.mean(); ry -= ry.mean()
    denom = np.sqrt((rx**2).sum() * (ry**2).sum())
    return float((rx * ry).sum() / denom) if denom > 0 else 0.0


def anova_F(x, labels):
    grand = x.mean()
    ss_between = 0.0
    ss_within = 0.0
    groups = {}
    for v, l in zip(x, labels):
        groups.setdefault(l, []).append(v)
    k = len(groups); n = len(x)
    for l, vals in groups.items():
        vals = np.array(vals)
        ss_between += len(vals) * (vals.mean() - grand) ** 2
        ss_within += ((vals - vals.mean()) ** 2).sum()
    df_b = k - 1; df_w = n - k
    if df_w <= 0 or ss_within == 0:
        return float("inf") if ss_between > 0 else 0.0
    return float((ss_between / df_b) / (ss_within / df_w))


def perm_p_continuous(x, y, nperm, rng):
    obs = abs(spearman_rho(x, y))
    cnt = 1
    for _ in range(nperm):
        if abs(spearman_rho(x, rng.permutation(y))) >= obs - 1e-12:
            cnt += 1
    return cnt / (nperm + 1), spearman_rho(x, y)


def perm_p_categorical(x, labels, nperm, rng):
    labels = np.array(labels, dtype=object)
    obs = anova_F(x, labels)
    cnt = 1
    for _ in range(nperm):
        if anova_F(x, rng.permutation(labels)) >= obs - 1e-12:
            cnt += 1
    return cnt / (nperm + 1), obs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("design")
    ap.add_argument("--qc", default=os.path.join(os.path.dirname(__file__), "per_sample_strand_qc.tsv"))
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--nperm", type=int, default=20000)
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "silentloss_confound_result.tsv"))
    a = ap.parse_args()
    rng = np.random.default_rng(SEED)

    qh, qc = load_tsv(a.qc)
    # derive s0_amb_rate
    s0A = np.array([to_float(v) for v in qc["s0_Assigned"]])
    s0Amb = np.array([to_float(v) for v in qc["s0_Amb"]])
    qc["s0_amb_rate"] = list(s0Amb / (s0A + s0Amb))
    qc_by_sample = {s: i for i, s in enumerate(qc["sample"])}

    dh, dd = load_tsv(a.design)
    if "sample" not in dh:
        sys.exit("ERROR: design matrix must have a 'sample' column. See design_matrix.TEMPLATE.tsv")
    design_vars = [c for c in dh if c != "sample"]
    if not design_vars:
        sys.exit("ERROR: design matrix has no design columns besides 'sample'.")

    # align
    idx = []
    missing = []
    for s in dd["sample"]:
        if s in qc_by_sample:
            idx.append(qc_by_sample[s])
        else:
            missing.append(s)
    if missing:
        print(f"WARNING: {len(missing)} design samples not in QC table (ignored): {missing[:5]}...", file=sys.stderr)
    keep = [i for i, s in enumerate(dd["sample"]) if s in qc_by_sample]
    if len(keep) < 4:
        sys.exit(f"ERROR: only {len(keep)} samples matched; need >=4.")

    metric_vals = {}
    for m in QC_METRICS:
        col = qc[m] if m in qc else None
        if col is None:
            continue
        metric_vals[m] = np.array([to_float(col[qc_by_sample[dd['sample'][i]]]) for i in keep])

    missing_metrics = [m for m in QC_METRICS if m not in metric_vals]
    if missing_metrics:
        print(f"WARNING: QC table is missing metric column(s) {missing_metrics} -- "
              f"rebuild with build_per_sample_qc.py.", file=sys.stderr)

    results = []
    print(f"\n=== technical-fraction condition-linkage pre-check (3 metrics) "
          f"(n={len(keep)} samples, alpha={a.alpha}, nperm={a.nperm}) ===")
    print(f"    QC: {os.path.basename(a.qc)}   design: {os.path.basename(a.design)}")
    print(f"    fractions: silent-loss | multimapper-rate | strand-capture  (one design-matrix pass)\n")
    any_flag_primary = False
    flag_by_family = {}
    for dv in design_vars:
        dcol = [dd[dv][i] for i in keep]
        cat = is_categorical(dcol)
        for m, mv in metric_vals.items():
            if cat:
                p, stat = perm_p_categorical(mv, dcol, a.nperm, rng)
                test = "ANOVA-F(perm)"
            else:
                yv = np.array([to_float(v) for v in dcol])
                p, stat = perm_p_continuous(mv, yv, a.nperm, rng)
                test = "Spearman-rho(perm)"
            flag = "FLAG" if p < a.alpha else "ok"
            fam = METRIC_FAMILY.get(m, "")
            if dv.lower() in ("condition", "group", "treatment") and m == PRIMARY and flag == "FLAG":
                any_flag_primary = True
            if flag == "FLAG":
                flag_by_family.setdefault(fam, []).append((m, dv, p))
            note = fam + ("|PRIMARY" if m == PRIMARY else "")
            results.append((dv, m, test, round(stat, 4), round(p, 5), flag, note))
            mark = " <<< PRIMARY" if m == PRIMARY else ""
            print(f"  [{flag:4}] {dv:>14} x {m:<16} {test:<18} stat={stat:+.4f} p={p:.4f}  ({fam}){mark}")

    with open(a.out, "w") as f:
        f.write("design_var\tmetric\ttest\tstatistic\tperm_p\tverdict\tnote\n")
        for r in results:
            f.write("\t".join(str(x) for x in r) + "\n")

    print(f"\n  wrote {a.out}")
    print("\n=== OVERALL (per technical fraction) ===")
    FRACTIONS = [
        ("silent-loss", ["silent-loss", "silent-loss/strand-capture"],
         "include net_offset_frac (or a surrogate) as a DE covariate and/or down-weight the "
         "ERV/SINE/satellite overlap tail; re-examine those subfamilies' calls."),
        ("multimapper-rate", ["multimapper-rate"],
         "the gene-no-`-M` vs TE-`-M` mismatch tilts with this design -> do NOT apply gene size "
         "factors to TE unqualified; model multimapper_rate (or use a -M-consistent gene pass) "
         "for any gene-vs-TE comparison."),
        ("strand-capture", ["strand-capture", "silent-loss/strand-capture"],
         "strand specificity / net strand-offset tilts with this design -> add net_offset_frac "
         "(strand_capture) as a covariate; check library-prep batch."),
    ]
    flagged = [r for r in results if r[5] == "FLAG"]
    for label, fams, action in FRACTIONS:
        hits = [h for fam in fams for h in flag_by_family.get(fam, [])]
        # de-dup (net_offset_frac belongs to two families)
        hits = sorted(set(hits))
        if not hits:
            print(f"  GREEN  [{label:16}] no associated design variable -> flat, condition-independent tax (cancels in DE).")
        else:
            print(f"  FLAG   [{label:16}] {len(hits)} association(s):")
            for m, dv, p in hits:
                print(f"            - {m} differs by {dv} (p={p}).")
            print(f"            ACTION: {action}")
    print("\n=== SUMMARY ===")
    if not flagged:
        print("  GREEN: none of the 3 technical fractions is associated with any design variable.")
        print("  All behave as flat, condition-independent taxes -> standard cross-sample DE is safe.")
    else:
        print(f"  {len(flagged)} association(s) FLAGGED across the 3 fractions (see per-fraction block above).")
        if any_flag_primary:
            print("  NOTE: the PRIMARY test (net_offset_frac x condition) is flagged -> "
                  "treat silent loss / strand-capture as a confound.")
    print()


if __name__ == "__main__":
    main()
