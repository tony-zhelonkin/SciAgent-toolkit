#!/usr/bin/env python3
"""03 — strand-split invariant + per-sample excess meter + per-subfamily residual.

Generalized from build_residual.py (per-subfamily residual_frac), build_persample.py
(per-sample scalars + estimators), and p_stability.py (3x3 solve + P-stability +
net-offset SD). Parameterized: all inputs via flags, no hard-coded paths/sample list.
Dependency-free (stdlib only).

Operationalizes lessons 1 (strand-split invariant / A1) and 2 (excess/s0_Amb directional
meter / A4). Feeds the DE_precheck (qc/de_precheck): the per_sample_strand_qc.tsv it writes
carries net_offset_frac, s0_amb_rate, P_pct_s0A.

Inputs (scalars, REQUIRED for the per-sample meter):
    --s0 --s2 --s1 --gene   the four featureCounts .summary files (s0/sense-s2/anti-s1 TE +
                            gene). Sample columns are derived from the BAM basenames in the
                            Status header (so any sample set works).
Inputs (matrices, OPTIONAL, for the per-subfamily residual table):
    --s0-matrix --s2-matrix --s1-matrix   the three TE count matrices (GeneID + per-sample cols).
                            If given, writes residual_subfamily.tsv.

Outputs (into --outdir):
    per_sample_strand_qc.tsv  net_offset_frac, s0_amb_rate, P_pct_s0A, multimapper_rate*,
                              strand_capture, plus the raw scalars (* multimapper_rate filled
                              from the gene summary).
    residual_subfamily.tsv    (only if matrices given) per-subfamily s0/sense/anti sums +
                              residual_frac + per-sample SD.
    summary.json              excess/s0_Amb meter mean/SD (directional-only), the 3x3 solve
                              (flagged solve-not-witness), residual distribution.

Thresholds / interpretation:
    meter excess/s0_Amb ~ 0.78 is DIRECTIONAL ONLY (FLAG-METER-NOT-ESTIMATOR); never report
    the 3x3 solve (AP/M/P) as the split -- it omits the A/A/A leak (printed, labelled
    "solve, refuted by witness"). The per-subfamily residual must be one-sided: a POSITIVE
    subfamily residual_frac (> +5%) signals a bug, not loss (FLAG-RESIDUAL-NOT-LOSS context;
    silent loss is residual-invisible, so a clean/negative residual does NOT certify no loss).

Usage:
    03_strand_invariant.py --s0 s0.summary --s2 sense.summary --s1 anti.summary \\
        --gene gene.summary --outdir OUT [--s0-matrix M0 --s2-matrix M2 --s1-matrix M1]
"""
import argparse
import json
import os
import statistics as st
import sys


def parse_summary(path):
    """Return (samples, {status: [int per sample]}) from a featureCounts .summary."""
    rows = {}
    samples = None
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "Status":
                samples = [_basename_sample(p) for p in parts[1:]]
                continue
            if len(parts) < 2:
                continue
            rows[parts[0]] = [int(x) for x in parts[1:]]
    return samples, rows


def _basename_sample(p):
    b = os.path.basename(p)
    for suf in (".markdup.sorted.bam", ".sorted.bam", ".bam"):
        if b.endswith(suf):
            return b[: -len(suf)]
    return b


def col(d, key, n):
    return d.get(key, [0] * n)


def verdict_line(msg, green):
    tag = "GREEN" if green else "RED"
    print(f"  VERDICT [{tag}] {msg}")


# ----------------------------------------------------------------------------
# per-sample scalars + meter (build_persample.py + p_stability.py generalized)
# ----------------------------------------------------------------------------
def per_sample(s0, se, an, gene, samples, outdir):
    n = len(samples)
    s0A = col(s0, "Assigned", n); s0M = col(s0, "Unassigned_Ambiguity", n)
    seA = col(se, "Assigned", n); seM = col(se, "Unassigned_Ambiguity", n)
    anA = col(an, "Assigned", n); anM = col(an, "Unassigned_Ambiguity", n)
    geneMM = col(gene, "Unassigned_MultiMapping", n)
    # gene library total per sample (== TE N_fragments per sample)
    gene_tot = [0] * n
    for vals in gene.values():
        for i, v in enumerate(vals):
            gene_tot[i] += v

    records = []
    meter = []
    for i, s in enumerate(samples):
        excess = (seA[i] + anA[i]) - s0A[i]
        ambsum = seM[i] + anM[i]
        # 3x3 solve (refuted by the witness; printed for transparency, never the split)
        AP = s0M[i] - ambsum
        M = excess - 2 * AP
        P = ambsum - M
        ratio = excess / s0M[i] if s0M[i] else float("nan")
        net_off = (anA[i] - seA[i]) / (anA[i] + seA[i]) if (anA[i] + seA[i]) else float("nan")
        s0_amb_rate = s0M[i] / (s0A[i] + s0M[i]) if (s0A[i] + s0M[i]) else float("nan")
        strand_capture = seA[i] / (seA[i] + anA[i]) if (seA[i] + anA[i]) else float("nan")
        mm_rate = geneMM[i] / gene_tot[i] if gene_tot[i] else float("nan")
        P_pct_s0A = 100.0 * P / s0A[i] if s0A[i] else float("nan")
        if ratio == ratio:
            meter.append(ratio)
        records.append(dict(
            sample=s, s0_Assigned=s0A[i], s0_Amb=s0M[i],
            sense_Assigned=seA[i], sense_Amb=seM[i], anti_Assigned=anA[i], anti_Amb=anM[i],
            excess=excess, ambsum=ambsum, AP_solve=AP, M_solve=M, P_solve=P,
            ratio_excess_s0Amb=ratio, net_offset_frac=net_off, s0_amb_rate=s0_amb_rate,
            P_pct_s0A=P_pct_s0A, multimapper_rate=mm_rate, strand_capture=strand_capture))

    cols = ["sample", "s0_Assigned", "s0_Amb", "sense_Assigned", "sense_Amb",
            "anti_Assigned", "anti_Amb", "excess", "ambsum", "net_offset_frac",
            "s0_amb_rate", "P_pct_s0A", "multimapper_rate", "strand_capture",
            "ratio_excess_s0Amb", "AP_solve", "M_solve", "P_solve"]
    tsv = os.path.join(outdir, "per_sample_strand_qc.tsv")
    with open(tsv, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in records:
            fh.write("\t".join(
                (f"{r[c]:.6f}" if isinstance(r[c], float) else str(r[c])) for c in cols) + "\n")

    # global solve (labelled refuted-by-witness)
    GsM = sum(s0M); Gex = sum((seA[i] + anA[i] - s0A[i]) for i in range(n))
    Gamb = sum(seM) + sum(anM)
    GAP = GsM - Gamb; GM = Gex - 2 * GAP; GP = Gamb - GM

    def stats(vals):
        vals = [v for v in vals if v == v]
        if not vals:
            return dict(mean=float("nan"), sd=0.0, min=float("nan"), max=float("nan"))
        return dict(mean=st.mean(vals), sd=st.pstdev(vals), min=min(vals), max=max(vals))

    meter_stats = stats(meter)
    print(f"  N samples: {n}")
    print(f"  excess/s0_Amb DIRECTIONAL METER: mean={meter_stats['mean']:.4f} "
          f"sd={meter_stats['sd']:.4f} range[{meter_stats['min']:.4f},{meter_stats['max']:.4f}]")
    verdict_line("excess/s0_Amb is a directional meter, NOT an estimator (FLAG-METER-NOT-ESTIMATOR); "
                 "do not report the 3x3 solve as the AP/M/P split", True)
    print(f"  3x3 SOLVE (global, REFUTED BY WITNESS -- omits the A/A/A leak): "
          f"AP={GAP:,} M={GM:,} P={GP:,}  [solve, not the split]")
    print(f"  wrote {tsv}")
    return records, meter_stats, dict(AP_solve=GAP, M_solve=GM, P_solve=GP, s0_Amb=GsM,
                                      excess=Gex, ambsum=Gamb)


# ----------------------------------------------------------------------------
# per-subfamily residual (build_residual.py generalized)
# ----------------------------------------------------------------------------
def load_matrix(path):
    """Load a TE count matrix keyed by GeneID -> {sample: count}.

    Handles both layouts:
      * collapsed matrix: header "Geneid <s1> <s2> ..."   (1 leading id column)
      * raw featureCounts: optional "# Program:" comment line, then header
        "Geneid Chr Start End Strand Length <bam1> <bam2> ..."  (6 leading columns)
    """
    with open(path) as f:
        first = f.readline()
        while first.startswith("#"):
            first = f.readline()
        hdr = first.rstrip("\n").split("\t")
        # raw featureCounts has the fixed 6-column lead-in
        if len(hdr) >= 6 and hdr[:6] == ["Geneid", "Chr", "Start", "End", "Strand", "Length"]:
            lead = 6
        else:
            lead = 1
        cols = [_basename_sample(c) for c in hdr[lead:]]
        d = {}
        for line in f:
            p = line.rstrip("\n").split("\t")
            d[p[0]] = {cols[i]: int(float(p[lead + i])) for i in range(len(cols))}
    return d, cols


def residual_table(m0, m2, m1, outdir):
    d0, cols = load_matrix(m0)
    d2, _ = load_matrix(m2)
    d1, _ = load_matrix(m1)
    genes = [g for g in d0 if g in d2 and g in d1]
    rows = []
    for g in genes:
        parts = g.split(":")
        cls = parts[-1] if len(parts) >= 3 else "NA"
        fam = parts[-2] if len(parts) >= 3 else "NA"
        S0 = sum(d0[g].values()); SS = sum(d2[g].values()); SA = sum(d1[g].values())
        R = S0 - SS - SA
        RF = R / S0 if S0 > 0 else 0.0
        persamp = [(d0[g][s] - d2[g][s] - d1[g][s]) / d0[g][s] for s in cols if d0[g][s] > 0]
        sd = st.pstdev(persamp) if len(persamp) >= 2 else 0.0
        rows.append((g, cls, fam, S0, SS, SA, R, RF, sd, len(persamp)))

    tsv = os.path.join(outdir, "residual_subfamily.tsv")
    with open(tsv, "w") as f:
        f.write("Geneid\tClass\tFamily\ts0_sum\tsense_sum\tanti_sum\tresidual_sum\t"
                "residual_frac\tpersample_sd\tn_samples_nonzero\n")
        for r in sorted(rows, key=lambda x: -x[7]):
            g, cls, fam, S0, SS, SA, R, RF, sd, nn = r
            f.write(f"{g}\t{cls}\t{fam}\t{S0}\t{SS}\t{SA}\t{R}\t{RF:.6f}\t{sd:.6f}\t{nn}\n")

    nz = [r for r in rows if r[3] > 0]
    rfs = sorted(r[7] for r in nz)
    max_rf = rfs[-1] if rfs else 0.0
    leaky_pos = sum(1 for x in rfs if x >= 0.05)
    tot_s0 = sum(r[3] for r in rows); tot_r = sum(r[6] for r in rows)
    print(f"  residual: {len(nz)} subfamilies with s0>0; max residual_frac={max_rf:+.4f}; "
          f"net residual_frac={tot_r / tot_s0:+.4f}" if tot_s0 else "  residual: no signal")
    green = max_rf < 0.05
    verdict_line(f"per-subfamily residual one-sided (max residual_frac {max_rf:+.4f} < +0.05); "
                 f"{leaky_pos} subfamilies with positive residual >= +5%"
                 + ("" if green else " -> FLAG-RESIDUAL-NOT-LOSS: positive residual signals a bug"),
                 green)
    print(f"  wrote {tsv}")
    return dict(max_residual_frac=max_rf, leaky_pos=leaky_pos,
                net_residual_frac=(tot_r / tot_s0) if tot_s0 else 0.0,
                n_nonzero=len(nz))


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--s0", required=True, help="s0 (unstranded) TE .summary")
    ap.add_argument("--s2", required=True, help="sense (-s2) TE .summary")
    ap.add_argument("--s1", required=True, help="antisense (-s1) TE .summary")
    ap.add_argument("--gene", required=True, help="gene .summary (for multimapper_rate)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--s0-matrix", default=None)
    ap.add_argument("--s2-matrix", default=None)
    ap.add_argument("--s1-matrix", default=None)
    a = ap.parse_args(argv[1:])
    os.makedirs(a.outdir, exist_ok=True)

    print("== 03 strand-invariant + per-sample meter ==")
    samples, s0 = parse_summary(a.s0)
    _, se = parse_summary(a.s2)
    _, an = parse_summary(a.s1)
    _, gene = parse_summary(a.gene)
    records, meter_stats, gsolve = per_sample(s0, se, an, gene, samples, a.outdir)

    resid_summary = None
    if a.s0_matrix and a.s2_matrix and a.s1_matrix:
        print("== per-subfamily residual ==")
        resid_summary = residual_table(a.s0_matrix, a.s2_matrix, a.s1_matrix, a.outdir)
    else:
        print("  (per-subfamily residual skipped: pass --s0-matrix/--s2-matrix/--s1-matrix to enable)")

    summary = dict(
        n_samples=len(samples),
        meter_excess_over_s0Amb=meter_stats,
        global_3x3_solve_REFUTED_BY_WITNESS=gsolve,
        residual=resid_summary,
        flags=["FLAG-METER-NOT-ESTIMATOR", "FLAG-RESIDUAL-NOT-LOSS"],
    )
    with open(os.path.join(a.outdir, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"  wrote {os.path.join(a.outdir, 'summary.json')}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
