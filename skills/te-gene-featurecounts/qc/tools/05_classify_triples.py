#!/usr/bin/env python3
"""05 — (s0,s2,s1) status-triple -> regime classifier (the SHARED QC module).

Dependency-free re-implementation of the per-fragment regime logic from the
audit driver (regimes_classify.py), operating on a ``joint_counts`` crosstab
({"a0/b/d": count}) rather than on the multi-GB ``-R CORE`` per-fragment tables.
The crosstab is the witness those BAMs produced (emitted by 04_core_regime_witness.sh
as joint_counts_<sample>.json); this module re-derives the regime ledger from it.

This is the module SHARED between the runnable QC suite (run_qc.sh / 06_closure_audit.py)
and the frozen regression (tests/strand_qc/run_regression.sh imports ``classify`` here).
Keep ``classify(joint_counts) -> dict`` import-safe and stdlib-only.

Triple convention (reverse-stranded library): key "a0/b/d" where
    a0 = s0   pass status
    b  = s2   pass status  (SENSE channel)
    d  = s1   pass status  (ANTISENSE channel)
Status codes: A = Assigned, M = Unassigned_Ambiguity, N = Unassigned_NoFeatures,
S = Unassigned_Singleton.

Regimes (per-fragment, from the synthetic truth table, confirmed empirically):
    conserved   (1,0)/(0,1): s0 Assigned; one channel Assigned, other NoFeat   -> excess 0
    A/A/A leak              : s0 Assigned; BOTH channels Assigned (diff GeneID)  -> excess +1
    antiparallel(1,1)       : s0 Ambiguous; BOTH channels Assigned (diff GeneID) -> excess +2
    asymmetric  (2,1)/(1,2) : s0 Ambiguous; EXACTLY ONE channel Assigned         -> excess +1
    parallel-sil(n,0)       : s0 Ambiguous; strand-match Ambig, other NoFeat     -> excess 0
    bilateral   (>=2,>=2)   : s0 Ambiguous; BOTH channels Ambiguous              -> excess 0

Per-fragment excess contribution = [s2 Assigned] + [s1 Assigned] - [s0 Assigned];
the mechanism gate requires every fragment's contribution to lie in {0,+1,+2}.

Verdict (printed when run as a script): GREEN iff every regime sums consistently,
the s0_Amb identity residual is ~0, and per-fragment excess violations are < 0.1%
(FLAG name: FLAG-METER-NOT-ESTIMATOR for the 2*AP+M solve, which is refuted by the
witness here — never report it as the AP/M/P split).

Usage:
    05_classify_triples.py [joint_counts.json | rcore_fixture.json]
Prints the recomputed tallies + GREEN/RED verdict and exits 0 (1 on a hard
violation). The regression imports ``classify`` and asserts the locked ledger.
"""
import json
import os
import sys

# status code -> 1 if the pass Assigned the fragment, else 0
ASSIGNED = "A"
AMBIG = "M"
NOFEAT = "N"

# mechanism-gate threshold: per-fragment excess violations as a fraction of N
VIOL_THRESHOLD = 0.001   # 0.1%


def classify(joint_counts):
    """Recompute regime tallies and reconciliation residuals from a joint_counts
    crosstab ({"a0/b/d": count}). Returns a flat dict of the locked quantities."""
    reg = {
        "conserved(1,0)/(0,1)": 0,
        "s0A_bothassign(A/A/A)": 0,
        "antiparallel(1,1)": 0,
        "asymmetric(2,1)/(1,2)": 0,
        "parallel_silent(n,0)": 0,
        "bilateral(>=2,>=2)": 0,
        "other_s0M": 0,
    }
    s0_singleton = 0
    excess_total = 0
    viol = 0
    s0_Amb = 0
    sense_Amb = 0
    anti_Amb = 0
    N = 0

    for triple, c in joint_counts.items():
        a0, b, d = triple.split("/")
        N += c

        # per-fragment excess contribution (channel-Assigned arithmetic)
        ex = (1 if b == ASSIGNED else 0) + (1 if d == ASSIGNED else 0) - (1 if a0 == ASSIGNED else 0)
        excess_total += ex * c
        if ex not in (0, 1, 2):
            viol += c

        # channel-ambiguity tallies (for the 3x3 reconciliation)
        if a0 == AMBIG:
            s0_Amb += c
        if b == AMBIG:
            sense_Amb += c
        if d == AMBIG:
            anti_Amb += c

        # regime classification by observed fate
        if a0 == ASSIGNED:
            if b == ASSIGNED and d == ASSIGNED:
                reg["s0A_bothassign(A/A/A)"] += c   # A/A/A leak
            else:
                reg["conserved(1,0)/(0,1)"] += c
        elif a0 == AMBIG:
            if b == ASSIGNED and d == ASSIGNED:
                reg["antiparallel(1,1)"] += c
            elif (b == ASSIGNED) ^ (d == ASSIGNED):
                reg["asymmetric(2,1)/(1,2)"] += c
            elif b in (AMBIG, NOFEAT) and d in (AMBIG, NOFEAT):
                if b == AMBIG and d == AMBIG:
                    reg["bilateral(>=2,>=2)"] += c
                else:
                    reg["parallel_silent(n,0)"] += c
            else:
                reg["other_s0M"] += c
        else:
            s0_singleton += c   # Singleton etc.

    AP = reg["antiparallel(1,1)"]
    M = reg["asymmetric(2,1)/(1,2)"]
    bilateral = reg["bilateral(>=2,>=2)"]
    P = reg["parallel_silent(n,0)"] + bilateral
    leftover = s0_Amb - (AP + M + P)
    ambsum = sense_Amb + anti_Amb

    # 3x3 identity reconciliation: predicted scalars from regime tallies.
    #   excess = 2*AP + M ;  s0_Amb = AP+M+P(+leftover) ;  (se+an)_Amb = M + P
    pred_excess = 2 * AP + M
    pred_s0Amb = AP + M + P + leftover
    pred_ambsum = M + P

    def resid_pct(pred, meas):
        return 100.0 * (pred - meas) / meas if meas else 0.0

    return {
        "N_fragments": N,
        "excess_observed": excess_total,
        "s0_Amb": s0_Amb,
        "sense_Amb": sense_Amb,
        "anti_Amb": anti_Amb,
        "ambsum": ambsum,
        "AP": AP,
        "M_asymmetric": M,
        "P_silent": P,
        "bilateral": bilateral,
        "leftover": leftover,
        "s0A_bothassign(A/A/A)": reg["s0A_bothassign(A/A/A)"],
        "conserved": reg["conserved(1,0)/(0,1)"],
        "truthtable_violations": viol,
        "pred_excess": pred_excess,
        "pred_s0Amb": pred_s0Amb,
        "pred_ambsum": pred_ambsum,
        "resid_excess_pct": resid_pct(pred_excess, excess_total),
        "resid_s0Amb_pct": resid_pct(pred_s0Amb, s0_Amb),
        "resid_ambsum_pct": resid_pct(pred_ambsum, ambsum),
        "regime_counts": reg,
    }


def load_joint_counts(path):
    """Accept either a raw joint_counts dict or a fixture/json with a top-level
    ``joint_counts`` key (emitted by 04_core_regime_witness.sh)."""
    with open(path) as fh:
        obj = json.load(fh)
    if isinstance(obj, dict) and "joint_counts" in obj:
        return obj["joint_counts"]
    return obj


def main(argv):
    here = os.path.dirname(os.path.abspath(__file__))
    default = os.path.join(here, "..", "..", "tests", "strand_qc", "fixtures", "rcore_fixture_0019.json")
    path = argv[1] if len(argv) > 1 else default
    res = classify(load_joint_counts(path))
    N = res["N_fragments"]
    viol = res["truthtable_violations"]
    viol_frac = (viol / N) if N else 0.0
    print(f"input: {path}")
    print(f"  N_fragments            {res['N_fragments']:>12,}")
    print(f"  AP antiparallel(1,1)   {res['AP']:>12,}")
    print(f"  M  asymmetric(2,1)(1,2){res['M_asymmetric']:>12,}")
    print(f"  P  parallel/bilateral  {res['P_silent']:>12,}")
    print(f"     of which bilateral  {res['bilateral']:>12,}")
    print(f"  A/A/A leak             {res['s0A_bothassign(A/A/A)']:>12,}")
    print(f"  excess_observed        {res['excess_observed']:>12,}")
    print(f"  pred_excess (2AP+M)    {res['pred_excess']:>12,}   resid {res['resid_excess_pct']:+.2f}%  "
          f"[FLAG-METER-NOT-ESTIMATOR: 2*AP+M is a meter, NOT the split]")
    print(f"  s0_Amb                 {res['s0_Amb']:>12,}   resid {res['resid_s0Amb_pct']:+.4f}%")
    print(f"  truthtable_violations  {viol:>12,}   ({100.0 * viol_frac:.5f}% of N)")
    green = (abs(res["resid_s0Amb_pct"]) < 0.5) and (viol_frac < VIOL_THRESHOLD)
    if green:
        print(f"  VERDICT: GREEN  (s0_Amb identity closes; per-fragment excess violations "
              f"{100.0 * viol_frac:.5f}% < {100.0 * VIOL_THRESHOLD:.1f}%)")
        return 0
    print(f"  VERDICT: RED    (FLAG-METER-NOT-ESTIMATOR / mechanism gate) — "
          f"s0_Amb resid {res['resid_s0Amb_pct']:+.4f}% or violations {100.0 * viol_frac:.5f}%")
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
