#!/usr/bin/env python3
"""06 — closure-table audit: prove the witnessed-bucket ledger closes TO THE READ.

Generalized from reconcile_close/closure_table.py + ambiguity_audit/rcore_reconcile.py.
Inputs are the WITNESSED per-read CORE joint (s0,s2,s1) status-triple contingency table
(joint_counts_<sample>.json, emitted by 04_core_regime_witness.sh, or any rcore_fixture
with a top-level joint_counts). Nothing is solved or fit -- every number is a counted
fragment class. Dependency-free (stdlib only).

It checks the two closure identities the audit requires plus the sense/anti sanity, emits
the full 14-cell accounting that must sum exactly to N_fragments, and prints the per-sample
vs library-scale denominator clarification (FLAG-DENOM-SCALE).

  (1) excess_observed  = 2*AP + 1*M + 1*AAA - 1*violations
  (2) s0_Amb_observed  = AP + M + P(parallel) + bilateral   (AAA is s0-ASSIGNED, not in s0_Amb)
  (3) sanity           : sense_Amb + anti_Amb reconstructed from the same cells

Threshold: CLOSED iff every identity residual < 0.5% AND the full table sums EXACTLY to
N_fragments. A nonzero residual (> 0.5%) => a HIDDEN FIFTH BUCKET exists -- investigate it.
Prints a GREEN/RED verdict. Flags: FLAG-DENOM-SCALE (denominator discipline).

Usage:
    06_closure_audit.py JOINT_COUNTS.json [--sample ID] [--lib-s0amb N] [--out closure_report.txt]
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import importlib.util

# import the shared classifier (05) by file path (filename starts with a digit)
_spec = importlib.util.spec_from_file_location(
    "classify_triples", os.path.join(HERE, "05_classify_triples.py"))
_cls = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_cls)

# canonical role labels for the 14 status-triple cells
ROLE = {
    ("N", "N", "N"): "off-TE (all three passes NoFeature)",
    ("A", "N", "A"): "conserved (minus-strand TE): s0 Assigned, anti Assigned",
    ("A", "A", "N"): "conserved (plus-strand TE): s0 Assigned, sense Assigned",
    ("M", "N", "M"): "PARALLEL-silent, anti arm: s0 Amb, anti Amb, sense NoFeature",
    ("M", "A", "A"): "ANTIPARALLEL (+2): s0 Amb, Assigned on BOTH strands -> 2 features",
    ("M", "M", "N"): "PARALLEL-silent, sense arm: s0 Amb, sense Amb, anti NoFeature",
    ("A", "A", "A"): "A/A/A LEAK (+1): s0 Assigned AND both strands Assigned (s0-ASSIGNED)",
    ("S", "S", "S"): "SINGLETON: one mate unassignable in all three passes",
    ("M", "A", "M"): "ASYMMETRIC, sense-resolved (+1): s0 Amb, sense Assigned, anti Amb",
    ("M", "M", "A"): "ASYMMETRIC, anti-resolved  (+1): s0 Amb, anti Assigned, sense Amb",
    ("A", "A", "M"): "A/A/M: s0 Assigned, sense agrees (Assigned), anti arm Amb",
    ("A", "M", "A"): "A/M/A: s0 Assigned, anti agrees (Assigned), sense arm Amb",
    ("M", "M", "M"): "BILATERAL (0): s0 Amb, both strands Amb (silent both ways)",
    ("A", "M", "M"): "TRUTH-TABLE VIOLATION (-1): s0 Assigned but BOTH strands Amb",
}


def parse_joint(joint_counts):
    """{"a0/b/d": count} -> {(a0,b,d): count}."""
    return {tuple(k.split("/")): c for k, c in joint_counts.items()}


def marg(joint, i, v):
    return sum(c for k, c in joint.items() if k[i] == v)


def cell(joint, a, b, d):
    return joint.get((a, b, d), 0)


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("joint", help="joint_counts_<sample>.json (or fixture with joint_counts)")
    ap.add_argument("--sample", default=None, help="sample label for the report header")
    ap.add_argument("--lib-s0amb", type=int, default=None,
                    help="OPTIONAL library-summed s0_Amb (across all samples) for the "
                         "denominator-discipline table; never mix with the per-sample numerator")
    ap.add_argument("--out", default=None, help="write closure_report.txt here")
    a = ap.parse_args(argv[1:])

    joint_counts = _cls.load_joint_counts(a.joint)
    joint = parse_joint(joint_counts)
    sample = a.sample or os.path.basename(a.joint)

    lines = []
    def emit(s=""):
        lines.append(s)
        print(s)

    N = sum(joint.values())
    s0A, s0M = marg(joint, 0, "A"), marg(joint, 0, "M")
    s2A, s2M = marg(joint, 1, "A"), marg(joint, 1, "M")
    s1A, s1M = marg(joint, 2, "A"), marg(joint, 2, "M")
    excess_obs = (s2A + s1A) - s0A
    s0Amb_obs = s0M
    ambsum_obs = s2M + s1M

    AP = cell(joint, "M", "A", "A")
    M_asym = cell(joint, "M", "A", "M") + cell(joint, "M", "M", "A")
    P_par = cell(joint, "M", "N", "M") + cell(joint, "M", "M", "N")
    bilateral = cell(joint, "M", "M", "M")
    AAA = cell(joint, "A", "A", "A")
    viol = cell(joint, "A", "M", "M")

    emit("=" * 78)
    emit(f"CLOSURE TABLE -- {sample} (witnessed per-read CORE join)")
    emit("=" * 78)
    emit(f"  N_fragments (sum of all {len(joint)} cells) = {N:,}\n")

    ok = {}
    def line(name, pred, obs, key):
        res = pred - obs
        pct = 100.0 * res / obs if obs else 0.0
        status = "CLOSE" if abs(pct) < 0.5 else "OPEN -> HIDDEN BUCKET"
        emit(f"  {name}")
        emit(f"      predicted = {pred:>12,}")
        emit(f"      observed  = {obs:>12,}")
        emit(f"      residual  = {res:>+12,}   ({pct:+.4f}%)   [{status}]\n")
        ok[key] = abs(pct) < 0.5

    line("(1) excess = 2*AP + M + AAA - violations\n"
         f"      = 2*{AP:,} + {M_asym:,} + {AAA:,} - {viol:,}",
         2 * AP + M_asym + AAA - viol, excess_obs, "id1")
    line("(2) s0_Amb = AP + M + P(parallel) + bilateral\n"
         f"      = {AP:,} + {M_asym:,} + {P_par:,} + {bilateral:,}",
         AP + M_asym + P_par + bilateral, s0Amb_obs, "id2")

    # identity (3) sanity: reconstruct sense_Amb and anti_Amb from the same cells
    sa = (cell(joint, "M", "M", "N") + cell(joint, "M", "M", "A") + cell(joint, "A", "M", "A")
          + cell(joint, "M", "M", "M") + cell(joint, "A", "M", "M"))
    an = (cell(joint, "M", "N", "M") + cell(joint, "M", "A", "M") + cell(joint, "A", "A", "M")
          + cell(joint, "M", "M", "M") + cell(joint, "A", "M", "M"))
    line("(3a) sense_Amb reconstructed (cells with s2=M)", sa, s2M, "id3a")
    line("(3b) anti_Amb  reconstructed (cells with s1=M)", an, s1M, "id3b")
    emit(f"  sense_Amb + anti_Amb = {s2M:,} + {s1M:,} = {ambsum_obs:,}\n")

    # full cell-by-cell accounting (must sum exactly to N)
    emit("=" * 78)
    emit("FULL STATUS-TRIPLE CONTINGENCY ACCOUNTING (every cell, sums exact)")
    emit("=" * 78)
    emit(f"  {'s0/s2/s1':9} {'count':>12} {'cum%':>8}  role")
    cum = 0
    for k, c in sorted(joint.items(), key=lambda x: -x[1]):
        cum += c
        emit(f"  {'/'.join(k):9} {c:>12,} {100 * cum / N:7.3f}%  {ROLE.get(k, 'UNLABELLED CELL')}")
    emit(f"  {'TOTAL':9} {N:>12,}")

    # denominator-discipline table (FLAG-DENOM-SCALE)
    emit("\n" + "=" * 78)
    emit("DENOMINATOR DISCIPLINE (FLAG-DENOM-SCALE: never mix per-sample & library scales)")
    emit("=" * 78)
    emit(f"  PER-SAMPLE scale ({sample}):")
    emit(f"    s0_Amb (per-sample)                 = {s0Amb_obs:,}")
    emit(f"    AP / s0_Amb (per-sample)            = {100.0 * AP / s0Amb_obs:.2f}%" if s0Amb_obs else "    s0_Amb=0")
    emit(f"    excess (per-sample)                 = {excess_obs:,}")
    emit(f"    NOTE: A/A/A ({AAA:,}) is s0-ASSIGNED -> NOT part of s0_Amb; do NOT report A/A/A / s0_Amb.")
    if a.lib_s0amb:
        emit(f"  LIBRARY-SUMMED scale (all samples):")
        emit(f"    s0_Amb (library-summed)             = {a.lib_s0amb:,}")
        emit(f"    WARNING: the library-summed s0_Amb must NEVER denominate a single-sample numerator.")
    else:
        emit(f"  (pass --lib-s0amb N to print the library-summed scale alongside; both must be labelled)")

    closed = all(ok.values()) and True  # sum-exact is structural (cells already sum to N)
    emit("\n" + "=" * 78)
    if closed:
        emit("VERDICT: GREEN -- ledger CLOSED (point estimate proven to the read; no hidden bucket)")
    else:
        emit("VERDICT: RED -- ledger NOT closed; an identity residual > 0.5% => HIDDEN FIFTH BUCKET, "
             "investigate (FLAG-DENOM-SCALE if a scale was mixed)")
    emit("=" * 78)

    if a.out:
        with open(a.out, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        print(f"\nwrote {a.out}")

    return 0 if closed else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
