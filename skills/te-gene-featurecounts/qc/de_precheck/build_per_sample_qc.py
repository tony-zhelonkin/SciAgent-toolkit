#!/usr/bin/env python3
"""
build_per_sample_qc.py -- (re)build/augment per_sample_strand_qc.tsv for the DE_precheck.

Reusable, dataset-agnostic version (relocated into the skill's qc/de_precheck). Ensures the
per-sample QC table carries the two technical-fraction columns the 3-metric condition-linkage
check needs, keeping all existing columns byte-stable:

  multimapper_rate : gene Unassigned_MultiMapping / gene library total, per sample.
                     The GENE pass ran WITHOUT -M (unique-only), so its Unassigned_MultiMapping
                     row carries the dropped multimappers; the TE pass ran WITH -M so its
                     MultiMapping row is all-zero. The gene-no-`-M` vs TE-`-M` kernel mismatch
                     (FLAG-KERNEL-MISMATCH); this column is its per-sample size.
                     Source: the gene .summary (--gene-summary); gene library total per sample
                     == TE N_fragments per sample (same library).

  strand_capture   : sense_Assigned / (sense_Assigned + anti_Assigned), per sample.
                     = (1 - net_offset_frac)/2 exactly. A sign-flipped view of the net
                     strand-offset axis. Computed from the QC table's own sense_Assigned /
                     anti_Assigned columns (written by qc/tools/03_strand_invariant.py), so this
                     tool needs no external per-sample table.

Parameterized: --qc-tsv (the per_sample_strand_qc.tsv to augment, in place) and --gene-summary
(the gene featureCounts .summary). Sample IDs are read from the QC table's `sample` column and
matched to the gene-summary BAM-basename columns (no dataset-specific ID regex). Deterministic,
dependency-free (stdlib only).

Usage:
    build_per_sample_qc.py --qc-tsv per_sample_strand_qc.tsv --gene-summary gene.summary
"""
import argparse
import os
import re
import sys

NEW_COLS = ["multimapper_rate", "strand_capture"]


def load_tsv(path):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = [ln.rstrip("\n").split("\t") for ln in f if ln.strip()]
    return header, rows


def basename_sample(col):
    """Reduce a gene-summary column (a BAM path) to its sample id."""
    b = os.path.basename(col)
    for suf in (".markdup.sorted.bam", ".sorted.bam", ".bam"):
        if b.endswith(suf):
            return b[: -len(suf)]
    return b


def gene_multimapper_rate(gene_summary):
    """per-sample gene multimapper-rate = Unassigned_MultiMapping / column total."""
    with open(gene_summary) as f:
        gh = f.readline().rstrip("\n").split("\t")
        rows = {}
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            rows[p[0]] = [int(x) for x in p[1:]]
    samples = [basename_sample(x) for x in gh[1:]]
    totals = [0] * len(samples)
    for vals in rows.values():
        for i, v in enumerate(vals):
            totals[i] += v
    mm = rows.get("Unassigned_MultiMapping", [0] * len(samples))
    return {s: (mm[i] / totals[i] if totals[i] else float("nan")) for i, s in enumerate(samples)}


def main(argv):
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser()
    ap.add_argument("--qc-tsv", default=os.path.join(here, "per_sample_strand_qc.tsv"),
                    help="per_sample_strand_qc.tsv to augment in place (default: alongside this script)")
    ap.add_argument("--gene-summary", required=True,
                    help="gene featureCounts .summary (for multimapper_rate)")
    a = ap.parse_args(argv[1:])

    if not os.path.exists(a.qc_tsv):
        sys.exit(f"ERROR: --qc-tsv not found: {a.qc_tsv} (run qc/tools/03_strand_invariant.py first)")

    qh, qr = load_tsv(a.qc_tsv)
    if "sample" not in qh:
        sys.exit("ERROR: QC table has no 'sample' column.")
    have_sense = "sense_Assigned" in qh and "anti_Assigned" in qh
    have_capture = "strand_capture" in qh

    mm = gene_multimapper_rate(a.gene_summary)

    # strip any prior run of the new columns so re-running is idempotent
    keep_idx = [i for i, c in enumerate(qh) if c not in NEW_COLS]
    base_h = [qh[i] for i in keep_idx]
    out_h = base_h + NEW_COLS
    lines = ["\t".join(out_h)]
    missing = 0
    for r in qr:
        d = dict(zip(qh, r))
        s = d["sample"]
        base = [r[i] for i in keep_idx]
        # strand_capture: prefer the table's existing column; else derive from sense/anti
        if have_capture and d.get("strand_capture", "") not in ("", "nan"):
            sc = float(d["strand_capture"])
        elif have_sense:
            se = float(d["sense_Assigned"]); an = float(d["anti_Assigned"])
            sc = se / (se + an) if (se + an) else float("nan")
        else:
            sc = float("nan")
        mr = mm.get(s, float("nan"))
        if s not in mm:
            missing += 1
        lines.append("\t".join(base + [f"{mr:.6f}", f"{sc:.6f}"]))
    with open(a.qc_tsv, "w") as f:
        f.write("\n".join(lines) + "\n")
    if missing:
        print(f"WARNING: {missing} QC sample(s) had no gene-summary column (multimapper_rate=nan)",
              file=sys.stderr)
    print(f"wrote {a.qc_tsv}  ({len(qr)} samples; +{NEW_COLS})")


if __name__ == "__main__":
    main(sys.argv)
