#!/usr/bin/env python3
"""08 — young-autonomous assignable-evidence gate + concordance read-in.

Generalized from ambiguity_audit/gate_young.py + green_close/build_tables.py. Parameterized:
young-regex via flag (default ^(L1MdT|L1MdGf|L1MdA|IAPEz)); CORE paths via flags; no hard-coded
sample. Dependency-free (stdlib only). Operationalizes lesson 8 (the GREEN gate) + the
silent-attribution / geometry concordance read-in.

Each fragment's regime is attributed to a SUBFAMILY via the Assigned-channel GeneID:
  conserved fragments carry the s0-assigned GeneID; antiparallel + asymmetric carry the
  stranded-assigned GeneID(s). Silent-loss fragments carry NO GeneID in default-kernel CORE
  (target=NA) -> NOT attributable here; the gate is therefore on ASSIGNABLE evidence and the
  silent-loss young share is read in separately from 07 (young_silent_summary.tsv) and 09
  (geometry_vs_witness.tsv) when present.

Inputs:
    --s0-core --s2-core --s1-core   the three default-kernel CORE per-fragment files
                                    (readname<TAB>status<TAB>n_targets<TAB>targets), plain or .gz.
    --young-regex                   default ^(L1MdT|L1MdGf|L1MdA|IAPEz) (case-insensitive)
    --silent-summary                OPTIONAL young_silent_summary.tsv (from 07) to fold in
    --outdir
Outputs (into --outdir):
    young_gate.json        conserved% on assignable + non-conserved bands + verdict
    young_subfamily_regimes.tsv  per-young-subfamily regime breakdown

Thresholds: young conserved-fraction GREEN >= 90% (non-conserved < 5%); hard-RED at
non-conserved >= 10%; YELLOW in between. Reports the silent-loss floor (exclusively-young)
and upper bound (contains-young) when 07's summary is provided.

Usage:
    08_young_gate.py --s0-core s0.gz --s2-core s2.gz --s1-core s1.gz --outdir OUT \\
        [--young-regex '^(L1MdT|L1MdGf|L1MdA|IAPEz)'] [--silent-summary young_silent_summary.tsv]
"""
import argparse
import collections
import gzip
import json
import os
import re
import sys

ST = {"Assigned": "A", "Unassigned_Ambiguity": "M",
      "Unassigned_NoFeatures": "N", "Unassigned_Singleton": "S"}


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_full(path):
    d = {}
    with opener(path) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            d[c[0]] = (ST.get(c[1], "?"), c[3] if len(c) > 3 else "NA")
    return d


def sub(t):
    if t in ("NA", "-1", ""):
        return None
    return t.split(",")[0].split(":")[0]


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--s0-core", required=True)
    ap.add_argument("--s2-core", required=True)
    ap.add_argument("--s1-core", required=True)
    ap.add_argument("--young-regex", default=r"^(L1MdT|L1MdGf|L1MdA|IAPEz)")
    ap.add_argument("--silent-summary", default=None,
                    help="OPTIONAL young_silent_summary.tsv from 07 to fold in the silent bands")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args(argv[1:])
    os.makedirs(a.outdir, exist_ok=True)
    YOUNG = re.compile(a.young_regex, re.I)

    def is_young(sf):
        return bool(sf and YOUNG.match(sf))

    print("== 08 young-autonomous assignable-evidence gate ==")
    print(f"  young-regex: {a.young_regex}")
    s0 = load_full(a.s0_core)
    s2 = load_full(a.s2_core)
    s1 = load_full(a.s1_core)

    persub = collections.defaultdict(collections.Counter)
    for rid, (a0, t0) in s0.items():
        b, tb = s2.get(rid, ("?", "NA"))
        d, td = s1.get(rid, ("?", "NA"))
        if a0 == "A":
            sf = sub(t0)
            if sf:
                persub[sf]["conserved" if not (b == "A" and d == "A") else "s0A_bothassign"] += 1
        elif a0 == "M":
            if b == "A" and d == "A":
                for t in (tb, td):
                    sf = sub(t)
                    if sf:
                        persub[sf]["antiparallel"] += 1
            elif (b == "A") ^ (d == "A"):
                sf = sub(tb if b == "A" else td)
                if sf:
                    persub[sf]["asymmetric"] += 1
            # silent loss -> NA target, not attributable here

    young = collections.Counter()
    other = collections.Counter()
    for sf, r in persub.items():
        tgt = young if is_young(sf) else other
        for k, v in r.items():
            tgt[k] += v

    def report(name, r):
        tot = sum(r.values()) or 1
        cons = r["conserved"]
        nonc = r["antiparallel"] + r["asymmetric"] + r["s0A_bothassign"]
        print(f"  -- {name} (assignable={tot:,}) conserved={100.0 * cons / tot:.2f}% "
              f"non-conserved={100.0 * nonc / tot:.2f}%")
        return dict(tot=tot, conserved=cons, antiparallel=r["antiparallel"],
                    asymmetric=r["asymmetric"], s0A_both=r["s0A_bothassign"],
                    conserved_frac=cons / tot, noncons_frac=nonc / tot)

    yr = report("YOUNG autonomous candidates", young)
    ar = report("ALL OTHER subfamilies", other)

    green = yr["conserved_frac"] >= 0.90 and yr["noncons_frac"] < 0.05
    red = yr["noncons_frac"] > 0.10
    gate = "GREEN" if green else ("RED" if red else "YELLOW")
    print(f"  young conserved_frac = {100 * yr['conserved_frac']:.2f}%  (gate: >=90%)")
    print(f"  young non-conserved  = {100 * yr['noncons_frac']:.2f}%  (gate: <5% green / >=10% RED)")
    print(f"  VERDICT [{gate}] young-gate on assignable evidence (FLAG-RESIDUAL-NOT-LOSS: "
          f"silent loss is residual-invisible -- see 07/09 for the silent young share)")

    # fold in 07's silent-loss bands if provided
    silent = None
    if a.silent_summary and os.path.exists(a.silent_summary):
        silent = {}
        with open(a.silent_summary) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                p = line.rstrip("\n").split("\t")
                row = dict(zip(hdr, p))
                silent[row.get("metric", p[0])] = row
        print("  silent-loss young bands (from 07):")
        for k in ("young_gate_contains", "young_strict_contains",
                  "young_gate_exclusive", "young_strict_exclusive"):
            if k in silent:
                print(f"    {k:24s} {silent[k].get('pct_of_silent_loss', '?')}% of silent loss")

    out = dict(young=yr, other=ar, gate=gate,
               young_regex=a.young_regex, silent_bands=silent)
    with open(os.path.join(a.outdir, "young_gate.json"), "w") as fh:
        json.dump(out, fh, indent=2)

    yl = [(sf, sum(r.values()), r) for sf, r in persub.items() if is_young(sf)]
    with open(os.path.join(a.outdir, "young_subfamily_regimes.tsv"), "w") as fh:
        fh.write("subfamily\ttot_assignable\tconserved\ts0A_both\tantiparallel\t"
                 "asymmetric\tconserved_pct\n")
        for sf, tot, r in sorted(yl, key=lambda x: -x[1]):
            fh.write(f"{sf}\t{tot}\t{r['conserved']}\t{r['s0A_bothassign']}\t"
                     f"{r['antiparallel']}\t{r['asymmetric']}\t{100.0 * r['conserved'] / tot:.2f}\n")
    print(f"  wrote {os.path.join(a.outdir, 'young_gate.json')} and young_subfamily_regimes.tsv")
    return 0 if not red else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
