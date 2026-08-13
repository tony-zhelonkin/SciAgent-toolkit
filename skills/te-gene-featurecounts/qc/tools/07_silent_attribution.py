#!/usr/bin/env python3
"""07 — silent-loss set identify + -O target lookup + young-share attribution.

Generalized from green_close/identify_silent.py + attribute_silent.py + build_tables.py
(hard-coded sample/paths stripped; young-regex via flag). Dependency-free (stdlib only).
Operationalizes lesson 5 (the -O target-revealer / A7) + the silent-loss young attribution.

STEP 1 (identify): silent loss = the P bucket -- s0=M AND not Assigned in either stranded pass
  AND not both NoFeatures. i.e. s0=M, s2 in {M,N}, s1 in {M,N}, >=1 of {s2,s1} == M.
  These carry NO GeneID in the default-kernel CORE (target=NA), which is why the GREEN gate's
  silent-loss leg rested on read-INDEPENDENT geometry.
STEP 2 (lookup): for each silent fragment ID, look up its AMBIGUITY TARGET list (candidate
  GeneIDs it overlapped) from the -O CORE (07_silent_attribution.sh). -O is a TARGET-REVEALER
  ONLY -- it does NOT redefine the silent set (fixed by STEP 1) and is never used on the matrix.
STEP 3 (attribute): young-autonomous share of the silent pile, by --young-regex
  (default ^(L1MdT|L1MdGf|L1MdA|IAPEz), case-insensitive). Bands: contains-young (>=1 young
  target = upper bound) and exclusively-young (all targets young = floor).

Inputs:
    --s0-core --s2-core --s1-core   default-kernel CORE files (from 04), plain or .gz
    --sO-core                       the -O CORE (from 07_silent_attribution.sh), plain or .gz
    --young-regex                   default ^(L1MdT|L1MdGf|L1MdA|IAPEz)
    --outdir
Outputs (into --outdir):
    silent_targets.tsv         silent frag -> comma-sep candidate GeneID list
    young_silent_summary.tsv   contains-young / exclusively-young % of silent loss
    silent_class_table.tsv     per-class share of silent-loss target lists
    silent_witness.json        machine-readable

Threshold: 100% of identified silent fragments SHOULD find a multi-target record in the -O CORE
(coverage check). The young silent-loss share feeds 08's gate as the silent-loss band.

Usage:
    07_silent_attribution.py --s0-core s0.gz --s2-core s2.gz --s1-core s1.gz \\
        --sO-core sO.gz --outdir OUT [--young-regex RE]
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


def load_status(path):
    d = {}
    with opener(path) as fh:
        for line in fh:
            i = line.find("\t"); j = line.find("\t", i + 1)
            d[line[:i]] = ST.get(line[i + 1:j], "?")
    return d


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--s0-core", required=True)
    ap.add_argument("--s2-core", required=True)
    ap.add_argument("--s1-core", required=True)
    ap.add_argument("--sO-core", required=True, help="the -O CORE from 07_silent_attribution.sh")
    ap.add_argument("--young-regex", default=r"^(L1MdT|L1MdGf|L1MdA|IAPEz)")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args(argv[1:])
    os.makedirs(a.outdir, exist_ok=True)
    YOUNG = re.compile(a.young_regex, re.I)

    def subf(g):
        return g.split(":")[0]

    def cls(g):
        p = g.split(":")
        return p[2] if len(p) >= 3 else "UNK"

    def is_young(g):
        return bool(YOUNG.match(subf(g)))

    print("== 07 silent-loss identify + -O target lookup + young attribution ==")
    print(f"  young-regex: {a.young_regex}")

    # STEP 1: identify silent-loss fragment IDs
    s2 = load_status(a.s2_core)
    s1 = load_status(a.s1_core)
    silent = {}   # rid -> bucket
    n_par = n_bil = 0
    with opener(a.s0_core) as fh:
        for line in fh:
            i = line.find("\t"); j = line.find("\t", i + 1)
            rid = line[:i]; a0 = ST.get(line[i + 1:j], "?")
            if a0 != "M":
                continue
            b = s2.get(rid, "?"); d = s1.get(rid, "?")
            if b in ("M", "N") and d in ("M", "N") and not (b == "N" and d == "N"):
                if b == "M" and d == "M":
                    silent[rid] = "bilateral"; n_bil += 1
                else:
                    silent[rid] = "parallel_silent"; n_par += 1
    n_silent = len(silent)
    print(f"  silent-loss total = {n_silent:,}  (parallel_silent={n_par:,} bilateral={n_bil:,})")
    if n_silent == 0:
        print("  VERDICT [GREEN] no silent-loss fragments to attribute")
        return 0

    # STEP 2: look up the -O ambiguity target list for each silent fragment
    n_found = 0
    targets_path = os.path.join(a.outdir, "silent_targets.tsv")
    class_present = collections.Counter()
    first_class = collections.Counter()
    young_subf_hits = collections.Counter()
    ntargets_hist = collections.Counter()
    cnt = dict(gate=0, gate_excl=0, strict=0, strict_excl=0)

    def y_strict(g):
        s = subf(g)
        return s != "L1MDa" and is_young(g)

    with opener(a.sO_core) as fh, open(targets_path, "w") as out:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            rid = c[0]
            if rid not in silent:
                continue
            # Accept both the full -O CORE (id, status, n_targets, targets) and a reduced
            # (id, targets) form. The target list is the LAST tab field; reject NA/-1.
            tgts = c[-1] if len(c) >= 2 else ""
            if tgts in ("NA", "-1", "", rid):
                continue
            gids = tgts.split(",")
            out.write(f"{rid}\t{tgts}\n")
            n_found += 1
            ntargets_hist[len(gids)] += 1
            classes = set(cls(g) for g in gids)
            for cc in classes:
                class_present[cc] += 1
            first_class[cls(gids[0])] += 1
            g_gate = [g for g in gids if is_young(g)]
            if g_gate:
                cnt["gate"] += 1
                for g in g_gate:
                    young_subf_hits[subf(g)] += 1
                if all(is_young(g) for g in gids):
                    cnt["gate_excl"] += 1
            if any(y_strict(g) for g in gids):
                cnt["strict"] += 1
                if all(y_strict(g) for g in gids):
                    cnt["strict_excl"] += 1

    coverage = 100.0 * n_found / n_silent if n_silent else 0.0
    print(f"  -O target coverage: {n_found:,}/{n_silent:,} silent frags found a multi-target record "
          f"({coverage:.2f}%)")

    n = n_found or 1
    # STEP 3: tables
    with open(os.path.join(a.outdir, "silent_class_table.tsv"), "w") as fh:
        fh.write("class\tsilent_frags_overlapping_class\tpct_of_silent\n")
        for cc, k in class_present.most_common():
            fh.write(f"{cc}\t{k}\t{100.0 * k / n:.3f}\n")

    summ_path = os.path.join(a.outdir, "young_silent_summary.tsv")
    with open(summ_path, "w") as fh:
        fh.write("metric\tdefinition\tcount\tpct_of_silent_loss\n")
        fh.write(f"young_gate_contains\tgate set contains >=1 young target\t{cnt['gate']}\t{100.0 * cnt['gate'] / n:.3f}\n")
        fh.write(f"young_strict_contains\tstrict (drop older L1MDa) contains >=1\t{cnt['strict']}\t{100.0 * cnt['strict'] / n:.3f}\n")
        fh.write(f"young_gate_exclusive\tALL targets young (gate)\t{cnt['gate_excl']}\t{100.0 * cnt['gate_excl'] / n:.3f}\n")
        fh.write(f"young_strict_exclusive\tALL targets young (strict)\t{cnt['strict_excl']}\t{100.0 * cnt['strict_excl'] / n:.3f}\n")

    out = dict(silent_total=n_silent, silent_found=n_found, coverage_pct=coverage,
               parallel_silent=n_par, bilateral=n_bil, counts=cnt,
               pct={k: 100.0 * v / n for k, v in cnt.items()},
               class_present=dict(class_present),
               young_subf_hits=dict(young_subf_hits),
               ntargets_hist=dict(ntargets_hist))
    json.dump(out, open(os.path.join(a.outdir, "silent_witness.json"), "w"), indent=2)

    print(f"  young (gate) contains-young = {100.0 * cnt['gate'] / n:.3f}% of silent loss "
          f"(upper bound); exclusively-young = {100.0 * cnt['gate_excl'] / n:.3f}% (floor)")
    green = coverage >= 99.0
    tag = "GREEN" if green else "RED"
    print(f"  VERDICT [{tag}] -O target-revealer attributed the silent pile "
          f"(coverage {coverage:.2f}%; -O is a target-revealer ONLY, never on the matrix)")
    print(f"  wrote {targets_path}, {summ_path}, silent_class_table.tsv, silent_witness.json")
    return 0 if green else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
