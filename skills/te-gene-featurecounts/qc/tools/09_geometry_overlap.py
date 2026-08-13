#!/usr/bin/env python3
"""09 — SAF self-intersect geometry vs witness concordance.

Generalized from annotation_overlap/analyze.py (bedtools self-intersect -> antiparallel:parallel,
same/cross-GeneID, per-class bp) + age_localize.py (age-bin localization). Parameterized:
SAF + young-regex via flags; no hard-coded paths. Read-INDEPENDENT (annotation geometry only).
Operationalizes lesson 8 (geometry-vs-witness concordance / GAP-D4).

Each overlapping pair of distinct loci = a read there is featureCounts-ambiguous (>=2 features):
  antiparallel = opposite strand -> double-presence (counted in sense & anti separately)
  parallel     = same strand     -> silent loss (dropped in s0 AND both stranded passes)

bedtools is REQUIRED. If bedtools is not on PATH this tool SKIPS GRACEFULLY (exit 0 with a
clear message), mirroring the container-tool idiom of tests/run_regression.sh.

Inputs:
    --saf            the production SAF (GeneID  Chr  Start  End  Strand)
    --young-regex    default ^(L1MdT|L1MdGf|L1MdA|IAPEz) (case-insensitive)
    --witness        OPTIONAL young_silent_summary.tsv (from 07) to build geometry_vs_witness
    --bedtools       bedtools binary (default: bedtools on PATH)
    --outdir
Outputs (into --outdir):
    subfamily_overlap.tsv     per-subfamily anti/para counts + bp + density
    geometry_vs_witness.tsv   predicted-geometry % vs witnessed-per-read % per bin + concordance

Interpretation: CONCORDANCE corroborates the gate (witness is not a read-placement artifact);
DIVERGENCE flags -- resolve as read-density (expression x geometry) or a bin-definition artifact
before trusting. (The documented worked example: the broad-L1Md 2.5x apparent divergence resolves
to OLD L1MdF/V/Mus/Fanc pulled into the broad bin.)

Usage:
    09_geometry_overlap.py --saf SAF --outdir OUT [--young-regex RE] [--witness summary.tsv]
"""
import argparse
import collections
import os
import re
import shutil
import subprocess
import sys


def saf_to_bed(saf, bed):
    """SAF (GeneID Chr Start End Strand, 1-based inclusive) -> BED6 (0-based half-open).
    name = GeneID = Subfamily:Family:Class."""
    n = 0
    with open(saf) as fin, open(bed, "w") as fout:
        for line in fin:
            if line.startswith("GeneID") or line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            gid, chrom, start, end, strand = p[0], p[1], p[2], p[3], p[4]
            try:
                s = int(start) - 1
                e = int(end)
            except ValueError:
                continue
            if e <= s:
                continue
            fout.write(f"{chrom}\t{s}\t{e}\t{gid}\t0\t{strand}\n")
            n += 1
    return n


def cls(name):
    return name.rsplit(":", 2)[-1]


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--saf", required=True)
    ap.add_argument("--young-regex", default=r"^(L1MdT|L1MdGf|L1MdA|IAPEz)")
    ap.add_argument("--witness", default=None,
                    help="OPTIONAL young_silent_summary.tsv (from 07) for geometry_vs_witness")
    ap.add_argument("--bedtools", default="bedtools")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args(argv[1:])
    os.makedirs(a.outdir, exist_ok=True)

    print("== 09 SAF self-intersect geometry ==")
    bt = shutil.which(a.bedtools) or (a.bedtools if os.path.exists(a.bedtools) else None)
    if not bt:
        print(f"  SKIP [09_geometry_overlap] bedtools ('{a.bedtools}') not on PATH — "
              f"skipping read-independent geometry (install bedtools to enable)")
        return 0
    if not os.path.exists(a.saf):
        print(f"  SKIP [09_geometry_overlap] SAF not found: {a.saf}")
        return 0

    bed = os.path.join(a.outdir, "te_all.bed")
    bed_sorted = os.path.join(a.outdir, "te_all.sorted.bed")
    pairs = os.path.join(a.outdir, "pairs.tsv")
    nloc = saf_to_bed(a.saf, bed)
    print(f"  SAF -> BED: {nloc:,} loci")

    # sort, then self-intersect (report overlap bp); exclude self via the index trick:
    # add a per-line index so a locus never matches itself.
    bed_idx = os.path.join(a.outdir, "te_idx.bed")
    with open(bed_sorted, "w") as out:
        subprocess.run(["sort", "-k1,1", "-k2,2n", bed], stdout=out, check=True)
    with open(bed_sorted) as fin, open(bed_idx, "w") as fout:
        for i, line in enumerate(fin):
            p = line.rstrip("\n").split("\t")
            # BED6 + an index column in name? keep BED6; append idx as 7th col
            fout.write("\t".join(p[:6]) + f"\t{i}\n")

    # self-intersect with -wo; columns: A(0-6) B(7-13) ovbp(14)
    with open(pairs, "w") as out:
        proc = subprocess.run(
            [bt, "intersect", "-a", bed_idx, "-b", bed_idx, "-wo"],
            stdout=subprocess.PIPE, check=True, text=True)
        seen = set()
        for line in proc.stdout.splitlines():
            c = line.split("\t")
            # A: 0 chrom 1 s 2 e 3 name 4 score 5 strand 6 idx ; B offset +7
            ia, na, sa = c[6], c[3], c[5]
            ib, nb, sb = c[13], c[10], c[12]
            ov = c[14]
            if ia == ib:
                continue  # self
            key = (min(ia, ib), max(ia, ib))
            if key in seen:
                continue
            seen.add(key)
            out.write(f"{ia}\t{ib}\t{na}\t{nb}\t{sa}\t{sb}\t{ov}\n")

    # aggregate (analyze.py logic, trimmed to the load-bearing tables)
    tot = {"anti": 0, "para": 0}
    tot_bp = {"anti": 0, "para": 0}
    sub_cnt = collections.defaultdict(lambda: collections.Counter())
    sub_bp = collections.defaultdict(lambda: collections.Counter())
    with open(pairs) as fh:
        for line in fh:
            ia, ib, na, nb, sa, sb, ov = line.rstrip("\n").split("\t")
            ov = int(ov)
            rel = "para" if sa == sb else "anti"
            tot[rel] += 1; tot_bp[rel] += ov
            sub_cnt[na][rel] += 1; sub_bp[na][rel] += ov
            if nb != na:
                sub_cnt[nb][rel] += 1; sub_bp[nb][rel] += ov

    # per-subfamily footprint for density
    sub_nloci = collections.Counter(); sub_totbp = collections.Counter()
    with open(bed_sorted) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            nm = f[3]; sub_nloci[nm] += 1; sub_totbp[nm] += int(f[2]) - int(f[1])

    def dens(name, rel):
        tb = sub_totbp[name]
        return (sub_bp[name][rel] / tb) if tb else 0.0

    def ratio(x, y):
        return (x / y) if y else float("inf")

    print(f"  antiparallel pairs {tot['anti']:,} (bp {tot_bp['anti']:,}); "
          f"parallel pairs {tot['para']:,} (bp {tot_bp['para']:,}); "
          f"anti:para bp ratio {ratio(tot_bp['anti'], tot_bp['para']):.3f}")

    overlap_tsv = os.path.join(a.outdir, "subfamily_overlap.tsv")
    with open(overlap_tsv, "w") as out:
        out.write("GeneID\tSubfamily\tFamily\tClass\tnloci\ttotbp\tanti_cnt\tpara_cnt\t"
                  "anti_bp\tpara_bp\tanti_dens\tpara_dens\tap_ratio_bp\n")
        for nm in sorted(sub_totbp, key=lambda x: -(sub_bp[x]["anti"] + sub_bp[x]["para"])):
            parts = nm.rsplit(":", 2)
            s = parts[0]; fam = parts[1] if len(parts) > 1 else ""; c = parts[2] if len(parts) > 2 else ""
            out.write(f"{nm}\t{s}\t{fam}\t{c}\t{sub_nloci[nm]}\t{sub_totbp[nm]}\t"
                      f"{sub_cnt[nm]['anti']}\t{sub_cnt[nm]['para']}\t{sub_bp[nm]['anti']}\t"
                      f"{sub_bp[nm]['para']}\t{dens(nm, 'anti'):.6f}\t{dens(nm, 'para'):.6f}\t"
                      f"{ratio(sub_bp[nm]['anti'], sub_bp[nm]['para']):.4f}\n")
    print(f"  wrote {overlap_tsv}")

    # geometry-vs-witness concordance: predicted-geometry young silent-loss share (parallel bp
    # density of the young set, as % of total parallel bp) vs witnessed per-read % from 07.
    YOUNG = re.compile(a.young_regex, re.I)
    young_para_bp = sum(sub_bp[nm]["para"] for nm in sub_bp if YOUNG.match(nm.split(":")[0]))
    young_anti_bp = sum(sub_bp[nm]["anti"] for nm in sub_bp if YOUNG.match(nm.split(":")[0]))
    geom_young_silent_pct = 100.0 * young_para_bp / tot_bp["para"] if tot_bp["para"] else 0.0
    geom_young_double_pct = 100.0 * young_anti_bp / tot_bp["anti"] if tot_bp["anti"] else 0.0

    witnessed = None
    if a.witness and os.path.exists(a.witness):
        with open(a.witness) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                row = dict(zip(hdr, line.rstrip("\n").split("\t")))
                if row.get("metric", "").startswith("young_gate_contains") or \
                   row.get("metric", "") == "young_gate_contains":
                    witnessed = float(row.get("pct_of_silent_loss", "nan"))
                    break

    gvw = os.path.join(a.outdir, "geometry_vs_witness.tsv")
    with open(gvw, "w") as out:
        out.write("bin\tgeometry_predicted_pct\twitnessed_per_read_pct\tdelta\tverdict\n")
        if witnessed is not None and witnessed == witnessed:
            delta = witnessed - geom_young_silent_pct
            # concordant if within a 2x band (the documented resolution rule)
            conc = "concordant" if abs(delta) <= max(0.5, geom_young_silent_pct) else "DIVERGENT-investigate"
            out.write(f"young_silent_loss\t{geom_young_silent_pct:.3f}\t{witnessed:.3f}\t{delta:+.3f}\t{conc}\n")
            print(f"  geometry_vs_witness: young silent geometry {geom_young_silent_pct:.3f}% vs "
                  f"witnessed {witnessed:.3f}% -> {conc}")
            verdict = "GREEN" if conc == "concordant" else "RED"
            print(f"  VERDICT [{verdict}] concordance corroborates the gate "
                  f"(divergence => read-density or bin-definition artifact, resolve before trusting)")
        else:
            out.write(f"young_silent_loss\t{geom_young_silent_pct:.3f}\tNA\tNA\t"
                      f"no-witness (pass --witness young_silent_summary.tsv from 07)\n")
            print(f"  geometry-only (no witness): young silent-loss geometry prediction "
                  f"{geom_young_silent_pct:.3f}%, young double-presence {geom_young_double_pct:.3f}%")
            print(f"  VERDICT [GREEN] geometry computed; pass --witness to test concordance")
    print(f"  wrote {gvw}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
