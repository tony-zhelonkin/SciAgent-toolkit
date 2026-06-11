#!/usr/bin/env bash
# tests/run_skill_tests.sh — smoke tests for the te-gene-featurecounts skill.
#
# Verifies the locked te-fc:2.0.2 container:
#   (a) image exists,
#   (b) featureCounts reports v2.0.2,
#   (c) a synthetic mini run produces an INTEGER count matrix of expected shape.
#
# Skips gracefully (exit 0) when docker or the image is unavailable, so the
# toolkit-wide tests/run-all.sh stays green on machines without docker.
set -u

IMAGE="te-fc:2.0.2"
NAME="te-gene-featurecounts"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

skip() { echo "SKIP [$NAME] $*"; exit 0; }
fail() { echo "FAIL [$NAME] $*" >&2; exit 1; }
pass() { echo "PASS [$NAME] $*"; }

command -v docker >/dev/null 2>&1 || skip "docker not on PATH — skipping container smoke tests"
docker image inspect "$IMAGE" >/dev/null 2>&1 || \
  skip "image '$IMAGE' not built — run env/build.sh to enable smoke tests"

# (a)+(b) version check.
ver="$(docker run --rm "$IMAGE" featureCounts -v 2>&1 | grep -oE 'v2\.0\.2' | head -1 || true)"
[[ "$ver" == "v2.0.2" ]] || fail "featureCounts version is '$ver', expected v2.0.2"
pass "image $IMAGE present; featureCounts $ver"

# (c) synthetic mini run: 3 SAF features, a 2-read SAM -> BAM, integer counts.
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# Tiny SAF: 3 features on contig chrT.
cat > "$WORK/mini.saf" <<'EOF'
GeneID	Chr	Start	End	Strand
featA	chrT	1	100	+
featB	chrT	201	300	+
featC	chrT	401	500	+
EOF

# Tiny SAM: header + 2 single-end reads, one in featA, one in featB.
cat > "$WORK/mini.sam" <<'EOF'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chrT	LN:600
r1	0	chrT	10	60	20M	*	0	0	ACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIII
r2	0	chrT	210	60	20M	*	0	0	ACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIII
EOF

# Run entirely inside the container (no samtools needed on host).
# featureCounts accepts SAM directly, so we feed mini.sam.
docker run --rm -u "$(id -u):$(id -g)" -v "$WORK":"$WORK" "$IMAGE" \
  bash -lc "cd '$WORK' && featureCounts -F SAF -a mini.saf -o out.txt mini.sam > fc.log 2>&1" \
  || { cat "$WORK/fc.log" 2>/dev/null >&2; fail "synthetic featureCounts run failed"; }

[[ -f "$WORK/out.txt" ]] || fail "no output matrix produced"

# Column 7 = counts. Assert integer, 3 feature rows, and featA/featB each got 1.
res="$(awk '
  BEGIN{rows=0; bad=0}
  /^#/ {next}
  $1=="Geneid" {next}                  # skip the column-header line
  {
    rows++
    c=$NF
    if (c !~ /^[0-9]+$/) bad++
    counts[$1]=c
  }
  END{
    printf "rows=%d bad=%d featA=%s featB=%s featC=%s\n", rows, bad, counts["featA"], counts["featB"], counts["featC"]
  }' "$WORK/out.txt")"
echo "  mini result: $res"

eval "$res"   # exports rows, bad, featA, featB, featC
[[ "${rows:-0}" -eq 3 ]]   || fail "expected 3 feature rows, got ${rows:-0}"
[[ "${bad:-1}" -eq 0 ]]    || fail "non-integer count(s) in output ($bad)"
[[ "${featA:-x}" == "1" ]] || fail "featA count expected 1, got ${featA:-<none>}"
[[ "${featB:-x}" == "1" ]] || fail "featB count expected 1, got ${featB:-<none>}"
[[ "${featC:-x}" == "0" ]] || fail "featC count expected 0, got ${featC:-<none>}"
pass "synthetic mini run: 3 integer features, counts featA=1 featB=1 featC=0"

echo "PASS [$NAME] all smoke tests passed"
exit 0
