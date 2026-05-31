#!/usr/bin/env bash
# golden_render.sh — Render the golden fixture and assert correctness.
#
# Assertions:
#   1. validate_components.py exits 0 on the golden fixture.
#   2. render_treemap.py exits 0.
#   3. Output HTML exists.
#   4. Self-contained: no http-scheme <script src=...> or <link href=...> pointing externally.
#   5. Data injected: output contains the fixture's project name.
#   6. Byte-stable across two runs (audit_date slot normalised before diff).
#
# Usage (from the skill directory):
#   bash checks/golden_render.sh
# Or from anywhere:
#   bash skills/architecture-treemap/checks/golden_render.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SKILL_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

FIXTURE="$SKILL_DIR/references/example/components.json"
VALIDATE="$SKILL_DIR/scripts/validate_components.py"
RENDER="$SKILL_DIR/scripts/render_treemap.py"
OUT1="/tmp/golden_render_1.html"
OUT2="/tmp/golden_render_2.html"

PASS=0
FAIL=0

ok()   { echo "  PASS: $1"; PASS=$((PASS + 1)); }
fail() { echo "  FAIL: $1"; FAIL=$((FAIL + 1)); }

echo "=== golden_render.sh ==="
echo "Skill dir: $SKILL_DIR"
echo "Fixture:   $FIXTURE"
echo ""

# Guard: fixture must exist
if [ ! -f "$FIXTURE" ]; then
  echo "ERROR: fixture not found at $FIXTURE"
  exit 1
fi

# ── 1. Validate fixture ───────────────────────────────────────────────────────
echo "[1] validate_components.py exits 0 on golden fixture"
if python3 "$VALIDATE" "$FIXTURE" > /dev/null 2>&1; then
  ok "validator exit 0"
else
  fail "validator returned non-zero"
  python3 "$VALIDATE" "$FIXTURE" || true
fi

# ── 2. Render (first run) ─────────────────────────────────────────────────────
echo "[2] render_treemap.py exits 0"
if python3 "$RENDER" "$FIXTURE" --output "$OUT1" > /dev/null 2>&1; then
  ok "renderer exit 0"
else
  fail "renderer returned non-zero"
  python3 "$RENDER" "$FIXTURE" --output "$OUT1" || true
fi

# ── 3. Output file exists ─────────────────────────────────────────────────────
echo "[3] Output HTML exists"
if [ -f "$OUT1" ]; then
  ok "output file present"
else
  fail "output file missing: $OUT1"
fi

# ── 4. Self-contained: no external script/link src ───────────────────────────
echo "[4] Self-contained (no http-scheme script src or link href)"
if [ -f "$OUT1" ]; then
  # Match src="http or href="http inside tag attribute context.
  # Use grep -P if available (Linux), else extended regex.
  if grep -qP '(?:src|href)\s*=\s*"https?://' "$OUT1" 2>/dev/null; then
    fail "found external http-scheme src/href attribute"
  elif grep -qE '(src|href)=["\047]https?://' "$OUT1" 2>/dev/null; then
    fail "found external http-scheme src/href attribute"
  else
    ok "no external script/link references in HTML attributes"
  fi
else
  fail "output file missing; cannot check self-containment"
fi

# ── 5. Data injected: known fixture value present ─────────────────────────────
echo "[5] Injected data contains fixture project name"
if [ -f "$OUT1" ] && grep -q 'example-service' "$OUT1"; then
  ok "fixture project name found in output"
else
  fail "fixture project name not found in output"
fi

# ── 6. Byte-stable across two runs ───────────────────────────────────────────
echo "[6] Byte-stable across two renders"
python3 "$RENDER" "$FIXTURE" --output "$OUT2" > /dev/null 2>&1

if [ -f "$OUT1" ] && [ -f "$OUT2" ]; then
  # Normalise audit_date value before diffing to isolate any date-derived slot.
  # In practice the renderer is fully deterministic (no Date.now() calls), but
  # this guard makes the check resilient to future additions of a render timestamp.
  normalize() {
    sed 's/"audit_date":"[^"]*"/"audit_date":"NORMALIZED"/g' "$1"
  }

  DIFF_OUT=$(diff <(normalize "$OUT1") <(normalize "$OUT2") || true)
  if [ -z "$DIFF_OUT" ]; then
    ok "outputs are byte-identical (after normalization)"
  else
    fail "outputs differ across two runs"
    echo "--- diff (first 20 lines) ---"
    echo "$DIFF_OUT" | head -20
  fi
else
  fail "one or both output files missing; cannot diff"
fi

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "=== Result: $PASS passed, $FAIL failed ==="

if [ "$FAIL" -gt 0 ]; then
  exit 1
fi
echo "golden_render.sh: ALL ASSERTIONS PASSED"
exit 0
