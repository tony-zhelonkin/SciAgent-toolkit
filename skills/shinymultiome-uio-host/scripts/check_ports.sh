#!/usr/bin/env bash
# check_ports.sh — Decision Pause S3 collision probe (shinymultiome flavour)
#
# Reads a space-separated list of TCP ports as args (or defaults to 8090, the
# shinymultiome-uio-host nginx default — one above scrna-cxg-host's 8080+).
# Prints FREE/TAKEN per port. Always exits 0; the agent reads stdout to decide
# whether to fire Decision Pause S3.

set -euo pipefail

PORTS=("${@:-8090}")
echo "Probing TCP listen ports: ${PORTS[*]}"
echo "----"

probe() {
  local p=$1
  if command -v ss >/dev/null 2>&1; then
    ss -tlnp "( sport = :$p )" 2>/dev/null | tail -n +2 || true
  elif command -v lsof >/dev/null 2>&1; then
    lsof -nP -iTCP:"$p" -sTCP:LISTEN 2>/dev/null | tail -n +2 || true
  else
    echo "  (no ss or lsof available; cannot probe)" >&2
    return 0
  fi
}

any_taken=0
for p in "${PORTS[@]}"; do
  hit=$(probe "$p" || true)
  if [ -z "$hit" ]; then
    echo "  $p: FREE"
  else
    echo "  $p: TAKEN — $hit"
    any_taken=1
  fi
done

if [ "$any_taken" -eq 1 ]; then
  echo "----"
  echo "At least one port is in use. The skill should fire Decision Pause S3 (port collision)."
fi
