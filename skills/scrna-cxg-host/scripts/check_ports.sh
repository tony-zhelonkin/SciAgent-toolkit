#!/usr/bin/env bash
# check_ports.sh — Decision Pause 2 collision probe.
#
# Reads a space-separated list of TCP ports as args (or defaults to the
# scrna-cxg-host nginx defaults 8080..8082) and prints FREE/TAKEN per port.
# Always exits 0; the agent reads stdout to decide whether to fire the
# port-collision Decision Pause.
#
# Uses `ss` (iproute2; usually already installed). Falls back to `lsof` if
# `ss` is missing.

set -euo pipefail

PORTS=("${@:-8080 8081 8082}")
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
  echo "At least one port is in use. The skill should fire Decision Pause 2 (port collision)."
fi
