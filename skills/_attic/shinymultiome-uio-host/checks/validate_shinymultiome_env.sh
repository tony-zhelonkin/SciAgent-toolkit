#!/usr/bin/env bash
# validate_shinymultiome_env.sh — Phase B post-up smoke test
#
# Runs after `docker compose up -d`. Confirms the container is healthy at
# multiple layers: image build, R package availability, R startup of the app,
# nginx reverse proxy, websocket upgrade, and basic-auth.
#
# Usage:
#     bash checks/validate_shinymultiome_env.sh [hostname] [port] [user] [password]
#
# Environment fall-back: if positional args are missing, reads from .env in CWD.
# The function returns non-zero on any failure and writes a checklist to stdout.

set -uo pipefail

HOSTNAME_ARG="${1:-}"
PORT_ARG="${2:-}"
USER_ARG="${3:-}"
PASS_ARG="${4:-}"

# Read .env if present (no-op if vars already set)
if [ -f ".env" ]; then
  # shellcheck disable=SC1091
  set -a; . ./.env; set +a
fi

HOST="${HOSTNAME_ARG:-${HOSTNAME:-localhost}}"
PORT="${PORT_ARG:-${NGINX_PORT_SHINYMULTIOME:-8090}}"
U="${USER_ARG:-}"
P="${PASS_ARG:-}"

URL="http://${HOST}:${PORT}/"

ok=0
fail=0
report() {
  local tag=$1; shift
  if [ "$tag" = "ok" ]; then
    echo "  [OK  ] $*"; ok=$((ok+1))
  else
    echo "  [FAIL] $*"; fail=$((fail+1))
  fi
}

echo "Validating shinymultiome stack at ${URL}"
echo "----"

# 1. Docker compose service up
if command -v docker >/dev/null 2>&1; then
  state=$(docker compose ps --format json 2>/dev/null \
            | python3 -c "import sys,json; xs=[json.loads(l) for l in sys.stdin if l.strip()]; \
                          [print(x.get('Service'), x.get('State')) for x in xs]" 2>/dev/null \
            | grep -E '^(shinymultiome|nginx-shinymultiome) ' || true)
  if echo "$state" | grep -q "shinymultiome running" && \
     echo "$state" | grep -q "nginx-shinymultiome running"; then
    report ok "docker compose: both services running"
  else
    report fail "docker compose state:\n$state"
  fi
else
  report fail "docker CLI not found"
fi

# 2. R package availability inside the container
if docker compose exec -T shinymultiome \
   R --no-save -e "stopifnot(requireNamespace('Signac', quietly = TRUE)); \
                    stopifnot(requireNamespace('Seurat', quietly = TRUE)); \
                    stopifnot(requireNamespace(Sys.getenv('BSGENOME_PKG'), quietly = TRUE))" \
   >/dev/null 2>&1; then
  report ok "R packages available in container (Seurat, Signac, BSgenome)"
else
  report fail "R package check failed inside container"
fi

# 3. RDS path resolves inside container
if docker compose exec -T shinymultiome \
   sh -c 'test -f "$RDS_PATH"' >/dev/null 2>&1; then
  report ok "RDS_PATH resolves inside container"
else
  report fail "RDS_PATH not accessible inside container"
fi

# 4. Fragment dir resolves inside container
if docker compose exec -T shinymultiome \
   sh -c 'test -d "$FRAGMENTS_DIR" && ls "$FRAGMENTS_DIR"/*.tbi >/dev/null 2>&1' \
   >/dev/null 2>&1; then
  report ok "FRAGMENTS_DIR populated with .tbi files inside container"
else
  report fail "FRAGMENTS_DIR check failed (no .tbi files visible)"
fi

# 5. nginx 401 without auth, 200 with auth
code_no_auth=$(curl -sS -o /dev/null -w "%{http_code}" "$URL" 2>/dev/null || echo 000)
if [ "$code_no_auth" = "401" ]; then
  report ok "nginx returns 401 without basic-auth"
else
  report fail "nginx returned $code_no_auth without auth (expected 401)"
fi

if [ -n "$U" ] && [ -n "$P" ]; then
  code_auth=$(curl -sS -o /dev/null -w "%{http_code}" --user "${U}:${P}" "$URL" 2>/dev/null || echo 000)
  if [ "$code_auth" = "200" ]; then
    report ok "nginx returns 200 with basic-auth"
  else
    report fail "nginx returned $code_auth with auth (expected 200)"
  fi

  # 6. WebSocket upgrade
  ws_code=$(curl -sS -o /dev/null -w "%{http_code}" \
    -H "Connection: Upgrade" -H "Upgrade: websocket" \
    -H "Sec-WebSocket-Version: 13" \
    -H "Sec-WebSocket-Key: $(openssl rand -base64 16 2>/dev/null || echo testkey)" \
    --user "${U}:${P}" \
    "${URL}websocket" 2>/dev/null || echo 000)
  if [ "$ws_code" = "101" ]; then
    report ok "websocket upgrade returns 101 Switching Protocols"
  else
    report fail "websocket upgrade returned $ws_code (expected 101 — Shiny will disconnect)"
  fi
else
  echo "  [skip] basic-auth + websocket checks (no user/pass provided)"
fi

echo "----"
echo "  ${ok} OK / ${fail} FAIL"
[ "$fail" -eq 0 ]
