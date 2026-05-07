# Port-collision detection and remap

Two failure modes the deploy must catch *before* `docker compose up` rather than after. First, a default cellxgene/nginx port is taken — the deploy crashes with `port is already allocated`. Second, the deploy *succeeds* but the wet-lab URL points at someone else's service that happened to be on the same port — the worst kind of silent failure.

## Defaults

```
cellxgene internal:  5005, 5006, 5007, ...   (one per dataset)
nginx external:      8080, 8081, 8082, ...   (one per dataset, browser-facing)
```

The cellxgene ports are container-network-only — collisions there only matter inside the docker network, which means they never collide with host services. The nginx ports are host-bound and *do* collide with anything else on the host.

## Detection — `ss -tlnp`

```bash
#!/usr/bin/env bash
# scripts/check_ports.sh — Decision Pause 2 collision probe
set -euo pipefail

PORTS=("${@:-8080 8081 8082}")
echo "Probing TCP listen ports: ${PORTS[*]}"
echo "----"
for p in "${PORTS[@]}"; do
  hit=$(ss -tlnp "( sport = :$p )" 2>/dev/null | tail -n +2 || true)
  if [ -z "$hit" ]; then
    echo "  $p: FREE"
  else
    echo "  $p: TAKEN — $hit"
  fi
done
```

The script returns 0 always; the agent reads stdout and decides whether to fire Decision Pause 2 (collision branch).

`ss` is preferred over `netstat` (`netstat` is sometimes missing in minimal Linux installs; `ss` is in `iproute2`, almost universally present). The `( sport = :$p )` filter scopes to the listen-port TCP socket; `tail -n +2` strips the column header.

## Auto-remap heuristic

When at least one default is taken, the skill proposes a remap:

```python
def find_free_port(start: int) -> int:
    p = start
    while is_port_taken(p):
        p += 1
    return p

# Per-dataset remap, contiguous block above default
new_ports = {}
nginx_p = find_free_port(8080)
for ds in datasets:
    new_ports[f"NGINX_PORT_{ds.upper()}"] = nginx_p
    nginx_p = find_free_port(nginx_p + 1)
```

Cellxgene-internal ports rarely need remapping (no host-side conflict); the skill keeps `5005+` unless the user explicitly overrides.

## What the user sees in Pause 2

```
Port collision check (host:port):
  8080: TAKEN — users:(("apache2",pid=1234,fd=4))
  8081: FREE
  8082: FREE

Decision Pause — Port collision

Question: Default port 8080 is taken by apache2. How should I remap?

Option A (Use defaults)        — fails immediately on `docker compose up`
Option B (Auto-remap)          — proposed: 8083 → full, 8081 → t_cells, 8082 → b_cells
Option C (User-specified)      — name ports explicitly
```

The user picks; the skill writes back to `.env` and `analysis_config.yaml::decisions::scrna-cxg-host::ports`.

## Strict-policy servers

Some institutional servers restrict outbound listening to a specific range (e.g., `8000–8999` or even `30000+`). Pause 2 Option C handles this — the user names the range, the skill picks free ports inside it. The skill records the policy hint in config:

```yaml
decisions:
  scrna-cxg-host:
    ports:
      policy: "host_range_8000_8999"   # advisory; agent reads on re-deploy
      cxg:
        full: 5005
      nginx:
        full: 8087
        t_cells: 8088
```

Re-deploys re-probe and warn if the recorded ports are now taken (e.g., another lab's stack started in the meantime).

## Diagnosing a *silent* port collision

The dangerous case: `docker compose up` succeeds because no one else's container holds 8080, but the host's apache2 (or another reverse proxy) is binding on the same interface at a slightly higher precedence. Symptom: the wet lab opens the URL and sees apache2's default page, not CXG.

Check after `docker compose up`:

```bash
curl -s -o /dev/null -w "%{http_code}\n" --user "$U:$P" http://<HOSTNAME>:8080/
# expect 200
curl -s -I http://<HOSTNAME>:8080/ | grep -i server
# expect 'Server: nginx', not 'Server: Apache'
```

If apache2 (or anything not nginx) responds, the collision is at the host firewall / iptables / interface-binding layer — outside docker's port allocation, but resolvable by remap.
