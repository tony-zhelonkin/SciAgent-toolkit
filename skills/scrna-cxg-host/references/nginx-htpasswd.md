# nginx reverse proxy + htpasswd basic auth

The reference deployment puts nginx in front of every cellxgene container — for two reasons. First, basic-auth: cellxgene itself has no auth. Second, WebSocket pass-through: CXG uses upgrade-style WebSocket connections for diff-exp and the embedding picker, and nginx's default proxy strips them. The shipped config sets the upgrade headers correctly.

## The canonical nginx.conf (per dataset)

```nginx
worker_processes  1;

events {
  worker_connections 1024;
}

http {
  sendfile on;
  client_max_body_size 32m;

  upstream cellxgene_upstream {
    server cellxgene-<<DATASET>>:<<CXG_PORT>>;
    keepalive 64;
  }

  server {
    listen <<NGINX_PORT>>;
    server_name _;

    auth_basic            "Restricted";
    auth_basic_user_file  /etc/nginx/htpasswd;

    location / {
      proxy_set_header Host                $host;
      proxy_set_header X-Forwarded-For     $remote_addr;
      proxy_set_header X-Forwarded-Proto   $scheme;

      proxy_http_version 1.1;
      proxy_set_header   Upgrade           $http_upgrade;
      proxy_set_header   Connection        "upgrade";

      proxy_pass http://cellxgene_upstream;
    }
  }
}
```

Three load-bearing details:

- **`upstream cellxgene_upstream { server cellxgene-<dataset>:<port>; }`** — the hostname `cellxgene-<dataset>` resolves via Docker's internal DNS to the matching service. This is why the cellxgene service name in `compose.yml` must match `<<DATASET>>`.
- **`proxy_http_version 1.1` + `Upgrade`/`Connection` headers** — without these, WebSocket upgrades fail; the client falls back to long-poll, the embedding picker becomes laggy, diff-exp fails.
- **`auth_basic_user_file /etc/nginx/htpasswd`** — the file is mounted read-only from the host. Per `compose.yml`, both nginx services share the same htpasswd by default; per-dataset auth uses one file per nginx service.

## htpasswd creation

Apache `htpasswd` is the canonical tool. Bcrypt (`-B`) is preferred over the older `crypt` (default).

```bash
# First user (and create the file)
htpasswd -B -c <annotations_dir>/htpasswd <username>
# Subsequent users
htpasswd -B    <annotations_dir>/htpasswd <username2>
```

If `htpasswd` is not installed (Alpine: `apk add apache2-utils`), use docker:

```bash
docker run --rm -ti httpd:alpine htpasswd -nbB <user> <password> >> <annotations_dir>/htpasswd
```

The `-nbB` form prints the line to stdout (`-n`), takes the password from the command line (`-b`), and uses bcrypt (`-B`).

## Permission caveat

The htpasswd file must be *readable* from inside the nginx container. The container runs as a non-privileged user; if the host file is `0600`, nginx returns `Permission denied` and 401 for every request.

```bash
chmod 0644 <annotations_dir>/htpasswd
```

`0644` is fine — bcrypt-hashed passwords are not crackable in offline-feasible time at default rounds. If your threat model demands `0600`, run nginx as the file's UID (set `user:` in compose) — but that's a more complex change.

## WebSocket diagnostics

When CXG's diff-exp panel hangs:

```bash
curl -i --user "$U:$P" \
     -H "Connection: Upgrade" \
     -H "Upgrade: websocket" \
     -H "Sec-WebSocket-Version: 13" \
     -H "Sec-WebSocket-Key: $(openssl rand -base64 16)" \
     http://<HOSTNAME>:<NGINX_PORT_FULL>/api/v0.2/...
```

A correct setup returns `HTTP/1.1 101 Switching Protocols`. A misconfigured nginx returns `200 OK` (silent failure — the upgrade was stripped).

## Per-dataset auth (Pause 3 Option B)

For compartmentalised access, each nginx config points at a different htpasswd:

```nginx
auth_basic_user_file /etc/nginx/htpasswd-<<DATASET>>;
```

And `compose.yml` mounts:

```yaml
volumes:
  - ./nginx/nginx-<<DATASET>>.conf:/etc/nginx/nginx.conf:ro
  - ./nginx/htpasswd-<<DATASET>>:/etc/nginx/htpasswd-<<DATASET>>:ro
```

The skill's `render_compose.py` accepts `--htpasswd-strategy=shared|per-dataset` and renders the right paths.

## Health-check endpoint

CXG exposes `/health` (200 + `{"status":"ok"}` JSON) when ready. nginx upstream healthchecks are an enterprise feature — for the shipped Open Source nginx, just retry from the client. After `docker compose up -d`, expect 30–60s before the first request returns 200.

## Defence-in-depth

`auth_basic` over plain HTTP is acceptable on a closed intranet because the threat is "casual unauthorised browsing", not "credential theft over the wire". For deployments where the wire matters (VPN-less remote access, multi-tenant clouds), terminate TLS at a fronting reverse proxy (Caddy, Traefik) and forward to this stack on a private network. Out of scope for the skill but documented as the deliberate scope cut.
