---
name: container-port-tunnel
description: >-
  Expose a browser tool (marimo, httpgd, dashboard) running inside the
  devcontainer to the user's laptop. Use whenever the user needs a URL for a
  server started in this container, or reports that a served URL fails.
license: MIT
---

# Container port tunnel

## Topology
A project's sessions share one devcontainer, and every server they start listens
inside it. The compose template ships its `ports:` block commented out, so the
container publishes nothing and its bridge IP is reachable only from the lab
host. Access path:
laptop -> ssh -> lab host -> docker bridge -> container IP -> port.

## Scan
```bash
hostname -I                                    # container IP (changes on rebuild)
ps aux | grep -Ei "marimo|httpgd|serve" | grep -v grep   # running servers, ports, notebooks
```
Listener check — `ss`, `netstat` and `lsof` are all absent from the image, so
probe the ports directly:
```bash
python3 -c "import socket; [print(p, socket.socket().connect_ex(('127.0.0.1', p))==0) for p in (2718,2719,2782,8787)]"
```
Tokens live in the launcher log: `grep -o "access_token=[A-Za-z0-9_-]*" <log>`.

## Serve
```bash
marimo edit --headless --host 0.0.0.0 --port <PORT> <notebook.py> > <log> 2>&1 &
```
Pick a free port. Keep the default access token.

## Hand the user exactly two lines, in this order
1. Tunnel, run on the laptop; the target is the container IP:
```
ssh -L <PORT>:<CONTAINER_IP>:<PORT> <user>@<lab-host>
```
2. Browser:
```
http://localhost:<PORT>/?access_token=<TOKEN>
```

## Facts
- A tunnel of the form `-L PORT:localhost:PORT` lands on the host loopback,
  which has zero listeners for these servers.
- The container IP is assigned at container creation; rescan `hostname -I`
  after every rebuild and reissue the tunnel command.
- Durable option: publish `127.0.0.1:<range>:<range>` under `ports:` in
  `.devcontainer/docker-compose.yml`; it applies at the next rebuild and makes
  `-L PORT:localhost:PORT` tunnels work permanently.
