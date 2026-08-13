#!/usr/bin/env python3
"""render_compose.py — substitute .env values into shinymultiome templates.

Reads `.env` from the current directory (or --env-file), then:

  1. docker-compose.shinymultiome.yml.template
       -> docker-compose.shinymultiome.yml
  2. nginx-shinymultiome.conf.template
       -> nginx-shinymultiome.conf

Substitutions: every `<<KEY>>` in the templates is replaced with `.env::KEY`.
Unsubstituted `<<...>>` after the pass triggers a non-zero exit so the user
notices missing env vars before `docker compose up`.

Usage:
    cp scripts/env.template .env
    $EDITOR .env
    python scripts/render_compose.py [--env-file .env] [--out-dir .]
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

PLACEHOLDER_RE = re.compile(r"<<([A-Z_][A-Z0-9_]*)>>")

DEFAULT_TEMPLATES = [
    ("scripts/docker-compose.shinymultiome.yml.template", "docker-compose.shinymultiome.yml"),
    ("scripts/nginx-shinymultiome.conf.template",         "nginx-shinymultiome.conf"),
]


def parse_env(path: Path) -> dict[str, str]:
    env: dict[str, str] = {}
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        if "=" not in s:
            continue
        k, v = s.split("=", 1)
        env[k.strip()] = v.strip().strip('"').strip("'")
    return env


def render(text: str, env: dict[str, str]) -> tuple[str, set[str]]:
    missing: set[str] = set()

    def repl(m: re.Match[str]) -> str:
        key = m.group(1)
        if key in env:
            return env[key]
        missing.add(key)
        return m.group(0)

    return PLACEHOLDER_RE.sub(repl, text), missing


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    p.add_argument("--env-file", default=".env")
    p.add_argument("--out-dir", default=".")
    p.add_argument("--templates", nargs="*", default=None,
                   help="Override default template list as src:dst pairs")
    args = p.parse_args()

    env_path = Path(args.env_file)
    if not env_path.exists():
        print(f"ERROR: {env_path} not found. Copy scripts/env.template first.", file=sys.stderr)
        return 1
    env = parse_env(env_path)

    if args.templates:
        templates = []
        for t in args.templates:
            src, dst = t.split(":", 1)
            templates.append((src, dst))
    else:
        templates = DEFAULT_TEMPLATES

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    overall_missing: set[str] = set()
    for src, dst in templates:
        src_path = Path(src)
        if not src_path.exists():
            # try resolving relative to this script's grandparent (project root)
            here = Path(__file__).resolve().parent
            alt = here.parent / src
            if alt.exists():
                src_path = alt
            else:
                print(f"ERROR: template not found: {src}", file=sys.stderr)
                return 1
        rendered, missing = render(src_path.read_text(), env)
        out_path = out_dir / dst
        out_path.write_text(rendered)
        print(f"  rendered: {src_path.name} -> {out_path}")
        overall_missing |= missing

    if overall_missing:
        print(f"\nERROR: unsubstituted placeholders: {sorted(overall_missing)}", file=sys.stderr)
        print("Add these keys to .env and re-run.", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
