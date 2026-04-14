#!/usr/bin/env python3
"""
Measure context overhead from a role's skill list.

Two metrics:
  BASELINE: name + description per skill — always loaded at session start so
            the harness knows which skills exist and what they do.
  FULL:     entire SKILL.md — loaded only when a specific skill is invoked.

Token counting uses tiktoken (cl100k_base) if available; otherwise falls back
to a chars/4 approximation (Claude's tokenizer is not cl100k_base but produces
roughly similar token counts for English text).

Usage:
  python measure_skill_overhead.py [--role ROLE_NAME]

The script resolves the toolkit root by walking up from this file, so it works
regardless of the current working directory.
"""
from __future__ import annotations

import argparse
import pathlib
import re
import sys

try:
    import yaml
except ImportError:
    sys.exit("PyYAML is required. Install with: pip install pyyaml")

try:
    import tiktoken

    _ENC = tiktoken.get_encoding("cl100k_base")

    def count_tokens(text: str) -> int:
        return len(_ENC.encode(text))

    TOKENIZER = "tiktoken cl100k_base (GPT-4 tokenizer; Claude uses a different one but counts are in the same order of magnitude)"
except ImportError:
    def count_tokens(text: str) -> int:
        return max(1, len(text) // 4)

    TOKENIZER = "char/4 approximation (install tiktoken for more accurate counts)"


def find_toolkit_root(start: pathlib.Path) -> pathlib.Path:
    """Walk up from `start` until we find a directory containing roles/ and skills/."""
    for parent in [start, *start.parents]:
        if (parent / "roles").is_dir() and (parent / "skills").is_dir():
            return parent
    raise RuntimeError(f"Could not locate toolkit root from {start}")


TOOLKIT = find_toolkit_root(pathlib.Path(__file__).resolve())


def parse_role_skills(role_file: pathlib.Path) -> list[str]:
    """Extract the skill list from a role YAML (simple, no yq dependency)."""
    text = role_file.read_text()
    skills: list[str] = []
    in_skills = False
    for line in text.splitlines():
        if re.match(r"^skills:\s*(\[\])?\s*(#.*)?$", line):
            in_skills = True
            continue
        if in_skills:
            if re.match(r"^[a-z]", line):
                # Hit the next top-level YAML key
                break
            match = re.match(r"^\s*-\s*([a-z0-9-]+)", line)
            if match:
                skills.append(match.group(1))
    return skills


def load_skill(skill_name: str, skills_dir: pathlib.Path):
    """Return (full_text, frontmatter_dict, body_text) or (None, None, None)."""
    skill_md = skills_dir / skill_name / "SKILL.md"
    if not skill_md.exists():
        flat = skills_dir / f"{skill_name}.md"
        if not flat.exists():
            return None, None, None
        skill_md = flat

    content = skill_md.read_text()
    match = re.match(r"^---\n(.*?)\n---\n?", content, re.DOTALL)
    if not match:
        return content, {}, content

    try:
        frontmatter = yaml.safe_load(match.group(1)) or {}
    except yaml.YAMLError:
        frontmatter = {}

    body = content[match.end():]
    return content, frontmatter, body


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--role", default="base", help="Role name (default: base)")
    args = parser.parse_args()

    role_file = TOOLKIT / "roles" / f"{args.role}.yaml"
    skills_dir = TOOLKIT / "skills"
    if not role_file.exists():
        sys.exit(f"Role file not found: {role_file}")

    skills = parse_role_skills(role_file)
    print(f"Toolkit:    {TOOLKIT}")
    print(f"Role:       {args.role} ({len(skills)} skills)")
    print(f"Tokenizer:  {TOKENIZER}\n")

    baseline_rows: list[tuple[str, int, int]] = []
    full_rows: list[tuple[str, int, int]] = []
    skipped: list[str] = []

    for name in skills:
        content, fm, _body = load_skill(name, skills_dir)
        if content is None:
            skipped.append(name)
            continue

        desc = fm.get("description", "") if isinstance(fm, dict) else ""
        baseline_text = f"name: {name}\ndescription: {desc}"
        baseline_rows.append((name, count_tokens(baseline_text), len(desc)))
        full_rows.append((name, count_tokens(content), len(content)))

    n = max(1, len(baseline_rows))

    baseline_total = sum(row[1] for row in baseline_rows)
    print("=== BASELINE OVERHEAD (always loaded — name + description) ===")
    print(f"  Total:             {baseline_total:>8,} tokens across {len(baseline_rows)} skills")
    print(f"  Mean per skill:    {baseline_total / n:>8.0f} tokens")
    print(f"  % of 200K window:  {100 * baseline_total / 200_000:>8.2f}%")
    print(f"  % of 1M window:    {100 * baseline_total / 1_000_000:>8.3f}%")

    print("\n  Heaviest descriptions:")
    for name, tok, chars in sorted(baseline_rows, key=lambda r: r[1], reverse=True)[:5]:
        print(f"    {tok:>5} tok  ({chars:>4} chars)  {name}")
    print("\n  Lightest descriptions (may under-trigger):")
    for name, tok, chars in sorted(baseline_rows, key=lambda r: r[1])[:5]:
        print(f"    {tok:>5} tok  ({chars:>4} chars)  {name}")

    full_total = sum(row[1] for row in full_rows)
    print("\n=== FULL-LOAD UPPER BOUND (every SKILL.md loaded — not typical) ===")
    print(f"  Total:             {full_total:>8,} tokens across {len(full_rows)} skills")
    print(f"  Mean per skill:    {full_total / n:>8.0f} tokens")
    print(f"  % of 200K window:  {100 * full_total / 200_000:>8.1f}%")
    print(f"  % of 1M window:    {100 * full_total / 1_000_000:>8.2f}%")

    print("\n  Heaviest SKILL.md files:")
    for name, tok, chars in sorted(full_rows, key=lambda r: r[1], reverse=True)[:5]:
        print(f"    {tok:>5} tok  ({chars:>6} chars)  {name}")

    if skipped:
        print(f"\n  Skipped (skill not found in library): {skipped}")


if __name__ == "__main__":
    main()
