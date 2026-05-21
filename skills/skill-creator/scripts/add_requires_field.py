#!/usr/bin/env python3
"""
Migration helper (M1, per ADR-0002 §5): ensure every SKILL.md under
SciAgent-toolkit/skills/ carries `metadata.requires: []` and
`metadata.scope: atomic`.

Rules:
  * If the SKILL.md has no YAML frontmatter, the file is skipped with a
    warning (these are usually fixtures or templates).
  * If the frontmatter has no `metadata:` block, one is inserted right
    after `license:` (or after `description:` if `license:` is absent).
  * If `metadata:` exists but is missing `requires:`, the key is appended
    at the top of the metadata block.
  * Likewise for `scope:` — default value is `atomic`.
  * Existing keys/values are left untouched. Formatting / comments are
    preserved by operating on the raw text rather than round-tripping
    through PyYAML.

Run from the repository root:

    python3 skills/skill-creator/scripts/add_requires_field.py

The script is idempotent — re-running it on an already-migrated tree
makes zero changes.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
SKILLS_DIR = REPO_ROOT / "skills"


def split_frontmatter(text: str) -> tuple[str, str, str] | None:
    """Return (head_marker, frontmatter_body, rest_of_file) or None."""
    if not text.startswith("---\n"):
        return None
    end = text.find("\n---\n", 4)
    if end == -1:
        # Allow trailing form `---` at EOF or `---\n` w/o trailing newline.
        end = text.find("\n---", 4)
        if end == -1:
            return None
    head = text[:4]                      # "---\n"
    fm = text[4:end] + "\n"              # body, normalised to end with \n
    tail_start = end + len("\n---\n")
    if tail_start > len(text):
        tail_start = len(text)
    rest = "---\n" + text[tail_start:]
    return head, fm, rest


def fm_has_top_key(fm: str, key: str) -> bool:
    return re.search(rf"(?m)^{re.escape(key)}:", fm) is not None


def fm_metadata_block(fm: str) -> tuple[int, int] | None:
    """Return (start_line_idx, end_line_idx) of the metadata: block."""
    lines = fm.splitlines(keepends=True)
    start = None
    for i, ln in enumerate(lines):
        if re.match(r"^metadata:\s*(#.*)?$", ln):
            start = i
            break
    if start is None:
        return None
    # End = next top-level key (no leading whitespace, ends with colon) or EOF.
    end = len(lines)
    for j in range(start + 1, len(lines)):
        if re.match(r"^[A-Za-z_][A-Za-z0-9_-]*:", lines[j]):
            end = j
            break
    return start, end


def metadata_has_nested_key(fm: str, key: str) -> bool:
    block = fm_metadata_block(fm)
    if block is None:
        return False
    s, e = block
    lines = fm.splitlines(keepends=True)
    pat = re.compile(rf"^[ \t]+{re.escape(key)}:")
    for ln in lines[s + 1 : e]:
        if pat.match(ln):
            return True
    return False


def insert_top_key(fm: str, key: str, value: str) -> str:
    """Insert `key: value\n` either after `license:` or after `description:`."""
    lines = fm.splitlines(keepends=True)
    anchor = None
    for i, ln in enumerate(lines):
        if re.match(r"^license:", ln):
            anchor = i
            break
    if anchor is None:
        for i, ln in enumerate(lines):
            if re.match(r"^description:", ln):
                anchor = i
                break
    if anchor is None:
        anchor = len(lines) - 1
    lines.insert(anchor + 1, f"{key}:\n{value}")
    return "".join(lines)


def add_nested_key(fm: str, key: str, value: str) -> str:
    """Insert `  key: value` right after the `metadata:` line."""
    block = fm_metadata_block(fm)
    if block is None:
        return fm
    s, _ = block
    lines = fm.splitlines(keepends=True)
    lines.insert(s + 1, f"  {key}: {value}\n")
    return "".join(lines)


def migrate_one(path: Path) -> str:
    """Return 'changed', 'skipped:<reason>', or 'noop'."""
    text = path.read_text()
    split = split_frontmatter(text)
    if split is None:
        return "skipped:no-frontmatter"

    head, fm, rest = split
    original_fm = fm

    # Ensure metadata: block exists.
    if not fm_has_top_key(fm, "metadata"):
        fm = insert_top_key(fm, "metadata", "")  # blank value, body added below

    # Now ensure metadata.requires and metadata.scope exist.
    if not metadata_has_nested_key(fm, "requires"):
        fm = add_nested_key(fm, "requires", "[]")
    if not metadata_has_nested_key(fm, "scope"):
        fm = add_nested_key(fm, "scope", "atomic")

    if fm == original_fm:
        return "noop"

    path.write_text(head + fm + rest if rest != "---\n" else head + fm + "---\n")
    return "changed"


def main() -> int:
    if not SKILLS_DIR.is_dir():
        print(f"skills dir not found: {SKILLS_DIR}", file=sys.stderr)
        return 1
    changed = 0
    skipped = 0
    noop = 0
    for skill_md in sorted(SKILLS_DIR.glob("*/SKILL.md")):
        # Skip the fixture template — it has placeholder text that must
        # not be modified mechanically.
        if skill_md.parent.name == "_TEMPLATE":
            continue
        result = migrate_one(skill_md)
        if result == "changed":
            print(f"changed: {skill_md.relative_to(REPO_ROOT)}")
            changed += 1
        elif result.startswith("skipped"):
            print(f"skipped ({result.split(':',1)[1]}): {skill_md.relative_to(REPO_ROOT)}")
            skipped += 1
        else:
            noop += 1
    print(f"\nsummary: changed={changed} noop={noop} skipped={skipped}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
