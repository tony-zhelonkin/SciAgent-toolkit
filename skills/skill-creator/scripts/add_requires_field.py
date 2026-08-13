#!/usr/bin/env python3
"""
==============================================================================
DO NOT RUN THIS SCRIPT. HISTORICAL / INERT — kept as a record only.
==============================================================================

The `metadata.requires` / `metadata.scope` frontmatter fields this script
injects were part of the role/taxonomy layer that has since been demolished
(see `docs/architecture.md` and `docs/skills.md` "taxonomy removed" note).
`metadata:` is no longer populated on any SKILL.md, there is no resolver
that reads `requires:` or `scope:`, and `scio link` mounts the whole
catalog unconditionally regardless of these keys.

Running this script today would resurrect exactly the frontmatter block a
later refactor deliberately deleted across all 83 skills. It is retained
in the repo purely as a historical record of the M1 migration step, not as
a tool to invoke. `main()` hard-exits before touching any file — see below.

------------------------------------------------------------------------------
Original docstring (historical, describes what the script DID when it was
a live migration step; no longer applicable):

Migration helper (M1, per ADR-0002 §5): ensure every SKILL.md under
The Scio toolkit's skills carry `metadata.requires: []` and
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

The script was idempotent — re-running it on an already-migrated tree
made zero changes. That property is irrelevant now: it must not be run
at all.
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
    print(
        "This script is INERT and must not be run: it injects "
        "`metadata.requires` / `metadata.scope` frontmatter for a "
        "role/taxonomy layer that has been removed from the Scio toolkit. "
        "Running it would resurrect frontmatter a later refactor "
        "deliberately deleted from every SKILL.md. See the module "
        "docstring at the top of this file for details. Exiting without "
        "changing any file.",
        file=sys.stderr,
    )
    return 1


def _disabled_main() -> int:
    """Historical body of main(), preserved for reference. Not called."""
    if not SKILLS_DIR.is_dir():
        print(f"skills dir not found: {SKILLS_DIR}", file=sys.stderr)
        return 1
    changed = 0
    skipped = 0
    noop = 0
    for skill_md in sorted(SKILLS_DIR.glob("*/SKILL.md")):
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
