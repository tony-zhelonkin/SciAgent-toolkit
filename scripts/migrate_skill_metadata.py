#!/usr/bin/env python3
"""
scripts/migrate_skill_metadata.py — ADR-001/003 one-shot metadata codemod.

Applies the PR 1 migration to every skills/*/SKILL.md in the toolkit:

  1. Parse the YAML frontmatter block (between the two --- delimiters).
  2. Map legacy scope values:
       foundation  →  concept
       orchestrator →  concept
       atomic       →  implementation
     Skills missing the scope field entirely receive: implementation.
  3. Insert any of the five required metadata fields that are absent,
     with an empty-list default:
       scope, requires, complementary-skills, contraindications, tags
  4. Re-serialise the frontmatter with PyYAML; write the file with the
     body (everything after the closing ---) byte-for-byte unchanged.
  5. Emit one "AUDIT-ME" log line per skill where YAML comments were
     stripped by PyYAML so a human reviewer can restore them if needed.

PyYAML 6.0.3 is used; ruamel.yaml is NOT required and must NOT be installed.

Usage (from toolkit root):
    python3 scripts/migrate_skill_metadata.py [--dry-run]

Options:
    --dry-run   Print what would change without writing files.

Exit code: 0 on success; non-zero if any file fails to parse.
"""

import argparse
import os
import re
import sys

import yaml

SKILLS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "skills")

SCOPE_MAP = {
    "foundation": "concept",
    "orchestrator": "concept",
    "atomic": "implementation",
}

REQUIRED_METADATA_FIELDS = [
    "scope",
    "requires",
    "complementary-skills",
    "contraindications",
    "tags",
]

# Detect YAML comment lines inside the frontmatter block.
_COMMENT_RE = re.compile(r"^\s*#", re.MULTILINE)


def split_frontmatter(text):
    """Split a SKILL.md into (frontmatter_text, body_text).

    frontmatter_text is the raw YAML between the two --- delimiters
    (NOT including the delimiters themselves).  body_text is everything
    after the closing --- delimiter (including its trailing newline).

    Returns (None, text) when no valid frontmatter block is found.
    """
    if not text.startswith("---"):
        return None, text

    # Find the closing --- on its own line.
    match = re.search(r"\n---[ \t]*\n", text)
    if not match:
        return None, text

    fm_text = text[len("---\n") : match.start() + 1]  # content between the fences
    body_text = text[match.end() :]
    return fm_text, body_text


def migrate_frontmatter(fm_text, skill_name):
    """Parse, transform, and re-serialise the frontmatter YAML.

    Returns (new_fm_text, audit_comments) where audit_comments is a list
    of human-readable strings describing comment-stripping events.
    """
    audit_comments = []

    # Detect whether the raw frontmatter contained any comments.
    if _COMMENT_RE.search(fm_text):
        audit_comments.append(
            "AUDIT-ME [%s]: frontmatter contained YAML comments — stripped by PyYAML; "
            "review the diff and restore any prose that was in those comments." % skill_name
        )

    data = yaml.safe_load(fm_text)
    if data is None:
        data = {}

    # Ensure metadata block exists.
    if "metadata" not in data:
        data["metadata"] = {}

    meta = data["metadata"]

    # --- Step 2: scope mapping ---
    current_scope = meta.get("scope", None)
    if current_scope is None:
        meta["scope"] = "implementation"
    elif current_scope in SCOPE_MAP:
        meta["scope"] = SCOPE_MAP[current_scope]
    # else: already a valid new-vocabulary value — leave as-is.

    # --- Step 3: insert missing fields with empty-list default; reset tags ---
    for field in REQUIRED_METADATA_FIELDS:
        if field == "scope":
            continue  # handled above
        if field not in meta:
            meta[field] = []

    # Reset tags to empty list unconditionally.  Pre-migration tags were
    # free-form tool-level identifiers not governed by tags.yaml.  The
    # vocabulary test (test_tags_vocabulary.sh) now enforces that every tag
    # entry exists in tags.yaml; unrecognised tags are a FAIL.  Skills carry
    # tags: [] after migration; curators re-add vocabulary-compliant tags via
    # subsequent PRs as the vocabulary grows.
    meta["tags"] = []

    # Re-serialise.  PyYAML sorts keys by default when using safe_dump;
    # use sort_keys=False to preserve the insertion order of existing keys
    # (Python 3.7+ dicts are ordered).
    new_fm_text = yaml.dump(
        data,
        allow_unicode=True,
        default_flow_style=False,
        sort_keys=False,
        width=10000,
    )
    return new_fm_text, audit_comments


def process_skill(skill_dir, dry_run=False):
    """Migrate one SKILL.md; return list of audit strings."""
    skill_name = os.path.basename(skill_dir)
    filepath = os.path.join(skill_dir, "SKILL.md")

    with open(filepath, "r", encoding="utf-8") as fh:
        original = fh.read()

    fm_text, body_text = split_frontmatter(original)
    if fm_text is None:
        return ["AUDIT-ME [%s]: no frontmatter found — skipped." % skill_name]

    new_fm_text, audit = migrate_frontmatter(fm_text, skill_name)

    new_content = "---\n" + new_fm_text + "---\n" + body_text

    if new_content == original:
        return audit  # no change needed

    if not dry_run:
        with open(filepath, "w", encoding="utf-8") as fh:
            fh.write(new_content)
        print("migrated: %s" % skill_name)
    else:
        print("dry-run:  %s  (would change)" % skill_name)

    return audit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true", help="Print changes without writing files.")
    args = parser.parse_args()

    skill_dirs = sorted(
        [
            os.path.join(SKILLS_DIR, d)
            for d in os.listdir(SKILLS_DIR)
            if not d.startswith((".", "_"))
            and os.path.isdir(os.path.join(SKILLS_DIR, d))
        ]
    )

    all_audits = []
    errors = 0

    for skill_dir in skill_dirs:
        try:
            audits = process_skill(skill_dir, dry_run=args.dry_run)
            all_audits.extend(audits)
        except Exception as exc:
            print("ERROR [%s]: %s" % (os.path.basename(skill_dir), exc), file=sys.stderr)
            errors += 1

    if all_audits:
        print("")
        print("--- AUDIT LOG ---")
        for line in all_audits:
            print(line)

    if errors:
        print("\n%d skill(s) failed to process." % errors, file=sys.stderr)
        sys.exit(1)

    print("\nDone. %d skill(s) processed." % len(skill_dirs))


if __name__ == "__main__":
    main()
