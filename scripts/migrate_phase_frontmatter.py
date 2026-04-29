#!/usr/bin/env python3
"""
One-time migration: canonicalise `phase-NN.md` frontmatter across a portfolio.

Brings every `docs/{feature}/plan/phase-NN.md` into the canonical schema set by
`commands/implement.md` Phase 3 step 5, so `/status` and `/verify` can lean on
frontmatter instead of body-prose heuristics.

Canonical schema (target):

    ---
    phase: <N>
    feature: <slug>
    status: shipped | in-progress | not-started | deferred | blocked
    completed: YYYY-MM-DD          # may be empty
    commit: <SHA or empty>
    files_touched:
      - <path1>
    verification: pass | partial | deferred
    notes: <optional one-liner>
    ---

The script INFERS values rather than guessing:

1. If existing frontmatter has `status: DONE` → migrate to `status: shipped`.
2. Body markers `**Status: ✅** ... shipped <date>` or `**P0 ... landed <date>` →
   set status=shipped and parse the date into `completed:`.
3. `## Files to Create` and `## Files to Modify` body sections → seed
   `files_touched:` (user can refine after).
4. `commit:` is best-effort: `git log --follow -1 --pretty=%H -- <first file>`
   bounded by the completed date if present.
5. If nothing infers → status: not-started, leave fields empty.

DESTRUCTIVE: rewrites every phase-NN.md frontmatter block in place. Run on a
clean working tree so a `git diff` shows exactly what changed. The body of the
phase doc (everything after the `---` block) is never touched.

Usage:
    python migrate_phase_frontmatter.py --docs-dir /path/to/docs [--dry-run]
    python migrate_phase_frontmatter.py --docs-dir /path/to/docs --feature umap

Examples:
    cd /workspaces/DC_hum_verse/01_modules/pathway-explorer
    python ../SciAgent-toolkit/scripts/migrate_phase_frontmatter.py \\
        --docs-dir docs --dry-run
"""

import argparse
import re
import subprocess
from datetime import date
from pathlib import Path


FRONTMATTER_RE = re.compile(r"^---\s*\n(.*?\n)---\s*\n", re.DOTALL)
DATE_RE = re.compile(r"\b(20\d{2}-\d{2}-\d{2})\b")
SHIPPED_BODY_PATTERNS = [
    re.compile(
        r"\*\*(?:Status|P\d+|Subphase\s+\w+)[^*]*\*\*[^\n]*?"
        r"(?:shipped|landed|complete|done)[^\n]*?(20\d{2}-\d{2}-\d{2})",
        re.IGNORECASE,
    ),
    re.compile(
        r"\*\*\s*Subphase\s+\w+\s+(?:shipped|landed)\s*\(\s*(20\d{2}-\d{2}-\d{2})\s*\)[^*]*\*\*",
        re.IGNORECASE,
    ),
    re.compile(
        r"\*\*\s*Subphase\s+\w+\s+(?:shipped|landed)[^*]*\*\*",
        re.IGNORECASE,
    ),
    re.compile(
        r"\*\*\s*(?:Status|P\d+)[^*]*shipped[^*]*\*\*",
        re.IGNORECASE,
    ),
]


def find_shipped_in_body(body: str) -> tuple[bool, str]:
    """Return (is_shipped, completed_date_or_empty). Tries each pattern in order;
    the first one with a date wins for the date; the first one matching at all
    wins for the boolean."""
    is_shipped = False
    completed = ""
    for pattern in SHIPPED_BODY_PATTERNS:
        for match in pattern.finditer(body):
            is_shipped = True
            if match.groups() and match.group(1) and not completed:
                completed = match.group(1)
    return is_shipped, completed
CHECKBOX_RE = re.compile(r"^\s*-\s*\[(x| )\]", re.MULTILINE)
FILES_BLOCK_RE = re.compile(
    r"^##+\s+Files\s+to\s+(Create|Modify)\s*\n(.*?)(?=^##\s|\Z)",
    re.DOTALL | re.MULTILINE | re.IGNORECASE,
)
PATH_IN_BACKTICKS_RE = re.compile(r"`([^`\n]+\.[a-zA-Z]{1,5})`")


def parse_existing_frontmatter(text: str) -> tuple[dict, str]:
    """Return (frontmatter_dict, body_after_frontmatter). If no frontmatter,
    returns ({}, full_text)."""
    m = FRONTMATTER_RE.match(text)
    if not m:
        return {}, text
    block = m.group(1)
    body = text[m.end():]
    fm: dict = {}
    current_key = None
    for line in block.splitlines():
        if not line.strip():
            continue
        if line.startswith("  - ") and current_key == "files_touched":
            fm.setdefault("files_touched", []).append(line[4:].strip())
            continue
        if ":" in line:
            k, _, v = line.partition(":")
            k = k.strip()
            v = v.strip()
            current_key = k
            if v == "":
                fm[k] = []
            else:
                fm[k] = v
    return fm, body


def infer_status(existing_fm: dict, body: str) -> tuple[str, str, str]:
    """Returns (status, completed_date_or_empty, suggested_verification).
    Uses, in priority order:
    1. Existing frontmatter `status:`/`completed:` (with DONE→shipped migration).
    2. Body marker `**Status: ✅** ... shipped <date>` or Subphase variants.
    3. All checkboxes `[x]` in body → shipped (today as placeholder).
    suggested_verification is `partial` if existing status was PARTIAL (project
    convention: subphase-a shipped, subphase-b deferred), `pass` if shipped
    cleanly, `deferred` for status:deferred, empty otherwise."""
    raw_status = (existing_fm.get("status") or "").strip()
    if raw_status:
        normalised = raw_status.upper()
        completed_existing = (
            existing_fm.get("completed") or
            existing_fm.get("completed_date") or
            existing_fm.get("date") or ""
        ).strip()
        if normalised in {"DONE", "SHIPPED", "COMPLETED", "COMPLETE"}:
            return "shipped", completed_existing, "pass"
        if normalised == "PARTIAL":
            # Project convention: PARTIAL means subphase-a shipped, subphase-b
            # deferred (the headless-Chrome pattern). Treat as shipped at the
            # phase level; verification:partial captures the nuance.
            body_shipped, body_date = find_shipped_in_body(body)
            if body_shipped or completed_existing:
                return "shipped", body_date or completed_existing, "partial"
            return "in-progress", completed_existing, "partial"
        if normalised in {"IN-PROGRESS", "IN_PROGRESS"}:
            return "in-progress", completed_existing, ""
        if normalised in {"DEFERRED", "SKIPPED"}:
            return "deferred", completed_existing, "deferred"
        if normalised == "BLOCKED":
            return "blocked", "", ""
        if normalised in {"TODO", "NOT-STARTED", "NOT_STARTED"}:
            return "not-started", "", ""

    body_shipped, body_date = find_shipped_in_body(body)
    if body_shipped:
        return "shipped", body_date, "pass"

    boxes = CHECKBOX_RE.findall(body)
    if boxes and all(b == "x" for b in boxes):
        return "shipped", date.today().isoformat(), "pass"

    return "not-started", "", ""


def extract_files_touched(body: str) -> list[str]:
    """Pull every backticked path from `## Files to Create` and `## Files to
    Modify` sections."""
    paths: list[str] = []
    for m in FILES_BLOCK_RE.finditer(body):
        for p in PATH_IN_BACKTICKS_RE.findall(m.group(2)):
            if p not in paths and "/" in p:
                paths.append(p)
    return paths


def resolve_commit(repo_root: Path, files: list[str], completed: str) -> str:
    """Best-effort: most recent commit that touched any of `files`, optionally
    bounded by the completed date. Empty string on any failure."""
    if not files:
        return ""
    cmd = ["git", "-C", str(repo_root), "log", "-1", "--pretty=%H"]
    if completed:
        cmd += ["--until", f"{completed} 23:59:59"]
    cmd += ["--"] + files
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, check=False, timeout=10)
        return out.stdout.strip()
    except Exception:
        return ""


def render_frontmatter(
    phase: int,
    feature: str,
    status: str,
    completed: str,
    commit: str,
    files_touched: list[str],
    verification: str,
    notes: str,
) -> str:
    lines = ["---"]
    lines.append(f"phase: {phase}")
    lines.append(f"feature: {feature}")
    lines.append(f"status: {status}")
    lines.append(f"completed: {completed}")
    lines.append(f"commit: {commit}")
    if files_touched:
        lines.append("files_touched:")
        for p in files_touched:
            lines.append(f"  - {p}")
    else:
        lines.append("files_touched: []")
    lines.append(f"verification: {verification}")
    lines.append(f"notes: {notes}")
    lines.append("---")
    return "\n".join(lines) + "\n"


def migrate_phase_file(path: Path, repo_root: Path, dry_run: bool) -> str:
    text = path.read_text()
    existing_fm, body = parse_existing_frontmatter(text)

    parts = path.parts
    feature = parts[parts.index("docs") + 1]
    m_phase = re.search(r"phase-(\d+)", path.name)
    phase_n = int(m_phase.group(1)) if m_phase else 0

    status, completed, suggested_verification = infer_status(existing_fm, body)
    files_touched = (
        existing_fm.get("files_touched")
        if isinstance(existing_fm.get("files_touched"), list) and existing_fm["files_touched"]
        else extract_files_touched(body)
    )
    commit = (existing_fm.get("commit") or "").strip()
    if not commit and status == "shipped":
        commit = resolve_commit(repo_root, files_touched, completed)
    verification = (existing_fm.get("verification") or "").strip()
    if not verification:
        verification = suggested_verification or (
            "pass" if status == "shipped" else (
                "deferred" if status == "deferred" else "partial"
            )
        )
    notes = (existing_fm.get("notes") or "").strip()
    if not notes and suggested_verification == "partial" and status == "shipped":
        notes = "subphase-a shipped; subphase-b deferred (project pattern)"

    new_fm = render_frontmatter(
        phase=phase_n, feature=feature, status=status, completed=completed,
        commit=commit, files_touched=files_touched,
        verification=verification, notes=notes,
    )
    new_text = new_fm + body if body.startswith("\n") else new_fm + "\n" + body

    if dry_run:
        return f"[DRY] {path}: status={status} completed={completed!r} commit={commit[:7] or '-'} files={len(files_touched)}"
    path.write_text(new_text)
    return f"[OK ] {path}: status={status} completed={completed!r} commit={commit[:7] or '-'} files={len(files_touched)}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--docs-dir", required=True, type=Path,
                    help="Path to the docs/ root containing per-feature subdirs")
    ap.add_argument("--feature", default=None,
                    help="Limit migration to a single feature slug")
    ap.add_argument("--dry-run", action="store_true",
                    help="Print what would change; do not write")
    args = ap.parse_args()

    docs_dir: Path = args.docs_dir.resolve()
    repo_root = docs_dir
    while repo_root != repo_root.parent and not (repo_root / ".git").exists():
        repo_root = repo_root.parent

    pattern = f"{args.feature}/plan/phase-*.md" if args.feature else "*/plan/phase-*.md"
    files = sorted(p for p in docs_dir.glob(pattern) if "_meta" not in p.parts)
    if not files:
        print(f"No phase-NN.md files found under {docs_dir}/{pattern}")
        return 1

    for p in files:
        try:
            print(migrate_phase_file(p, repo_root, args.dry_run))
        except Exception as e:
            print(f"[ERR] {p}: {e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
