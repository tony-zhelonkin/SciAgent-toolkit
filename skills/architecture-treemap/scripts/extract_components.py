#!/usr/bin/env python3
"""Deterministic substrate extractor for the architecture-treemap skill.

Emits a components.json containing ONLY statically-derivable facts: the
physical file inventory, the Python intra-repo import graph (static
direct-call edges), and per-file metrics that AST, git, and optional tools
can prove. Judgment fields (logical grouping, classifications, smells,
background-knowledge edges, prune verdicts, the core/seam/removable
boundary) are NOT emitted here -- those are authored by /audit-slice and
/synthesize-audit. The honesty of this substrate (only-what-the-tools-prove)
is the trust scaffold the treemap UI renders behind every judgment badge.

CLI:
    python extract_components.py <repo-root> [--out <path>]
                                 [--since-days 90] [--lang python]

Optional dependencies (degrade gracefully when absent):
    radon       -> cyclomatic complexity (omitted if not importable)
    cloc/tokei  -> LOC counts (falls back to a stdlib Python counter)
    jsonschema  -> self-validation against components.schema.json
"""

import argparse
import ast
import datetime
import json
import os
import subprocess
import sys

TOOL_VERSION = "architecture-treemap 0.1.0"
SCHEMA_VERSION = "1.0"

# Extension -> physical-component kind. Location can override this (a file
# under a tests/ dir is a test regardless of extension).
EXTENSION_KIND = {
    ".py": "module",
    ".pyi": "module",
    ".r": "module",
    ".sh": "module",
    ".yaml": "config",
    ".yml": "config",
    ".toml": "config",
    ".cfg": "config",
    ".ini": "config",
    ".json": "config",
    ".md": "doc",
    ".rst": "doc",
    ".txt": "doc",
    ".png": "asset",
    ".jpg": "asset",
    ".jpeg": "asset",
    ".svg": "asset",
    ".gif": "asset",
    ".css": "asset",
    ".html": "asset",
}

# Directories never worth inventorying even if git happens to track them.
VENDOR_DIR_NAMES = {
    "node_modules",
    "dist",
    "build",
    ".venv",
    "venv",
    "__pycache__",
    ".git",
    ".egg-info",
    "vendor",
    "site-packages",
}


# ---------------------------------------------------------------------------
# File discovery
# ---------------------------------------------------------------------------

def list_repo_files(repo_root):
    """Return repo-relative paths of tracked + untracked-but-not-ignored files.

    Uses `git ls-files` so .gitignore is honored for free. Falls back to a
    plain walk (skipping vendor dirs) when the tree is not a git repo.
    """
    git_files = run_git(
        repo_root,
        ["ls-files", "--cached", "--others", "--exclude-standard"],
    )
    if git_files is not None:
        paths = [line for line in git_files.splitlines() if line.strip()]
        return [p for p in paths if not is_vendored(p)]
    return walk_files(repo_root)


def walk_files(repo_root):
    """Non-git fallback: walk the tree, skipping vendored directories."""
    collected = []
    for current_dir, subdirs, filenames in os.walk(repo_root):
        subdirs[:] = [d for d in subdirs if d not in VENDOR_DIR_NAMES]
        for name in filenames:
            absolute = os.path.join(current_dir, name)
            relative = os.path.relpath(absolute, repo_root)
            if not is_vendored(relative):
                collected.append(relative)
    return sorted(collected)


def is_vendored(relative_path):
    """True if any path segment names a vendored/build directory."""
    parts = relative_path.replace("\\", "/").split("/")
    for part in parts:
        if part in VENDOR_DIR_NAMES:
            return True
        if part.endswith(".egg-info"):
            return True
    return False


def classify_kind(relative_path):
    """Map a path to a physical-component kind. Location beats extension."""
    lowered = relative_path.lower()
    parts = lowered.split("/")
    name = parts[-1]

    is_in_test_dir = any(part in ("test", "tests", "testing") for part in parts)
    looks_like_test = name.startswith("test_") or name.endswith("_test.py")
    if is_in_test_dir or looks_like_test:
        return "test"

    _, ext = os.path.splitext(name)
    return EXTENSION_KIND.get(ext, "other")


def kebab_id(relative_path):
    """Kebab-case the path into a schema-legal id (^[a-z0-9-]+$)."""
    lowered = relative_path.lower()
    result_chars = []
    for ch in lowered:
        if ch.isalnum():
            result_chars.append(ch)
        else:
            result_chars.append("-")
    collapsed = "".join(result_chars)
    while "--" in collapsed:
        collapsed = collapsed.replace("--", "-")
    return collapsed.strip("-")


def build_id_map(repo_files):
    """Map repo-relative path -> unique schema-legal id.

    The naive kebab-case can collide (e.g. __main__.py and main.py both
    collapse to the same id). Edges reference these ids, so they must be
    unique and stable. Files are processed in sorted order so the
    disambiguation suffix (-2, -3, ...) is deterministic across runs.
    """
    id_map = {}
    used = {}
    for relative in sorted(repo_files):
        base = kebab_id(relative) or "file"
        if base not in used:
            used[base] = 1
            id_map[relative] = base
        else:
            used[base] += 1
            id_map[relative] = "{}-{}".format(base, used[base])
    return id_map


# ---------------------------------------------------------------------------
# LOC counting
# ---------------------------------------------------------------------------

def count_loc(absolute_path, kind):
    """Count lines of code for a file.

    Python files: non-blank, non-comment-only lines. Everything else: raw
    line count. cloc/tokei are not depended on; this stdlib counter is the
    deterministic floor and is recorded as the LOC tool used.
    """
    try:
        with open(absolute_path, "r", encoding="utf-8", errors="replace") as handle:
            lines = handle.readlines()
    except OSError:
        return 0

    if kind in ("module", "test") and absolute_path.endswith((".py", ".pyi")):
        return count_python_loc(lines)
    return sum(1 for line in lines if line.strip())


def count_python_loc(lines):
    """Non-blank, non-comment-only physical lines."""
    counted = 0
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            continue
        counted += 1
    return counted


# ---------------------------------------------------------------------------
# Python import graph
# ---------------------------------------------------------------------------

def find_python_roots(repo_root, python_files):
    """Discover import roots so module names resolve to repo files.

    Handles src-layout: if files live under src/<pkg>/, then src/ is an
    import root (so `import pkg.foo` resolves). The repo root is always a
    root. Returns repo-relative root prefixes ("" means the repo root).
    """
    roots = {""}
    for relative in python_files:
        parts = relative.split("/")
        if parts[0] == "src":
            roots.add("src")
    return roots


def module_name_for_file(relative_path, roots):
    """Return the dotted module name a .py file is importable as, given roots.

    Picks the longest matching root prefix (src beats repo-root) so a
    src-layout file maps to its package-qualified name.
    """
    if not relative_path.endswith(".py"):
        return None
    without_ext = relative_path[: -len(".py")]
    parts = without_ext.split("/")
    if parts[-1] == "__init__":
        parts = parts[:-1]

    best_module = None
    best_root_depth = -1
    for root in roots:
        root_parts = [p for p in root.split("/") if p]
        depth = len(root_parts)
        if parts[:depth] == root_parts and depth > best_root_depth:
            remaining = parts[depth:]
            best_module = ".".join(remaining)
            best_root_depth = depth
    return best_module or None


def build_module_index(python_files, roots):
    """Map dotted module name -> repo-relative file path."""
    index = {}
    for relative in python_files:
        module = module_name_for_file(relative, roots)
        if module:
            index[module] = relative
    return index


def collect_imports(absolute_path, source_module):
    """Parse a .py file and yield (candidate_names, line) tuples.

    `candidate_names` is a list of dotted module names to try, most-specific
    first. For `import a.b` it is `["a.b"]`. For `from pkg import a` it is
    `["pkg.a", "pkg"]` -- `a` may be a submodule (resolves to pkg/a.py) or a
    symbol inside pkg/__init__.py, so both candidates are offered and the
    first that maps to a repo file wins. Relative imports are resolved against
    the source module's package. Returns None on a syntax error so one bad
    file never crashes the run.
    """
    try:
        with open(absolute_path, "r", encoding="utf-8", errors="replace") as handle:
            tree = ast.parse(handle.read(), filename=absolute_path)
    except (SyntaxError, ValueError):
        return None

    package_parts = source_module.split(".")[:-1] if source_module else []
    imports = []

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                imports.append(([alias.name], node.lineno))
        elif isinstance(node, ast.ImportFrom):
            base = resolve_import_from_base(node, package_parts)
            if base is None and node.level == 0:
                continue
            for alias in node.names:
                candidates = import_from_candidates(base, alias.name)
                if candidates:
                    imports.append((candidates, node.lineno))
    return imports


def resolve_import_from_base(node, package_parts):
    """Resolve the base package of an ImportFrom to a dotted name (or '')."""
    if node.level == 0:
        return node.module or ""
    # Relative import: walk up `level` packages from the current package.
    base_parts = package_parts[: len(package_parts) - (node.level - 1)]
    if node.module:
        return ".".join(base_parts + node.module.split("."))
    return ".".join(base_parts)


def import_from_candidates(base, imported_symbol):
    """Build candidate dotted names for `from <base> import <symbol>`.

    Offers `base.symbol` (symbol-is-submodule) first, then bare `base`
    (symbol-is-a-name-in-base's-__init__). Drops empties.
    """
    candidates = []
    if base and imported_symbol and imported_symbol != "*":
        candidates.append("{}.{}".format(base, imported_symbol))
    if base:
        candidates.append(base)
    return candidates


def resolve_to_repo_module(candidate_names, module_index):
    """Map candidate dotted names to a repo module, or None if all external.

    Tries each candidate in order; for each, also trims trailing segments
    (a `from pkg.mod import thing` may import a symbol, not a submodule).
    """
    for name in candidate_names:
        parts = name.split(".")
        while parts:
            candidate = ".".join(parts)
            if candidate in module_index:
                return candidate
            parts = parts[:-1]
    return None


def build_import_graph(repo_root, python_files, module_index, roots, id_map):
    """Build static direct-call edges plus fan-in / fan-out counts.

    Returns (edges, fan_out_by_file, fan_in_by_file, parse_failures) where
    edges are dicts keyed by physical-component id and fan counts are keyed
    by repo-relative path.
    """
    edges = []
    seen_edges = set()
    fan_out = {relative: set() for relative in python_files}
    fan_in = {relative: set() for relative in python_files}
    parse_failures = 0

    for relative in python_files:
        absolute = os.path.join(repo_root, relative)
        source_module = module_name_for_file(relative, roots)
        imports = collect_imports(absolute, source_module)
        if imports is None:
            parse_failures += 1
            continue

        for candidate_names, line in imports:
            target_module = resolve_to_repo_module(candidate_names, module_index)
            if target_module is None:
                continue
            target_file = module_index[target_module]
            if target_file == relative:
                continue  # self-import (e.g. package __init__); ignore.

            fan_out[relative].add(target_file)
            fan_in[target_file].add(relative)

            edge_key = (relative, target_file)
            if edge_key in seen_edges:
                continue
            seen_edges.add(edge_key)
            edges.append(
                {
                    "from": id_map[relative],
                    "to": id_map[target_file],
                    "type": "direct-call",
                    "evidence_class": "static",
                    "evidence": "{}:{} -> import {}".format(
                        relative, line, target_module
                    ),
                }
            )

    edges.sort(key=lambda e: (e["from"], e["to"]))
    return edges, fan_out, fan_in, parse_failures


# ---------------------------------------------------------------------------
# Git-derived metrics
# ---------------------------------------------------------------------------

def run_git(repo_root, args):
    """Run a git subcommand; return stdout text, or None if git is unavailable."""
    try:
        completed = subprocess.run(
            ["git", "-C", repo_root] + args,
            capture_output=True,
            text=True,
            check=False,
        )
    except (OSError, FileNotFoundError):
        return None
    if completed.returncode != 0:
        return None
    return completed.stdout


def git_head_sha(repo_root):
    """Current HEAD SHA, or 'unknown' when not a git repo."""
    output = run_git(repo_root, ["rev-parse", "HEAD"])
    if output is None:
        return "unknown"
    return output.strip() or "unknown"


def churn_for_file(repo_root, relative_path, since_days):
    """Count commits touching a file in the last `since_days`, or None."""
    output = run_git(
        repo_root,
        [
            "log",
            "--since={}.days".format(since_days),
            "--format=%H",
            "--",
            relative_path,
        ],
    )
    if output is None:
        return None
    return sum(1 for line in output.splitlines() if line.strip())


# ---------------------------------------------------------------------------
# Cyclomatic complexity (optional, via radon)
# ---------------------------------------------------------------------------

def load_radon():
    """Return the radon cc_visit callable, or None if radon is absent."""
    try:
        from radon.complexity import cc_visit
    except ImportError:
        return None
    return cc_visit


def cyclomatic_for_file(absolute_path, cc_visit):
    """Mean cyclomatic complexity across blocks in a Python file, or None."""
    if cc_visit is None:
        return None
    try:
        with open(absolute_path, "r", encoding="utf-8", errors="replace") as handle:
            blocks = cc_visit(handle.read())
    except (SyntaxError, ValueError, OSError):
        return None
    if not blocks:
        return None
    total = sum(block.complexity for block in blocks)
    mean = total / len(blocks)
    return round(mean, 2)


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------

def compute_repo_test_ratio(physical_records):
    """Repo-level test LOC / source LOC.

    test_ratio is fundamentally a logical-grouping rollup (the synth step
    owns it per-component). At the physical level the only honest figure is
    the whole-repo ratio, which is attached to test files' metrics so the
    badge still reflects a real number. Source LOC counts module files only.
    """
    test_loc = sum(r["size_loc"] for r in physical_records if r["kind"] == "test")
    source_loc = sum(r["size_loc"] for r in physical_records if r["kind"] == "module")
    if source_loc == 0:
        return None
    return round(test_loc / source_loc, 3)


def build_physical_components(
    repo_root,
    repo_files,
    fan_out,
    fan_in,
    cc_visit,
    since_days,
    metric_tools,
    id_map,
    include_docs=False,
):
    """Assemble physical_components[] with their metrics blocks.

    Documentation files (kind "doc" — .md/.rst/.txt by extension) are excluded
    by default because they dominate the tile count and obscure architectural
    signal.  Pass include_docs=True to restore them.
    """
    records = []
    churn_available = False

    for relative in sorted(repo_files):
        absolute = os.path.join(repo_root, relative)
        kind = classify_kind(relative)

        # Skip doc files unless the caller explicitly opted in.
        if not include_docs and kind == "doc":
            continue

        loc = count_loc(absolute, kind)

        metrics = {"loc": loc}

        if relative in fan_out:
            metrics["fan_out"] = len(fan_out[relative])
        if relative in fan_in:
            metrics["fan_in"] = len(fan_in[relative])

        complexity = cyclomatic_for_file(absolute, cc_visit)
        if complexity is not None:
            metrics["cyclomatic"] = complexity

        churn = churn_for_file(repo_root, relative, since_days)
        if churn is not None:
            metrics["churn_90d"] = churn
            churn_available = True

        records.append(
            {
                "id": id_map[relative],
                "path": relative,
                "size_loc": loc,
                "logical_owners": [],
                "kind": kind,
                "metrics": metrics,
            }
        )

    if churn_available and "git-log" not in metric_tools:
        metric_tools.append("git-log")

    attach_test_ratio(records)
    return records


def attach_test_ratio(records):
    """Attach the repo-level test/code ratio to each test file's metrics."""
    ratio = compute_repo_test_ratio(records)
    if ratio is None:
        return
    for record in records:
        if record["kind"] == "test":
            record["metrics"]["test_ratio"] = ratio


def build_manifest(repo_root, args):
    """Produce the full deterministic-substrate manifest dict."""
    repo_files = list_repo_files(repo_root)
    python_files = sorted(f for f in repo_files if f.endswith(".py"))
    id_map = build_id_map(repo_files)

    roots = find_python_roots(repo_root, python_files)
    module_index = build_module_index(python_files, roots)
    edges, fan_out, fan_in, parse_failures = build_import_graph(
        repo_root, python_files, module_index, roots, id_map
    )

    metric_tools = ["ast"]
    cc_visit = load_radon()
    if cc_visit is not None:
        metric_tools.append("radon")
    metric_tools.append("stdlib-loc")

    physical = build_physical_components(
        repo_root,
        repo_files,
        fan_out,
        fan_in,
        cc_visit,
        args.since_days,
        metric_tools,
        id_map,
        include_docs=getattr(args, 'include_docs', False),
    )

    today = datetime.date.today().isoformat()
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "audit_date": today,
        "project": os.path.basename(os.path.abspath(repo_root)),
        "git_sha": git_head_sha(repo_root),
        "snapshot_id": "{}".format(today),
        "extractor": {
            "tool_version": TOOL_VERSION,
            "extracted_at": datetime.datetime.now().isoformat(timespec="seconds"),
            "metric_tools": metric_tools,
        },
        "logical_components": [],
        "physical_components": physical,
        "edges": edges,
    }

    summary = {
        "files_scanned": len(repo_files),
        "python_files": len(python_files),
        "edges": len(edges),
        "parse_failures": parse_failures,
        "metric_tools": metric_tools,
    }
    return manifest, summary


# ---------------------------------------------------------------------------
# Self-validation (best-effort)
# ---------------------------------------------------------------------------

def validate_manifest(manifest, schema_path):
    """Validate against the schema if jsonschema is importable. Best-effort."""
    try:
        import jsonschema
    except ImportError:
        return "skipped (jsonschema not installed)"
    try:
        with open(schema_path, "r", encoding="utf-8") as handle:
            schema = json.load(handle)
    except OSError:
        return "skipped (schema not found at {})".format(schema_path)
    try:
        jsonschema.validate(instance=manifest, schema=schema)
    except jsonschema.ValidationError as error:
        return "FAILED: {}".format(error.message)
    return "passed"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Extract the deterministic architecture substrate "
        "(physical files, static import edges, metrics) to a components.json."
    )
    parser.add_argument("repo_root", help="Path to the repository root to scan.")
    parser.add_argument(
        "--out",
        default=None,
        help="Output path for components.json (default: <repo-root>/components.json).",
    )
    parser.add_argument(
        "--since-days",
        type=int,
        default=90,
        help="Churn window in days (default: 90).",
    )
    parser.add_argument(
        "--lang",
        default="python",
        help="Source language for the import graph (only 'python' is supported).",
    )
    parser.add_argument(
        "--include-docs",
        action="store_true",
        default=False,
        dest="include_docs",
        help=(
            "Include documentation files (kind 'doc': .md/.rst/.txt) in the "
            "physical_components output.  By default they are excluded because "
            "they typically dominate the tile count and obscure architectural signal."
        ),
    )
    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)

    repo_root = os.path.abspath(args.repo_root)
    if not os.path.isdir(repo_root):
        sys.stderr.write("error: not a directory: {}\n".format(repo_root))
        return 2

    if args.lang != "python":
        sys.stderr.write(
            "warning: --lang={} unsupported; import graph is Python-only.\n".format(
                args.lang
            )
        )

    manifest, summary = build_manifest(repo_root, args)

    out_path = args.out or os.path.join(repo_root, "components.json")
    with open(out_path, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=False)
        handle.write("\n")

    schema_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..",
        "components.schema.json",
    )
    validation = validate_manifest(manifest, schema_path)

    sys.stderr.write(
        "extract_components: scanned {files} files ({py} python), "
        "{edges} static edges, {fail} unparseable file(s) skipped.\n".format(
            files=summary["files_scanned"],
            py=summary["python_files"],
            edges=summary["edges"],
            fail=summary["parse_failures"],
        )
    )
    sys.stderr.write(
        "metric tools used: {tools}\n".format(tools=", ".join(summary["metric_tools"]))
    )
    sys.stderr.write("schema validation: {}\n".format(validation))
    sys.stderr.write("wrote {}\n".format(out_path))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
