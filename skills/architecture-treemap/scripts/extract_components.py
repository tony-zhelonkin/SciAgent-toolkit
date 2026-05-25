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

# The metric model is declared once in metric_registry.py (single source of
# truth). The extractor derives WHAT to compute from it instead of hardcoding
# a compute-list; see build_physical_components and metric_registry.py.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import metric_registry  # noqa: E402
import pedagogy_registry  # noqa: E402

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

# The skill's own output lands under an `architecture-audit` directory. A tool
# must not render its own output as architecture, so any path with this segment
# is excluded by default (committed prior snapshots included). The current
# run's --out directory is excluded the same way (computed in build_manifest).
# Override with --include-audit-output when auditing the skill repo itself.
AUDIT_OUTPUT_SEGMENT = "architecture-audit"


# ---------------------------------------------------------------------------
# Physical-component inventory filter ("is this plausibly architectural?")
# ---------------------------------------------------------------------------
#
# `git ls-files --exclude-standard` already honors .gitignore for UNTRACKED
# files. This filter is the complement: it handles TRACKED scaffolding that
# gitignore cannot — license files, lockfiles, root-level tooling/dotfiles,
# harness dirs, and data files — that are real tracked files but are NOT
# architecture. Without it, LICENSE / package.json / .gitignore / .mcp.json /
# .claude/* / docs/**/*.csv all leak in as tiles and drown the signal.
#
# The rule honored by `is_nonsource_scaffolding`:
#   * Always-drop sets below (license / lock / dotfile / data-ext) drop the
#     file wherever it lives.
#   * EVERYTHING under a top-level `docs/` dir is dropped by PATH-PREFIX (this
#     supersedes the old extension-only doc filter and fixes the .csv-under-docs
#     leak that extension matching missed).
#   * Harness/tooling dirs (.claude, .github, ...) are dropped wherever they sit.
#   * A `config`/`other`/`asset`-kind file is KEPT only when it lives under a
#     source root (src/ or the package dir) AND is not in a drop set; otherwise
#     it is dropped. This keeps src-resident schemas/contracts (the
#     plugin_payload.schema.json wire contract) while dropping root-level
#     tooling JSON (package.json, .mcp.json).
#   * `module`/`test`-kind files are always kept (verify_*.py and one-off
#     scripts under scripts/ are real modules — do not over-filter them).
# Override the whole filter with --include-nonsource.

# License / legal files, by basename stem (case-insensitive), wherever they sit.
LICENSE_BASENAMES = ("license", "copying", "notice", "authors", "patents")

# Lockfiles — generated, never architectural.
LOCKFILE_NAMES = {
    "package-lock.json",
    "yarn.lock",
    "pnpm-lock.yaml",
    "poetry.lock",
    "pipfile.lock",
    "cargo.lock",
    "composer.lock",
    "gemfile.lock",
}
LOCKFILE_EXTS = {".lock"}

# Root-level tooling / dotfiles. Matched only at the REPO ROOT (depth 0) so a
# legitimately-source-resident file of the same name is not collateral.
ROOT_TOOLING_NAMES = {
    ".gitignore",
    ".gitattributes",
    ".mcp.json",
    ".editorconfig",
    ".pre-commit-config.yaml",
    ".dockerignore",
    ".npmignore",
    ".prettierrc",
    ".eslintrc",
    ".eslintrc.json",
    ".eslintrc.js",
    ".flake8",
    ".python-version",
    "package.json",
    "package-lock.json",
    "pyproject.toml",
    "setup.cfg",
    "setup.py",
    "tox.ini",
    "makefile",
    "dockerfile",
    "requirements.txt",
    "requirements-dev.txt",
}

# Harness / editor / CI directories: drop everything beneath them, wherever
# they appear (a .claude/ may be nested). These are tooling, not architecture.
TOOLING_DIR_NAMES = {
    ".claude",
    ".agents",
    ".sciagent",
    ".github",
    ".gitlab",
    ".vscode",
    ".idea",
    ".devcontainer",
    ".circleci",
}

# Data files by extension — payloads, not code.
DATA_EXTS = {
    ".csv", ".tsv", ".parquet", ".feather", ".h5", ".h5ad", ".hdf5",
    ".rds", ".rdata", ".npy", ".npz", ".pkl", ".pickle", ".joblib",
    ".arrow", ".db", ".sqlite", ".sqlite3",
}

# Source-root prefixes under which a config/other/asset file is "architectural"
# (e.g. a wire-contract schema living next to the code it governs). The package
# dir is also treated as a source root, discovered at runtime.
SOURCE_ROOT_PREFIXES = ("src/",)

# Documentation directory at the repo root. Dropped by path-prefix so .csv (and
# anything else) living under docs/ never leaks. This is the prefix-based
# superset of the old extension-only doc exclusion.
DOC_DIR_PREFIX = "docs/"


def _basename_stem(relative_path):
    """Lowercase (stem, ext) of a path's final segment."""
    name = relative_path.replace("\\", "/").split("/")[-1].lower()
    stem, ext = os.path.splitext(name)
    return stem, ext


def _is_root_level(relative_path):
    """True if the path sits directly at the repo root (no directory part)."""
    return "/" not in relative_path.replace("\\", "/")


def _under_source_root(relative_path, package_dirs):
    """True if the path lives under a recognised source root.

    A source root is `src/`, or any discovered top-level package directory (a
    dir containing an __init__.py). This is what distinguishes a src-resident
    schema (KEEP) from a root-level tooling JSON (DROP).
    """
    normalised = relative_path.replace("\\", "/")
    for prefix in SOURCE_ROOT_PREFIXES:
        if normalised.startswith(prefix):
            return True
    top = normalised.split("/")[0]
    return top in package_dirs


def discover_package_dirs(repo_files):
    """Top-level directories that look like a Python package (have __init__.py).

    Lets a flat-layout package (e.g. `mypkg/__init__.py` at the repo root)
    count as a source root so its in-package schemas/configs are kept.
    """
    package_dirs = set()
    for relative in repo_files:
        normalised = relative.replace("\\", "/")
        parts = normalised.split("/")
        if len(parts) == 2 and parts[1] == "__init__.py":
            package_dirs.add(parts[0])
    return package_dirs


def is_nonsource_scaffolding(relative_path, kind, package_dirs):
    """True if a TRACKED file is scaffolding, not architecture (default-drop).

    See the module-level commentary for the full rule. Returns True to EXCLUDE.
    Source modules and tests are never excluded here.
    """
    normalised = relative_path.replace("\\", "/")
    stem, ext = _basename_stem(normalised)

    # Everything under a top-level docs/ dir, by path-prefix (catches .csv etc).
    if normalised.startswith(DOC_DIR_PREFIX):
        return True

    # Anything beneath a harness/tooling dir, wherever it appears.
    parts = normalised.split("/")
    if any(part in TOOLING_DIR_NAMES for part in parts):
        return True

    # License / legal files, wherever they sit.
    if stem in LICENSE_BASENAMES:
        return True

    # Lockfiles.
    if normalised.split("/")[-1].lower() in LOCKFILE_NAMES or ext in LOCKFILE_EXTS:
        return True

    # Data files by extension.
    if ext in DATA_EXTS:
        return True

    # Root-level tooling / dotfiles (only at depth 0; a same-named file deeper
    # in a source tree is judged by the kind rule below).
    if _is_root_level(normalised) and normalised.lower() in ROOT_TOOLING_NAMES:
        return True

    # Real source: modules and tests are always architectural.
    if kind in ("module", "test"):
        return False

    # config / asset / other: keep ONLY when it lives under a source root (a
    # wire-contract schema next to its code). Otherwise it is scaffolding.
    if kind in ("config", "asset", "other"):
        return not _under_source_root(normalised, package_dirs)

    # doc and anything unforeseen: drop (docs are handled by build_physical too).
    return True


# ---------------------------------------------------------------------------
# File discovery
# ---------------------------------------------------------------------------

def list_repo_files(repo_root, exclude_audit=True, out_dir_rel=None):
    """Return repo-relative paths of tracked + untracked-but-not-ignored files.

    Uses `git ls-files` so .gitignore is honored for free. Falls back to a
    plain walk (skipping vendor dirs) when the tree is not a git repo. Audit
    output (this skill's own snapshots) is excluded unless exclude_audit=False.
    """
    git_files = run_git(
        repo_root,
        ["ls-files", "--cached", "--others", "--exclude-standard"],
    )
    if git_files is not None:
        paths = [line for line in git_files.splitlines() if line.strip()]
        return [
            p for p in paths
            if not is_vendored(p)
            and not (exclude_audit and is_audit_output(p, out_dir_rel))
        ]
    return walk_files(repo_root, exclude_audit=exclude_audit, out_dir_rel=out_dir_rel)


def walk_files(repo_root, exclude_audit=True, out_dir_rel=None):
    """Non-git fallback: walk the tree, skipping vendored + audit-output dirs."""
    collected = []
    for current_dir, subdirs, filenames in os.walk(repo_root):
        subdirs[:] = [d for d in subdirs if d not in VENDOR_DIR_NAMES]
        for name in filenames:
            absolute = os.path.join(current_dir, name)
            relative = os.path.relpath(absolute, repo_root)
            if is_vendored(relative):
                continue
            if exclude_audit and is_audit_output(relative, out_dir_rel):
                continue
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


def is_audit_output(relative_path, out_dir_rel=None):
    """True if the path is one of the skill's own output artifacts.

    Two signals, unioned: (1) any path segment equals AUDIT_OUTPUT_SEGMENT,
    which catches committed prior snapshots wherever they live; (2) the path
    sits under the resolved --out directory for this run, which catches the
    current destination even when redirected elsewhere. A tool must not render
    its own output as architecture.
    """
    normalised = relative_path.replace("\\", "/")
    parts = normalised.split("/")
    if AUDIT_OUTPUT_SEGMENT in parts:
        return True
    if out_dir_rel:
        out_norm = out_dir_rel.replace("\\", "/").strip("/")
        if out_norm and (normalised == out_norm or normalised.startswith(out_norm + "/")):
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
# Cyclic-dependency detection (Tarjan SCC over the static import graph)
# ---------------------------------------------------------------------------

def find_cycles(edges):
    """Return import cycles as a list of {id, members[]} via Tarjan's SCC.

    Deterministic: runs over the static direct-call edges already built by
    build_import_graph. A strongly-connected component of size > 1 is an import
    cycle — every member can reach every other, so no member can be understood
    or changed in isolation. Teaches the Acyclic Dependencies Principle. A
    single self-loop is ignored (build_import_graph already drops self-imports).

    Iterative (explicit stack) so a deep graph cannot blow the recursion limit.
    Node and member ordering is sorted so the output is byte-stable across runs.
    """
    adjacency = {}
    for edge in edges:
        adjacency.setdefault(edge["from"], set()).add(edge["to"])
        adjacency.setdefault(edge["to"], set())
    nodes = sorted(adjacency.keys())

    index_of = {}
    lowlink = {}
    on_stack = {}
    stack = []
    next_index = [0]
    sccs = []

    def strongconnect(start):
        # Explicit work stack of (node, neighbour-iterator) frames.
        work = [(start, iter(sorted(adjacency[start])))]
        index_of[start] = lowlink[start] = next_index[0]
        next_index[0] += 1
        stack.append(start)
        on_stack[start] = True

        while work:
            node, neighbours = work[-1]
            advanced = False
            for nxt in neighbours:
                if nxt not in index_of:
                    index_of[nxt] = lowlink[nxt] = next_index[0]
                    next_index[0] += 1
                    stack.append(nxt)
                    on_stack[nxt] = True
                    work.append((nxt, iter(sorted(adjacency[nxt]))))
                    advanced = True
                    break
                if on_stack.get(nxt):
                    lowlink[node] = min(lowlink[node], index_of[nxt])
            if advanced:
                continue
            # All neighbours processed: maybe close an SCC rooted at `node`.
            if lowlink[node] == index_of[node]:
                component = []
                while True:
                    w = stack.pop()
                    on_stack[w] = False
                    component.append(w)
                    if w == node:
                        break
                if len(component) > 1:
                    sccs.append(sorted(component))
            work.pop()
            if work:
                parent = work[-1][0]
                lowlink[parent] = min(lowlink[parent], lowlink[node])

    for node in nodes:
        if node not in index_of:
            strongconnect(node)

    sccs.sort(key=lambda members: (len(members), members))
    cycles = []
    for i, members in enumerate(sccs):
        cycles.append({"id": "cycle-{}".format(i + 1), "members": members})
    return cycles


def annotate_cycle_membership(physical, cycles):
    """Stamp each physical component in a cycle with cycle_id + cycle_size.

    A per-component flag the renderer can read directly (so a tile can show "in
    a cycle of size N") without re-deriving the SCC client-side. Deterministic;
    mutates the physical records in place.
    """
    member_to_cycle = {}
    for cycle in cycles:
        for member in cycle["members"]:
            member_to_cycle[member] = cycle
    for record in physical:
        cycle = member_to_cycle.get(record["id"])
        if cycle is not None:
            record["cycle_id"] = cycle["id"]
            record["cycle_size"] = len(cycle["members"])


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

def build_physical_components(
    repo_root,
    repo_files,
    fan_out,
    fan_in,
    cc_visit,
    since_days,
    metric_tools,
    id_map,
    coverage_by_file,
    package_dirs,
    include_docs=False,
    include_nonsource=False,
):
    """Assemble physical_components[] with their metrics blocks.

    Two layers of exclusion, both overridable:
      * The non-source scaffolding filter (is_nonsource_scaffolding): license
        files, lockfiles, root-level tooling, harness dirs, data files, and
        everything under docs/. These are tracked files that gitignore cannot
        catch but that are NOT architecture. Overridden by include_nonsource.
      * Documentation files (kind "doc") are also excluded by default because
        they dominate the tile count and obscure signal. Overridden by
        include_docs. (docs/ is already covered by the scaffolding filter; this
        catches stray .md/.rst/.txt outside docs/.)

    Base metrics are computed here (the extractor owns the AST/git handles
    those formulas need). Derived metrics (e.g. refactor_pressure) are then
    layered on via the metric registry — no per-metric branch lives here, so a
    new derived metric never touches this function.
    """
    records = []
    churn_available = False
    excluded_nonsource = 0

    for relative in sorted(repo_files):
        absolute = os.path.join(repo_root, relative)
        kind = classify_kind(relative)

        # Drop non-source scaffolding (license/lock/tooling/data/docs) unless
        # the caller opts in. gitignore handled untracked files; this handles
        # tracked scaffolding that is not architecture.
        if not include_nonsource and is_nonsource_scaffolding(relative, kind, package_dirs):
            excluded_nonsource += 1
            continue

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

        # Per-module coverage signal, derived from the static test->module
        # import graph (a module imported by >=1 test file is "covered").
        # Deterministic; no fragile filename heuristic. test_ratio carries the
        # magnitude (covered modules get the repo-level ratio; uncovered get 0).
        if kind == "module" and relative in coverage_by_file:
            metrics["test_ratio"] = coverage_by_file[relative]

        # Derived metrics (registry-driven). Adds only what the formulas can
        # honestly compute from the base metrics present.
        metrics = metric_registry.compute_derived(metrics)

        # Degraded "look here first" proxy when radon is absent (so the real,
        # cyclomatic-based refactor_pressure could not be computed). Emitted
        # under a distinct key + label so it is never conflated with the real
        # metric. When radon is present this returns None and nothing is added.
        if cc_visit is None:
            proxy = metric_registry.refactor_pressure_loc_proxy(metrics)
            if proxy is not None:
                metrics["refactor_pressure_loc_proxy"] = proxy

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

    return records, excluded_nonsource


def compute_module_coverage(repo_root, repo_files, fan_in):
    """Per-module test coverage, derived from the static import graph.

    A module is "covered" if at least one test-kind file imports it (a real
    fan_in edge whose source is a test file) — fully deterministic, no fragile
    test-filename heuristic. The repo-level test/source LOC ratio supplies the
    magnitude so the existing `test_ratio` field carries a real number: covered
    modules get that ratio, uncovered modules get 0.0 so the renderer's "low
    tests" badge lights honestly on whatever no test touches.

    Returns {repo_relative_path: test_ratio_value} for module files only.
    """
    test_paths = {p for p in repo_files if classify_kind(p) == "test"}
    module_paths = [p for p in repo_files if classify_kind(p) == "module"]

    test_loc = sum(
        count_loc(os.path.join(repo_root, p), "test") for p in test_paths
    )
    source_loc = sum(
        count_loc(os.path.join(repo_root, p), "module") for p in module_paths
    )
    repo_ratio = round(test_loc / source_loc, 3) if source_loc else 0.0

    coverage = {}
    for module_path in module_paths:
        importers = fan_in.get(module_path, set())
        covered = any(imp in test_paths for imp in importers)
        coverage[module_path] = repo_ratio if covered else 0.0
    return coverage


def relative_out_dir(repo_root, out_path):
    """Repo-relative directory of the --out target, or None if outside repo."""
    if not out_path:
        return None
    out_dir = os.path.dirname(os.path.abspath(out_path))
    try:
        rel = os.path.relpath(out_dir, repo_root)
    except ValueError:
        return None
    if rel.startswith(".."):
        return None  # output is outside the repo; nothing to exclude
    return rel


def build_manifest(repo_root, args):
    """Produce the full deterministic-substrate manifest dict."""
    exclude_audit = not getattr(args, "include_audit_output", False)
    out_dir_rel = relative_out_dir(repo_root, getattr(args, "out", None))
    repo_files = list_repo_files(
        repo_root, exclude_audit=exclude_audit, out_dir_rel=out_dir_rel
    )
    python_files = sorted(f for f in repo_files if f.endswith(".py"))
    id_map = build_id_map(repo_files)

    roots = find_python_roots(repo_root, python_files)
    module_index = build_module_index(python_files, roots)
    edges, fan_out, fan_in, parse_failures = build_import_graph(
        repo_root, python_files, module_index, roots, id_map
    )

    metric_tools = ["ast"]
    cc_visit = load_radon()
    radon_present = cc_visit is not None
    if radon_present:
        metric_tools.append("radon")
    metric_tools.append("stdlib-loc")

    # Honest degradation: stamp which metrics could NOT be computed because an
    # optional tool was absent, so the renderer can grey the relevant lens with
    # a reason rather than show an empty panel that reads as "nothing here".
    unavailable_metrics = []
    degradation_notes = {}
    if not radon_present:
        unavailable_metrics.extend(["cyclomatic", "refactor_pressure"])
        degradation_notes["cyclomatic"] = "radon not installed (pip install radon)"
        degradation_notes["refactor_pressure"] = (
            "depends on cyclomatic; install radon. A degraded LOC-proxy "
            "(refactor_pressure_loc_proxy) is emitted instead."
        )

    coverage_by_file = compute_module_coverage(repo_root, repo_files, fan_in)
    package_dirs = discover_package_dirs(repo_files)

    physical, excluded_nonsource = build_physical_components(
        repo_root,
        repo_files,
        fan_out,
        fan_in,
        cc_visit,
        args.since_days,
        metric_tools,
        id_map,
        coverage_by_file,
        package_dirs,
        include_docs=getattr(args, 'include_docs', False),
        include_nonsource=getattr(args, 'include_nonsource', False),
    )

    # Cyclic-dependency detection over the static import graph (deterministic).
    cycles = find_cycles(edges)
    annotate_cycle_membership(physical, cycles)

    today = datetime.date.today().isoformat()
    extractor_block = {
        "tool_version": TOOL_VERSION,
        "extracted_at": datetime.datetime.now().isoformat(timespec="seconds"),
        "metric_tools": metric_tools,
    }
    if unavailable_metrics:
        extractor_block["unavailable_metrics"] = unavailable_metrics
        extractor_block["degradation_notes"] = degradation_notes

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "audit_date": today,
        "project": os.path.basename(os.path.abspath(repo_root)),
        "git_sha": git_head_sha(repo_root),
        "snapshot_id": "{}".format(today),
        "extractor": extractor_block,
        # Self-describing metric model: the renderer derives badges + panel rows
        # from this block, so the rendered HTML is a pure function of the file.
        "metric_descriptors": metric_registry.descriptor_manifest_block(),
        # Self-describing pedagogy: classification/edge/evidence explainers +
        # glossary travel inside the file so the (?) affordances render from the
        # manifest, never from hardcoded help in the JS/HTML.
        "pedagogy": pedagogy_registry.pedagogy_manifest_block(),
        "logical_components": [],
        "physical_components": physical,
        "edges": edges,
        "cycles": cycles,
    }

    summary = {
        "files_scanned": len(repo_files),
        "python_files": len(python_files),
        "physical_components": len(physical),
        "excluded_nonsource": excluded_nonsource,
        "edges": len(edges),
        "cycles": len(cycles),
        "parse_failures": parse_failures,
        "metric_tools": metric_tools,
        "unavailable_metrics": unavailable_metrics,
        "regime": regime_hint(len(physical), len(edges)),
    }
    return manifest, summary


def regime_hint(n_components, n_edges):
    """A one-line scaling-regime hint keyed to component/edge counts.

    Tells a learner whether a single synthesis pass can hold the graph (small)
    or whether judgment quality will degrade into "enumerate everything" and
    concerns must be scoped across multiple passes (large). The deterministic
    substrate does not degrade with size — only the judgment layer does — so the
    hint is guidance about the JUDGMENT cost, not the extraction.
    """
    if n_components <= 120 and n_edges <= 400:
        return ("small: one /audit-slice + /synthesize-audit pass suffices; "
                "the whole import graph fits in one model's context.")
    return ("large: scope concerns narrowly and expect MULTIPLE passes "
            "(~1M tokens per the SKILL cost note). The facts scale; the "
            "verdicts do not — trust the structural pass, audit in slices.")


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
    parser.add_argument(
        "--include-nonsource",
        action="store_true",
        default=False,
        dest="include_nonsource",
        help=(
            "Include non-source scaffolding in physical_components: license "
            "files, lockfiles, root-level tooling/dotfiles (.gitignore, "
            ".mcp.json, package.json, ...), harness dirs (.claude/, .github/, "
            "...), data files (.csv/.parquet/.rds/...), and everything under a "
            "top-level docs/ directory. By DEFAULT these are excluded because "
            "they are tracked files that are not architecture and otherwise "
            "drown the signal. (git ls-files already honors .gitignore for "
            "UNTRACKED files; this filter handles TRACKED scaffolding.)"
        ),
    )
    parser.add_argument(
        "--include-audit-output",
        action="store_true",
        default=False,
        dest="include_audit_output",
        help=(
            "Include the skill's own output artifacts (anything under an "
            "'architecture-audit' directory, and the --out directory) in the "
            "physical_components output.  By default they are excluded so the "
            "tool never renders its own prior snapshots as architecture.  Use "
            "only when auditing the skill repository itself."
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
        "kept {phys} physical component(s) ({excl} non-source scaffolding "
        "excluded), {edges} static edges, {cyc} import cycle(s), "
        "{fail} unparseable file(s) skipped.\n".format(
            files=summary["files_scanned"],
            py=summary["python_files"],
            phys=summary["physical_components"],
            excl=summary["excluded_nonsource"],
            edges=summary["edges"],
            cyc=summary["cycles"],
            fail=summary["parse_failures"],
        )
    )
    sys.stderr.write(
        "metric tools used: {tools}\n".format(tools=", ".join(summary["metric_tools"]))
    )

    # LOUD, honest degradation: a missing optional tool must announce itself, not
    # silently produce an empty lens. Both warnings are inspectable by the test
    # agent that installs radon to verify the present-vs-absent paths.
    if summary["unavailable_metrics"]:
        sys.stderr.write(
            "WARNING: radon not installed — cyclomatic & refactor_pressure are "
            "UNAVAILABLE (not zero). A degraded LOC-proxy refactor pressure was "
            "emitted instead and flagged in the manifest "
            "(extractor.unavailable_metrics). Install radon (pip install radon) "
            "for true complexity & refactor-pressure: {keys}\n".format(
                keys=", ".join(summary["unavailable_metrics"])
            )
        )
    if validation.startswith("skipped (jsonschema"):
        sys.stderr.write(
            "WARNING: jsonschema not installed — self-validation was SKIPPED. "
            "Run scripts/validate_components.py (stdlib validator) to gate the "
            "manifest, or pip install jsonschema.\n"
        )

    sys.stderr.write("regime hint — {}\n".format(summary["regime"]))
    sys.stderr.write("schema validation: {}\n".format(validation))
    sys.stderr.write("wrote {}\n".format(out_path))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
