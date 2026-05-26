#!/usr/bin/env python3
"""
rollup_logical_edges.py — Deterministic rollup of physical static import edges
to logical-component edges in an architecture-treemap components.json manifest.

Usage:
    python3 rollup_logical_edges.py <components.json> [--out <path>]

The script is IDEMPOTENT: running it twice produces a byte-identical file.
Any existing derived edges (marked `"derived": true`) are removed before
re-deriving, so repeated runs accumulate nothing.

Algorithm:
  1. Build physical_id -> path from physical_components[].
  2. Build path -> logical_owner_id from physical_components[].logical_owners
     (the canonical ownership map; each physical file knows its logical owners).
  3. For each static edge whose from/to are physical ids:
     map both endpoints to their logical owners.
     If both have owners and the owners differ, record the logical pair.
  4. Collapse all physical imports between the same ordered logical pair into
     ONE logical edge, tagging it derived=true with an evidence summary.
  5. Write back in place (or to --out if specified).

Exit codes:
    0 — success
    1 — error (malformed JSON, missing keys, etc.)
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path


def _build_physical_id_to_path(data: dict) -> dict[str, str]:
    """Map physical component id -> file path."""
    return {
        pc["id"]: pc["path"]
        for pc in data.get("physical_components", [])
        if "id" in pc and "path" in pc
    }


def _build_path_to_logical_owners(data: dict) -> dict[str, list[str]]:
    """
    Map file path -> list of logical owner ids.

    Reads from physical_components[].logical_owners (the extractor stamps this).
    Falls back to scanning logical_components[].physical_files if a physical
    component has no logical_owners entry (hand-authored manifests).
    """
    path_to_owners: dict[str, list[str]] = defaultdict(list)

    # Primary source: physical_components[].logical_owners
    # Append each owner at most once per path: a duplicate owner id (whether
    # repeated within one logical_owners list or contributed by two physical
    # components sharing a path) would otherwise multiply a single physical
    # import into N appends of the same logical pair, inflating the derived
    # edge's evidence count and example list.
    for pc in data.get("physical_components", []):
        path = pc.get("path")
        owners = pc.get("logical_owners") or []
        if not path:
            continue
        for owner in owners:
            if owner not in path_to_owners.get(path, []):
                path_to_owners[path].append(owner)

    # Fallback: logical_components[].physical_files (catches hand-authored manifests
    # where the physical component list has no logical_owners annotations).
    logical_ids = {lc["id"] for lc in data.get("logical_components", []) if "id" in lc}
    for lc in data.get("logical_components", []):
        lc_id = lc.get("id")
        if not lc_id:
            continue
        for pf in lc.get("physical_files", []):
            path = pf.get("path")
            if path and lc_id not in path_to_owners.get(path, []):
                path_to_owners[path].append(lc_id)

    return dict(path_to_owners)


def _format_evidence(pair_imports: list[tuple[str, str, str, str]], total: int) -> str:
    """
    Format the evidence string for a derived logical edge.
    Shows up to 3 concrete from_path:to_path examples plus the total count.
    """
    cap = 3
    examples = pair_imports[:cap]
    parts = [f"{fp} -> {tp}" for fp, tp, _fi, _ti in examples]
    summary = "; ".join(parts)
    if total > cap:
        summary += f"; … ({total} total)"
    else:
        summary += f" ({total} total)"
    return f"rolled up from {total} static import(s): {summary}"


def derive_logical_edges(data: dict) -> tuple[list[dict], int]:
    """
    Derive logical-level edges from the physical static import graph.

    Returns:
        (derived_edges, count_physical_used)
        derived_edges: list of new edge dicts tagged derived=true
        count_physical_used: number of physical static edges that contributed
    """
    physical_id_to_path = _build_physical_id_to_path(data)
    path_to_owners = _build_path_to_logical_owners(data)
    physical_ids = set(physical_id_to_path.keys())

    # Bucket: (from_logical, to_logical) -> list of (from_path, to_path, from_id, to_id)
    pair_to_imports: dict[tuple[str, str], list[tuple[str, str, str, str]]] = defaultdict(list)
    physical_used = 0

    for edge in data.get("edges", []):
        # Only process static edges keyed by physical ids
        if edge.get("evidence_class") != "static":
            continue
        from_id = edge.get("from", "")
        to_id = edge.get("to", "")
        if from_id not in physical_ids or to_id not in physical_ids:
            continue

        from_path = physical_id_to_path[from_id]
        to_path = physical_id_to_path[to_id]
        from_owners = path_to_owners.get(from_path, [])
        to_owners = path_to_owners.get(to_path, [])

        if not from_owners or not to_owners:
            continue

        physical_used += 1

        # Emit one pair per distinct (from_owner, to_owner) combination
        # (rare case: a file owned by multiple logical components)
        for fl in from_owners:
            for tl in to_owners:
                if fl != tl:
                    pair_to_imports[(fl, tl)].append((from_path, to_path, from_id, to_id))

    derived: list[dict] = []
    for (from_logical, to_logical), imports in sorted(pair_to_imports.items()):
        total = len(imports)
        evidence = _format_evidence(imports, total)
        derived.append({
            "from": from_logical,
            "to": to_logical,
            "type": "direct-call",
            "evidence_class": "static",
            "epistemic_source": "measured",
            "evidence": evidence,
            "derived": True,
        })

    return derived, physical_used


def rollup(manifest_path: Path, out_path: Path) -> None:
    """Run the rollup on the manifest and write the result."""
    try:
        with manifest_path.open("r", encoding="utf-8") as fh:
            data = json.load(fh)
    except json.JSONDecodeError as exc:
        print(f"ERROR: {manifest_path}: JSON parse error: {exc}", file=sys.stderr)
        sys.exit(1)
    except OSError as exc:
        print(f"ERROR: cannot read {manifest_path}: {exc}", file=sys.stderr)
        sys.exit(1)

    if "edges" not in data:
        print(f"ERROR: {manifest_path}: missing 'edges' key", file=sys.stderr)
        sys.exit(1)

    # IDEMPOTENCY: remove any previously derived edges before re-deriving
    data["edges"] = [e for e in data["edges"] if not e.get("derived")]

    # Derive logical edges from the physical static graph
    derived_edges, physical_used = derive_logical_edges(data)

    # Append derived edges after the existing edges
    data["edges"].extend(derived_edges)

    # Write back (in place or to --out)
    try:
        with out_path.open("w", encoding="utf-8") as fh:
            json.dump(data, fh, indent=2, ensure_ascii=False)
            fh.write("\n")
    except OSError as exc:
        print(f"ERROR: cannot write {out_path}: {exc}", file=sys.stderr)
        sys.exit(1)

    n = len(derived_edges)
    m = physical_used
    print(
        f"rollup_logical_edges: derived {n} logical edge(s) from {m} physical static import(s)",
        file=sys.stderr,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Derive logical-component edges from physical static import edges "
            "in an architecture-treemap components.json manifest."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("components_json", help="Path to components.json")
    parser.add_argument(
        "--out",
        metavar="PATH",
        help="Output path (default: overwrite the input file in place).",
    )
    args = parser.parse_args()

    manifest_path = Path(args.components_json)
    if not manifest_path.exists():
        print(f"ERROR: file not found: {manifest_path}", file=sys.stderr)
        sys.exit(1)

    out_path = Path(args.out) if args.out else manifest_path

    rollup(manifest_path, out_path)


if __name__ == "__main__":
    main()
