#!/usr/bin/env python3
"""
validate_components.py — JSON Schema + referential-integrity validator
for architecture-treemap components.json manifests.

Usage:
    python validate_components.py <components.json> [--strict]

Exit codes:
    0  — valid (warnings may be printed)
    1  — schema failure (always fatal) or referential-integrity failure (--strict only)
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Schema file location (resolved relative to this script)
# ---------------------------------------------------------------------------
_SCHEMA_PATH = Path(__file__).parent.parent / "components.schema.json"

# ---------------------------------------------------------------------------
# Lightweight JSON-schema-like validator
# We implement only the subset actually used by components.schema.json rather
# than pulling in the full jsonschema package.  The validation logic mirrors
# the schema structure exactly.
# ---------------------------------------------------------------------------

VALID_CLASSIFICATIONS = {"core", "seam", "removable"}
VALID_EDGE_TYPES       = {"direct-call", "shared-state", "background-knowledge"}
VALID_EVIDENCE_CLASSES = {"static", "audit-asserted"}
VALID_KINDS            = {"module", "config", "test", "doc", "asset", "other"}
VALID_SEVERITIES       = {"info", "warn", "error"}
VALID_OWNERSHIPS       = {"primary", "partial"}
ID_PATTERN             = re.compile(r"^[a-z0-9-]+$")
LINE_RANGE_PATTERN     = re.compile(r"^[0-9]+(-[0-9]+)?$")


class ValidationError(Exception):
    """Represents one validation failure."""
    def __init__(self, path: str, message: str):
        self.path = path
        self.message = message
        super().__init__(f"{path}: {message}")


def _require_type(path: str, value: Any, expected_type, type_label: str):
    if not isinstance(value, expected_type):
        raise ValidationError(path, f"expected {type_label}, got {type(value).__name__}")


def _require_key(path: str, obj: dict, key: str):
    if key not in obj:
        raise ValidationError(path, f"required key '{key}' is missing")


def _validate_physical_file_ref(path: str, ref: Any):
    _require_type(path, ref, dict, "object")
    _require_key(path, ref, "path")
    _require_type(f"{path}.path", ref["path"], str, "string")
    if "line_range" in ref:
        _require_type(f"{path}.line_range", ref["line_range"], str, "string")
        if not LINE_RANGE_PATTERN.match(ref["line_range"]):
            raise ValidationError(f"{path}.line_range", f"does not match pattern ^[0-9]+(-[0-9]+)?$")
    if "ownership" in ref:
        if ref["ownership"] not in VALID_OWNERSHIPS:
            raise ValidationError(f"{path}.ownership", f"must be one of {VALID_OWNERSHIPS}")


def _validate_metrics(path: str, m: Any):
    """Validate the OPEN numeric metric map.

    The set of keys is not fixed (the metric registry is the single source of
    truth for which metrics exist). Every value must be a non-negative number.
    The historically-integer fields are still type-checked as integers for the
    extra rigour the substrate guarantees; all other keys are open numerics.
    """
    _require_type(path, m, dict, "object")
    int_fields = {"loc", "fan_in", "fan_out", "churn_90d"}
    for field, value in m.items():
        if field in int_fields:
            # bool is a subclass of int; reject it explicitly.
            if isinstance(value, bool) or not isinstance(value, int):
                raise ValidationError(f"{path}.{field}", "expected integer")
            if value < 0:
                raise ValidationError(f"{path}.{field}", "must be >= 0")
        else:
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise ValidationError(f"{path}.{field}", "expected number")
            if value < 0:
                raise ValidationError(f"{path}.{field}", "must be >= 0")


def _validate_logical_component(path: str, lc: Any):
    _require_type(path, lc, dict, "object")
    for key in ["id", "name", "description", "classification", "physical_files", "size_estimate_loc"]:
        _require_key(path, lc, key)

    _require_type(f"{path}.id", lc["id"], str, "string")
    if not ID_PATTERN.match(lc["id"]):
        raise ValidationError(f"{path}.id", "must match ^[a-z0-9-]+$")

    _require_type(f"{path}.name", lc["name"], str, "string")
    _require_type(f"{path}.description", lc["description"], str, "string")

    if lc["classification"] not in VALID_CLASSIFICATIONS:
        raise ValidationError(f"{path}.classification", f"must be one of {VALID_CLASSIFICATIONS}")

    _require_type(f"{path}.physical_files", lc["physical_files"], list, "array")
    for i, ref in enumerate(lc["physical_files"]):
        _validate_physical_file_ref(f"{path}.physical_files[{i}]", ref)

    _require_type(f"{path}.size_estimate_loc", lc["size_estimate_loc"], int, "integer")
    if lc["size_estimate_loc"] < 0:
        raise ValidationError(f"{path}.size_estimate_loc", "must be >= 0")

    if "smoke_findings" in lc:
        _require_type(f"{path}.smoke_findings", lc["smoke_findings"], list, "array")
        for i, s in enumerate(lc["smoke_findings"]):
            _require_type(f"{path}.smoke_findings[{i}]", s, str, "string")

    if "audit_slices" in lc:
        _require_type(f"{path}.audit_slices", lc["audit_slices"], list, "array")
        for i, s in enumerate(lc["audit_slices"]):
            _require_type(f"{path}.audit_slices[{i}]", s, str, "string")

    if "metrics" in lc:
        _validate_metrics(f"{path}.metrics", lc["metrics"])


def _validate_physical_component(path: str, pc: Any):
    _require_type(path, pc, dict, "object")
    for key in ["id", "path", "size_loc", "kind"]:
        _require_key(path, pc, key)

    _require_type(f"{path}.id", pc["id"], str, "string")
    if not ID_PATTERN.match(pc["id"]):
        raise ValidationError(f"{path}.id", "must match ^[a-z0-9-]+$")

    _require_type(f"{path}.path", pc["path"], str, "string")

    _require_type(f"{path}.size_loc", pc["size_loc"], int, "integer")
    if pc["size_loc"] < 0:
        raise ValidationError(f"{path}.size_loc", "must be >= 0")

    if pc["kind"] not in VALID_KINDS:
        raise ValidationError(f"{path}.kind", f"must be one of {VALID_KINDS}")

    if "logical_owners" in pc:
        _require_type(f"{path}.logical_owners", pc["logical_owners"], list, "array")
        for i, o in enumerate(pc["logical_owners"]):
            _require_type(f"{path}.logical_owners[{i}]", o, str, "string")

    if "metrics" in pc:
        _validate_metrics(f"{path}.metrics", pc["metrics"])


def _validate_edge(path: str, e: Any):
    _require_type(path, e, dict, "object")
    for key in ["from", "to", "type"]:
        _require_key(path, e, key)

    _require_type(f"{path}.from", e["from"], str, "string")
    _require_type(f"{path}.to",   e["to"],   str, "string")

    if e["type"] not in VALID_EDGE_TYPES:
        raise ValidationError(f"{path}.type", f"must be one of {VALID_EDGE_TYPES}")

    if "evidence_class" in e and e["evidence_class"] is not None:
        if e["evidence_class"] not in VALID_EVIDENCE_CLASSES:
            raise ValidationError(f"{path}.evidence_class", f"must be one of {VALID_EVIDENCE_CLASSES}")

    if "evidence" in e and e["evidence"] is not None:
        _require_type(f"{path}.evidence", e["evidence"], str, "string")

    if "smell" in e and e["smell"] is not None:
        _require_type(f"{path}.smell", e["smell"], str, "string")

    if "smoke_finding" in e and e["smoke_finding"] is not None:
        _require_type(f"{path}.smoke_finding", e["smoke_finding"], str, "string")


def _validate_prune_candidate(path: str, p: Any):
    _require_type(path, p, dict, "object")
    for key in ["logical_component", "bounded", "loc_to_remove"]:
        _require_key(path, p, key)

    _require_type(f"{path}.logical_component", p["logical_component"], str, "string")
    _require_type(f"{path}.bounded", p["bounded"], bool, "boolean")
    _require_type(f"{path}.loc_to_remove", p["loc_to_remove"], int, "integer")
    if p["loc_to_remove"] < 0:
        raise ValidationError(f"{path}.loc_to_remove", "must be >= 0")

    if "unblocks_simplification_of" in p:
        _require_type(f"{path}.unblocks_simplification_of", p["unblocks_simplification_of"], list, "array")
        for i, s in enumerate(p["unblocks_simplification_of"]):
            _require_type(f"{path}.unblocks_simplification_of[{i}]", s, str, "string")


def _validate_core_boundary(path: str, cb: Any):
    _require_type(path, cb, dict, "object")
    for key in ["core", "seam", "removable"]:
        if key in cb:
            _require_type(f"{path}.{key}", cb[key], list, "array")
            for i, s in enumerate(cb[key]):
                _require_type(f"{path}.{key}[{i}]", s, str, "string")
    if "rationale" in cb:
        _require_type(f"{path}.rationale", cb["rationale"], str, "string")


def _validate_finding_index(path: str, fi: Any):
    _require_type(path, fi, dict, "object")
    for tag, record in fi.items():
        _require_type(f"{path}.{tag}", record, dict, "object")
        if "title" in record:
            _require_type(f"{path}.{tag}.title", record["title"], str, "string")
        if "severity" in record:
            if record["severity"] not in VALID_SEVERITIES:
                raise ValidationError(f"{path}.{tag}.severity", f"must be one of {VALID_SEVERITIES}")
        if "url" in record:
            _require_type(f"{path}.{tag}.url", record["url"], str, "string")


def validate_schema(data: dict) -> list[str]:
    """
    Validate the structural schema of a components manifest.
    Returns a list of error strings.  Empty list means valid.
    """
    errors: list[str] = []

    # Top-level must be object
    if not isinstance(data, dict):
        return ["root: expected object, got " + type(data).__name__]

    # Required top-level keys
    for key in ["schema_version", "logical_components", "physical_components", "edges"]:
        if key not in data:
            errors.append(f"root: required key '{key}' is missing")

    if errors:
        return errors  # cannot continue without these

    # schema_version const
    if data.get("schema_version") != "1.0":
        errors.append(f"schema_version: must be '1.0', got {data.get('schema_version')!r}")

    # Optional string scalars
    for key in ["audit_date", "project", "frame", "git_sha", "snapshot_id"]:
        if key in data and not isinstance(data[key], str):
            errors.append(f"{key}: expected string")

    # extractor block
    if "extractor" in data:
        ext = data["extractor"]
        if not isinstance(ext, dict):
            errors.append("extractor: expected object")

    # metric_descriptors block (self-describing metric model; optional)
    if "metric_descriptors" in data:
        md = data["metric_descriptors"]
        if not isinstance(md, dict):
            errors.append("metric_descriptors: expected object")
        else:
            for key, desc in md.items():
                if not isinstance(desc, dict):
                    errors.append(f"metric_descriptors.{key}: expected object")

    # logical_components
    if isinstance(data.get("logical_components"), list):
        for i, lc in enumerate(data["logical_components"]):
            try:
                _validate_logical_component(f"logical_components[{i}]", lc)
            except ValidationError as exc:
                errors.append(str(exc))
    else:
        errors.append("logical_components: expected array")

    # physical_components
    if isinstance(data.get("physical_components"), list):
        for i, pc in enumerate(data["physical_components"]):
            try:
                _validate_physical_component(f"physical_components[{i}]", pc)
            except ValidationError as exc:
                errors.append(str(exc))
    else:
        errors.append("physical_components: expected array")

    # edges
    if isinstance(data.get("edges"), list):
        for i, e in enumerate(data["edges"]):
            try:
                _validate_edge(f"edges[{i}]", e)
            except ValidationError as exc:
                errors.append(str(exc))
    else:
        errors.append("edges: expected array")

    # prune_candidates
    if "prune_candidates" in data:
        if isinstance(data["prune_candidates"], list):
            for i, p in enumerate(data["prune_candidates"]):
                try:
                    _validate_prune_candidate(f"prune_candidates[{i}]", p)
                except ValidationError as exc:
                    errors.append(str(exc))
        else:
            errors.append("prune_candidates: expected array")

    # core_boundary
    if "core_boundary" in data:
        try:
            _validate_core_boundary("core_boundary", data["core_boundary"])
        except ValidationError as exc:
            errors.append(str(exc))

    # finding_index
    if "finding_index" in data:
        try:
            _validate_finding_index("finding_index", data["finding_index"])
        except ValidationError as exc:
            errors.append(str(exc))

    return errors


def validate_referential_integrity(data: dict) -> list[str]:
    """
    Check that edge endpoints and physical_owner references exist within the
    manifest.  The check works across BOTH id-spaces:
      - Logical ids  (logical_components[].id)
      - Physical ids (physical_components[].id)
    An edge endpoint is valid if it appears in either space.
    Logical-owner entries in physical_components[] must resolve to logical ids.

    Returns a list of warning strings (may be empty).
    """
    warnings: list[str] = []

    logical_ids  = {lc["id"] for lc in data.get("logical_components", []) if "id" in lc}
    physical_ids = {pc["id"] for pc in data.get("physical_components", []) if "id" in pc}
    all_ids = logical_ids | physical_ids

    for i, edge in enumerate(data.get("edges", [])):
        frm = edge.get("from", "")
        to  = edge.get("to",   "")
        if frm and frm not in all_ids:
            warnings.append(f"edges[{i}].from: '{frm}' not found in logical or physical ids")
        if to and to not in all_ids:
            warnings.append(f"edges[{i}].to: '{to}' not found in logical or physical ids")

    for i, pc in enumerate(data.get("physical_components", [])):
        for j, owner in enumerate(pc.get("logical_owners") or []):
            if owner not in logical_ids:
                warnings.append(
                    f"physical_components[{i}].logical_owners[{j}]: "
                    f"'{owner}' not found in logical_components ids"
                )

    return warnings


def run_validation(components_path: Path, strict: bool) -> tuple[list[str], list[str]]:
    """
    Load and validate a components.json file.
    Returns (schema_errors, ref_warnings).
    Raises SystemExit on JSON parse failure.
    """
    try:
        with components_path.open("r", encoding="utf-8") as fh:
            data = json.load(fh)
    except json.JSONDecodeError as exc:
        print(f"ERROR: {components_path}: JSON parse error: {exc}", file=sys.stderr)
        sys.exit(1)

    schema_errors = validate_schema(data)
    ref_warnings  = validate_referential_integrity(data) if not schema_errors else []

    return schema_errors, ref_warnings


def main():
    parser = argparse.ArgumentParser(
        description="Validate a components.json against the architecture-treemap schema.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("components_json", help="Path to components.json")
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Treat referential-integrity mismatches as errors (non-zero exit).",
    )
    args = parser.parse_args()

    components_path = Path(args.components_json)
    if not components_path.exists():
        print(f"ERROR: file not found: {components_path}", file=sys.stderr)
        sys.exit(1)

    schema_errors, ref_warnings = run_validation(components_path, strict=args.strict)

    if schema_errors:
        print(f"SCHEMA ERRORS in {components_path}:")
        for err in schema_errors:
            print(f"  ERROR: {err}")

    if ref_warnings:
        label = "ERROR" if args.strict else "WARNING"
        print(f"\nREFERENTIAL INTEGRITY in {components_path}:")
        for w in ref_warnings:
            print(f"  {label}: {w}")

    if not schema_errors and not ref_warnings:
        print(f"OK: {components_path} is valid.")

    if schema_errors:
        sys.exit(1)

    if args.strict and ref_warnings:
        sys.exit(1)

    sys.exit(0)


if __name__ == "__main__":
    main()
