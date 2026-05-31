"""Offline, dependency-free self-test for ``mllmct selftest`` and ``mllmct check-env``.

Exercises the core invariants WITHOUT pytest, an API, or AnnData, so it runs immediately
after ``uv sync`` (which installs runtime deps only). Returns 0 on success, 1 on any failure.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path
from types import SimpleNamespace


def run_selftest() -> int:
    failures: list[str] = []

    # 1. consensus metrics ---------------------------------------------------
    from .core.consensus import compute_consensus_metrics, recompute_lens_metrics

    lbl, cp, ent = compute_consensus_metrics({"a": "X", "b": "X", "c": "X"})
    if not (lbl == "X" and cp == 1.0 and ent == 0.0):
        failures.append(f"unanimous metrics wrong: {(lbl, cp, ent)}")
    lbl, cp, ent = compute_consensus_metrics({"a": "X", "b": "Y"})
    if not (cp == 0.5 and abs(ent - 1.0) < 1e-9):
        failures.append(f"split metrics wrong: {(lbl, cp, ent)}")
    cps, ents, majs = recompute_lens_metrics({"m1": {"0": "X", "1": "Y"}, "m2": {"0": "X", "1": "Z"}})
    if cps.get("0") != 1.0 or cps.get("1") != 0.5:
        failures.append(f"recompute_lens_metrics wrong: {cps}")

    # 2. harmonize -----------------------------------------------------------
    from .core.harmonize import harmonize_label, is_novel_label

    if harmonize_label("interferon high", ["ISG_High"], {"interferon": "ISG_High"}) != "ISG_High":
        failures.append("synonym harmonize failed")
    nov = harmonize_label("weird new state", ["A"], {}, mode="open")
    if not is_novel_label(nov):
        failures.append(f"open-vocab novel failed: {nov}")
    if harmonize_label("weird", ["A"], {}, mode="closed", fallback_label="Amb") != "Amb":
        failures.append("closed-vocab fallback failed")
    if harmonize_label("CD8 T cells", [], {}, novel_prefix="", clean_novel=False) != "CD8 T cells":
        failures.append("celltype passthrough failed")

    # 3. guard ---------------------------------------------------------------
    from .core.guards import make_reject_guard

    g = make_reject_guard([r"\bth1\b"])
    if not g("th1 cells") or g("dna damage response cells"):
        failures.append("reject guard wrong")

    # 4. determinism config rewrite -----------------------------------------
    from .core.capture import force_deterministic_config

    d = force_deterministic_config({"temperature": 0.7, "max_output_tokens": 4096})
    if d.get("temperature") != 0.0 or d.get("max_output_tokens") != 4096 or d.get("seed") != 0:
        failures.append(f"force_deterministic_config(dict) wrong: {d}")

    # 5. token capture → trace non-clobber → cost ---------------------------
    from .core.capture import GeminiTokenCapture
    from .core.cost import report_cost, write_cost_summary
    from .core.trace import TraceWriter

    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp) / "trace"
        cap = GeminiTokenCapture()
        cap._accumulate("gemini-2.5-flash",
                        SimpleNamespace(prompt_token_count=4000, candidates_token_count=200,
                                        total_token_count=4200))
        tw = TraceWriter("lensA", base_dir=base / "lensA")
        tw.write_tokens(cap.usage)
        rec = json.loads((base / "lensA" / "tokens.json").read_text())
        if rec.get("status") == "unavailable" or rec.get("totals", {}).get("total") != 4200:
            failures.append(f"write_tokens lost real capture: {rec}")
        cost = report_cost(cap.usage, ["gemini-2.5-flash"])
        if cost.get("fully_cached") or cost.get("total_usd", 0) <= 0:
            failures.append(f"report_cost on live capture wrong: {cost}")
        # fully-cached re-run must NOT clobber
        tw.write_tokens({})
        rec2 = json.loads((base / "lensA" / "tokens.json").read_text())
        if rec2.get("last_known", {}).get("totals", {}).get("total") != 4200:
            failures.append(f"cached re-run clobbered prior tokens: {rec2}")
        summary = write_cost_summary(["lensA"], ["gemini-2.5-flash"], base_dir=base)
        if summary["grand_total_tokens"]["total"] != 4200 or summary["grand_total_usd"] <= 0:
            failures.append(f"cost_summary wrong: {summary['grand_total_tokens']}")

    # 6. reconcile join branches --------------------------------------------
    from .reconcile import JoinConfig, canonical_two_axis, join_cell

    cfg = JoinConfig(distinct_labels=["DistinctA"], order=["DistinctA", "B"],
                     uninformative_labels=["Amb"], conf_floor=0.5)
    if join_cell("Amb", 1.0, False, "Amb", 1.0, False, cfg)[1] != "both_uninformative":
        failures.append("join both_uninformative failed")
    if join_cell("DistinctA", 1.0, False, "B", 1.0, False, cfg)[1] != "two_axis":
        failures.append("join two_axis failed")
    if canonical_two_axis("B", "DistinctA", cfg) != canonical_two_axis("DistinctA", "B", cfg):
        failures.append("canonical_two_axis not order-invariant")

    if failures:
        print("SELFTEST FAILED:")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("SELFTEST PASSED: consensus, harmonize, guard, determinism, token-capture "
          "non-clobber, cost, and reconcile invariants all hold.")
    return 0
