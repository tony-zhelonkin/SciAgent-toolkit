#!/usr/bin/env python3
"""smoke_check_versions.py — OFFLINE gate: exact pins + monkeypatch seams intact.

mllmct wraps VERSION-SENSITIVE internals of four packages (attribute names + call
signatures, not public APIs). A loose dependency bump can silently no-op a wrapper —
the run would then succeed but lose determinism, token capture, the custom prompt, or
DEBUG capture, with NO error. This script is the tripwire: run it after any
``uv lock`` and in ``mllmct check-env``. It asserts

  1. the four exact versions in mllmct._version.PINNED, and
  2. that every structural seam the monkeypatches bind to still EXISTS.

No network, no API key, no AnnData. Exit 0 = safe to use; exit 1 = a seam moved,
do NOT trust the wrappers until references/monkeypatch-internals.md is reconciled.

Run:  uv run python checks/smoke_check_versions.py
"""

from __future__ import annotations

import importlib
import importlib.metadata as ilm
import sys

# Resolve mllmct._version whether run as a script or a module.
try:
    from mllmct._version import PINNED
except ModuleNotFoundError:  # run directly without install
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from mllmct._version import PINNED


def _check(label: str, ok: bool, detail: str = "") -> bool:
    mark = "PASS" if ok else "FAIL"
    print(f"  [{mark}] {label}" + (f" — {detail}" if detail else ""))
    return ok


def main() -> int:
    print("smoke_check_versions: exact pins")
    failures = 0

    # 1. exact versions ------------------------------------------------------
    for pkg, want in PINNED.items():
        try:
            got = ilm.version(pkg)
        except ilm.PackageNotFoundError:
            failures += not _check(f"{pkg} installed", False, "not found")
            continue
        failures += not _check(f"{pkg}=={want}", got == want, f"got {got}")

    print("smoke_check_versions: monkeypatch seams")

    # 2a. native prompt_template seam (threaded, not monkeypatched since 2.0.7) ---
    #     Seam (i) is retired: the engine passes prompt_template= to
    #     interactive_consensus_annotation instead of swapping DEFAULT_PROMPT_TEMPLATE.
    #     Assert the native parameter exists on both create_prompt (preview) and the
    #     consensus entrypoint (the real call) — if either drops it, custom prompts
    #     silently revert to the stock template.
    try:
        import inspect

        prompts = importlib.import_module("mllmcelltype.prompts")
        create_prompt = getattr(prompts, "create_prompt", None)
        failures += not _check("mllmcelltype.prompts.create_prompt callable",
                               callable(create_prompt))
        if callable(create_prompt):
            failures += not _check(
                "create_prompt accepts 'prompt_template'",
                "prompt_template" in inspect.signature(create_prompt).parameters)
        mc0 = importlib.import_module("mllmcelltype")
        ica = getattr(mc0, "interactive_consensus_annotation", None)
        if callable(ica):
            failures += not _check(
                "interactive_consensus_annotation accepts 'prompt_template'",
                "prompt_template" in inspect.signature(ica).parameters)
    except Exception as exc:  # noqa: BLE001
        failures += not _check("import mllmcelltype.prompts", False, str(exc))

    # 2b. consensus entrypoint + logger seam --------------------------------
    try:
        mc = importlib.import_module("mllmcelltype")
        failures += not _check("mllmcelltype.interactive_consensus_annotation",
                               callable(getattr(mc, "interactive_consensus_annotation", None)))
        logger_mod = importlib.import_module("mllmcelltype.logger")
        failures += not _check("mllmcelltype.logger.setup_logging callable",
                               callable(getattr(logger_mod, "setup_logging", None)))
    except Exception as exc:  # noqa: BLE001
        failures += not _check("import mllmcelltype / logger", False, str(exc))

    # 2c. OpenRouter HTTP seam (requests.post patch) ------------------------
    try:
        import requests as _req

        or_mod = importlib.import_module("mllmcelltype.providers.openrouter")
        same = getattr(or_mod, "requests", None) is _req
        failures += not _check("mllmcelltype.providers.openrouter.requests IS requests", same)
    except Exception as exc:  # noqa: BLE001
        failures += not _check("import mllmcelltype.providers.openrouter", False, str(exc))

    # 2d. google-genai seams (generate_content patch + config/usage fields) -
    try:
        genai_models = importlib.import_module("google.genai.models")
        failures += not _check("google.genai.models.Models.generate_content",
                               callable(getattr(genai_models.Models, "generate_content", None)))
        genai_types = importlib.import_module("google.genai.types")
        cfg = getattr(genai_types, "GenerateContentConfig", None)
        cfg_fields = set(getattr(cfg, "model_fields", {})) if cfg is not None else set()
        failures += not _check("GenerateContentConfig has 'temperature'", "temperature" in cfg_fields)
        failures += not _check("GenerateContentConfig has 'seed'", "seed" in cfg_fields)
        um = getattr(genai_types, "GenerateContentResponseUsageMetadata", None)
        um_fields = set(getattr(um, "model_fields", {})) if um is not None else set()
        for f in ("prompt_token_count", "candidates_token_count", "total_token_count"):
            failures += not _check(f"usage_metadata has '{f}'", f in um_fields)
    except Exception as exc:  # noqa: BLE001
        failures += not _check("import google.genai", False, str(exc))

    print()
    if failures:
        print(f"smoke_check_versions: {failures} FAILURE(S) — wrappers may be unsafe. "
              "Reconcile references/monkeypatch-internals.md before use.")
        return 1
    print("smoke_check_versions: ALL PASS — pins + seams intact.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
