"""mllmct CLI — the thin interface. Outsource execution here; internals are in core/.

Sub-commands: annotate | reconcile | preview-prompt | backfill-metrics | selftest | check-env
Run ``mllmct <cmd> --help`` for arguments. API keys: see ``_resolve_api_keys``.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

from ._version import __version__
from .logging_setup import configure_cli_logging

log = logging.getLogger("mllmct")

# Packaged skill root (dir holding profiles/ + templates/): src/mllmct/cli.py → parents[2].
SKILL_ROOT = Path(__file__).resolve().parents[2]

# Provider routing: a bare model id maps to a vendor by prefix; an id containing '/' is an
# OpenRouter slug. Each provider lists the env vars (first found wins) holding its key.
_PROVIDER_ENV: dict[str, list[str]] = {
    "openrouter": ["OPENROUTER_API_KEY"],
    "gemini": ["GEMINI_API_KEY", "GOOGLE_API_KEY"],
    "openai": ["OPENAI_API_KEY"],
    "anthropic": ["ANTHROPIC_API_KEY"],
}


def _infer_provider(model: str) -> str:
    if "/" in model:
        return "openrouter"
    m = model.lower()
    if m.startswith("gemini"):
        return "gemini"
    if m.startswith(("gpt", "o1", "o3", "o4", "text-")):
        return "openai"
    if m.startswith("claude"):
        return "anthropic"
    return "unknown"


def _parse_api_key_args(items: list[str] | None) -> dict[str, str]:
    out: dict[str, str] = {}
    for it in items or []:
        if "=" not in it:
            raise SystemExit(f"--api-key must be PROVIDER=KEY, got {it!r}")
        prov, key = it.split("=", 1)
        out[prov.strip().lower()] = key.strip()
    return out


def _resolve_api_keys(models: list[str], cli_api_keys: list[str] | None,
                      env_file: str | None) -> dict[str, str]:
    """Resolve provider→key for the model panel. Precedence: --api-key > env (+ --env-file).
    Fail fast with an actionable message if a required key is missing. Keys are never logged.
    """
    if env_file:
        try:
            from dotenv import load_dotenv
            load_dotenv(env_file, override=False)
            log.info("loaded --env-file %s", env_file)
        except Exception as exc:  # noqa: BLE001
            log.warning("could not load --env-file %s: %s", env_file, exc)

    import os

    cli_keys = _parse_api_key_args(cli_api_keys)
    required = sorted({_infer_provider(m) for m in models})
    if "unknown" in required:
        bad = [m for m in models if _infer_provider(m) == "unknown"]
        raise SystemExit(f"cannot route provider for model(s) {bad}. Use an OpenRouter slug "
                         "(e.g. 'openai/gpt-5') or a known vendor prefix (gemini-/gpt-/claude-).")

    api_keys: dict[str, str] = {}
    missing: list[str] = []
    for prov in required:
        if prov in cli_keys:
            api_keys[prov] = cli_keys[prov]
            continue
        env_names = _PROVIDER_ENV.get(prov, [])
        val = next((os.environ[n] for n in env_names if os.environ.get(n)), None)
        if val:
            api_keys[prov] = val
        else:
            missing.append(prov)

    if missing:
        lines = ["Missing API key(s) for provider(s): " + ", ".join(missing), "Provide them by EITHER:",
                 "  1. runtime arg:  " + " ".join(f"--api-key {p}=<KEY>" for p in missing),
                 "  2. environment:  " + "; ".join(
                     f"export {(_PROVIDER_ENV.get(p) or ['<VAR>'])[0]}=<KEY>" for p in missing),
                 "  3. a dotenv:     --env-file /path/to/.env  (with those vars set)"]
        raise SystemExit("\n".join(lines))

    log.info("api keys resolved for provider(s): %s", ", ".join(sorted(api_keys)))
    return api_keys


# ---------------------------------------------------------------------------
# handlers
# ---------------------------------------------------------------------------
def _cmd_annotate(args: argparse.Namespace) -> int:
    from .engine import AnnotationEngine
    from .profile import load_profile

    profile = load_profile(args.profile, skill_root=SKILL_ROOT)
    models = [m.strip() for m in args.models.split(",") if m.strip()]
    api_keys = _resolve_api_keys(models, args.api_key, args.env_file)

    engine = AnnotationEngine(profile, git_cwd=Path.cwd())
    res = engine.run(
        markers_csv=args.markers, evidence_csv=args.evidence, aux_path=args.aux,
        species=args.species, tissue=args.tissue, models=models, api_keys=api_keys,
        lens=args.lens, out_dir=args.out, n_markers=args.n_markers,
        cache_dir=args.cache, consensus_model=args.consensus_model,
    )
    print(str(res.labels_csv))
    return 0


def _cmd_reconcile(args: argparse.Namespace) -> int:
    import pandas as pd

    from .profile import load_profile
    from .reconcile import JoinConfig, reconcile_table

    profile = load_profile(args.profile, skill_root=SKILL_ROOT)
    if not profile.join.get("enabled"):
        raise SystemExit(f"profile {profile.name!r} has join.enabled=false — nothing to reconcile.")
    cfg = JoinConfig.from_profile_join(profile.join, profile.fallback_label)
    df = pd.read_csv(args.axes_csv)
    out = reconcile_table(df, cfg)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "reconciled.csv"
    out.to_csv(out_path, index=False)
    summary = {
        "joint_counts": out["joint_label"].value_counts().to_dict(),
        "method_counts": out["joint_method"].value_counts().to_dict(),
        "n_cells": int(len(out)),
    }
    (out_dir / "reconcile_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    log.info("reconcile methods: %s", summary["method_counts"])
    print(str(out_path))
    return 0


def _cmd_preview_prompt(args: argparse.Namespace) -> int:
    from .core.evidence import load_markers_csv
    from .core.prompt import render_prompt_preview
    from .engine import _build_evidence_block, _fill_template
    from .profile import load_profile

    profile = load_profile(args.profile, skill_root=SKILL_ROOT)
    markers = load_markers_csv(args.markers, args.n_markers)
    evidence: dict[str, str] = {}
    if profile.inject_evidence:
        if not args.evidence:
            raise SystemExit("profile sets inject_evidence:true — pass --evidence")
        from .core.evidence import load_provider
        provider = load_provider(profile.evidence_provider, options=profile.evidence_options,
                                 aux_path=args.aux)
        evidence = provider.assemble(args.evidence)
    template = _fill_template(Path(profile.prompt_template).read_text(),
                              lens=args.lens, profile=profile,
                              evidence_block=_build_evidence_block(evidence))
    print(render_prompt_preview(args.species, args.tissue, template, markers))
    return 0


def _cmd_backfill_metrics(args: argparse.Namespace) -> int:
    import pandas as pd

    from .core.consensus import recompute_lens_metrics

    targets: list[Path] = []
    if args.labels_csv:
        targets.append(Path(args.labels_csv))
    if args.labels_dir:
        targets.extend(sorted(Path(args.labels_dir).rglob("labels.csv")))
    if not targets:
        raise SystemExit("pass --labels-csv or --labels-dir")

    for path in targets:
        df = pd.read_csv(path)
        if "model_annotations" not in df.columns:
            log.warning("%s: no model_annotations column — skipping", path)
            continue
        model_annotations: dict[str, dict[str, str]] = {}
        for _, row in df.iterrows():
            cid = str(row["cluster"])
            per = json.loads(row["model_annotations"]) if pd.notna(row["model_annotations"]) else {}
            for model, lab in per.items():
                model_annotations.setdefault(model, {})[cid] = lab
        cp, ent, maj = recompute_lens_metrics(model_annotations)
        df["py_consensus_proportion"] = df["cluster"].astype(str).map(cp)
        df["py_entropy"] = df["cluster"].astype(str).map(ent)
        df["py_majority_label"] = df["cluster"].astype(str).map(maj)
        df.to_csv(path, index=False)
        log.info("backfilled py_* metrics in %s (%d clusters)", path, len(df))
    return 0


def _cmd_selftest(args: argparse.Namespace) -> int:
    from .selftest import run_selftest
    return run_selftest()


def _cmd_check_env(args: argparse.Namespace) -> int:
    import os

    # 1. versions + monkeypatch seams
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "checks"))
    from smoke_check_versions import main as smoke_main
    rc = smoke_main()

    # 2. API key presence for the given panel (reported, never printed)
    if args.models:
        models = [m.strip() for m in args.models.split(",") if m.strip()]
        if args.env_file:
            try:
                from dotenv import load_dotenv
                load_dotenv(args.env_file, override=False)
            except Exception:  # noqa: BLE001
                pass
        cli_keys = _parse_api_key_args(args.api_key)
        print("check-env: API keys")
        for prov in sorted({_infer_provider(m) for m in models}):
            if prov == "unknown":
                print(f"  [WARN] unroutable model provider in {models}")
                continue
            have = prov in cli_keys or any(os.environ.get(n) for n in _PROVIDER_ENV.get(prov, []))
            print(f"  [{'PASS' if have else 'MISS'}] {prov} key {'found' if have else 'NOT found'}")
    return rc


# ---------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="mllmct", description="Version-locked multi-LLM consensus "
                                "cell-type / cell-state annotation.")
    p.add_argument("--version", action="version", version=f"mllmct {__version__}")
    p.add_argument("--log-file", default=None, help="also write logs to this file")
    p.add_argument("--log-level", default="INFO")
    sub = p.add_subparsers(dest="cmd", required=True)

    a = sub.add_parser("annotate", help="annotate clusters (cell-type or cell-state, per the profile)")
    a.add_argument("--profile", required=True)
    a.add_argument("--markers", required=True, help="CSV: cluster, markers (';'-joined symbols)")
    a.add_argument("--evidence", default=None, help="CSV of per-cluster evidence (cell-state profiles)")
    a.add_argument("--aux", default=None, help="optional provider aux table (e.g. program-decode map)")
    a.add_argument("--species", required=True)
    a.add_argument("--tissue", required=True)
    a.add_argument("--models", required=True, help="comma-separated model ids/slugs")
    a.add_argument("--consensus-model", default=None)
    a.add_argument("--lens", default="lens", help="name for trace dir + provenance")
    a.add_argument("--out", required=True)
    a.add_argument("--n-markers", type=int, default=15)
    a.add_argument("--cache", default=None, help="mllmcelltype response cache dir (default <out>/cache)")
    a.add_argument("--api-key", action="append", help="PROVIDER=KEY (repeatable)")
    a.add_argument("--env-file", default=None)
    a.set_defaults(func=_cmd_annotate)

    r = sub.add_parser("reconcile", help="two-axis semantic join of a per-cell axes table (no LLM)")
    r.add_argument("--profile", required=True)
    r.add_argument("--axes-csv", required=True,
                   help="CSV: axis_a_label, axis_a_conf, axis_b_label, axis_b_conf (per cell)")
    r.add_argument("--out", required=True)
    r.set_defaults(func=_cmd_reconcile)

    pv = sub.add_parser("preview-prompt", help="render the exact prompt (no API call)")
    pv.add_argument("--profile", required=True)
    pv.add_argument("--markers", required=True)
    pv.add_argument("--evidence", default=None)
    pv.add_argument("--aux", default=None)
    pv.add_argument("--species", required=True)
    pv.add_argument("--tissue", required=True)
    pv.add_argument("--lens", default="lens")
    pv.add_argument("--n-markers", type=int, default=15)
    pv.set_defaults(func=_cmd_preview_prompt)

    b = sub.add_parser("backfill-metrics", help="recompute py_* metrics from model_annotations (no API)")
    b.add_argument("--labels-csv", default=None)
    b.add_argument("--labels-dir", default=None)
    b.set_defaults(func=_cmd_backfill_metrics)

    s = sub.add_parser("selftest", help="offline core-logic smoke test (no API)")
    s.set_defaults(func=_cmd_selftest)

    c = sub.add_parser("check-env", help="assert pins + monkeypatch seams + key presence")
    c.add_argument("--models", default=None)
    c.add_argument("--api-key", action="append")
    c.add_argument("--env-file", default=None)
    c.set_defaults(func=_cmd_check_env)

    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    configure_cli_logging(level=getattr(args, "log_level", "INFO"),
                          log_file=getattr(args, "log_file", None))
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
