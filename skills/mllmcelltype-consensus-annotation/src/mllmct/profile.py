"""Profile = the one knob that turns cell-TYPE into cell-STATE (and back).

A profile is a YAML file (versionable, diffable, citable in provenance) that bundles
everything that differs between annotation modes: whether to inject evidence, which prompt
template to use, the vocabulary mode + terms + synonyms, the domain guards, and an optional
two-axis join config for ``reconcile``. Switching modes = passing a different ``--profile``.
No code edits.

``load_profile`` merges a YAML over built-in defaults and validates types/required keys at
load time (fail fast, before any API call), resolving ``prompt_template`` relative to the
skill root.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

_VOCAB_MODES = ("open", "closed")


@dataclass
class Profile:
    name: str
    # evidence
    inject_evidence: bool = False
    evidence_provider: str | None = None
    evidence_options: dict[str, Any] = field(default_factory=dict)
    aux_table_key: str | None = None  # which CLI input feeds the provider's aux_path (e.g. "cnmf")
    # prompt
    prompt_template: Path | None = None  # resolved absolute path
    # harmonization
    vocab_mode: str = "open"
    vocab: list[str] = field(default_factory=list)
    synonyms: dict[str, str] = field(default_factory=dict)
    fallback_label: str = "Ambiguous_LowSignal"
    novel_prefix: str = "Novel:"
    clean_novel: bool = True
    # guards
    domain_guards: dict[str, Any] = field(default_factory=dict)
    # models (CLI --models overrides)
    models: list[str] = field(default_factory=list)
    consensus_model: str | None = None
    # consensus thresholds
    consensus_threshold: float = 0.7
    entropy_threshold: float = 1.0
    max_discussion_rounds: int = 2
    # two-axis join (consumed only by `reconcile`)
    join: dict[str, Any] = field(default_factory=lambda: {"enabled": False})
    # provenance
    source_path: Path | None = None

    @property
    def reject_patterns(self) -> list[str]:
        return list(self.domain_guards.get("reject_patterns", []))

    @property
    def restricted_states(self) -> dict[str, list[str]]:
        return dict(self.domain_guards.get("restricted_states", {}))


_DEFAULTS: dict[str, Any] = {
    "inject_evidence": False,
    "evidence_provider": None,
    "evidence_options": {},
    "aux_table_key": None,
    "vocab_mode": "open",
    "vocab": [],
    "synonyms": {},
    "fallback_label": "Ambiguous_LowSignal",
    "novel_prefix": "Novel:",
    "clean_novel": True,
    "domain_guards": {},
    "models": [],
    "consensus_model": None,
    "consensus_threshold": 0.7,
    "entropy_threshold": 1.0,
    "max_discussion_rounds": 2,
    "join": {"enabled": False},
}


def load_profile(path: str | Path, skill_root: str | Path | None = None) -> Profile:
    """Load and validate a profile YAML.

    Args:
        path: profile YAML file.
        skill_root: base dir for resolving ``prompt_template``. Defaults to the profile
            file's parent's parent (profiles live directly under the skill root).
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"profile not found: {path}")
    raw = yaml.safe_load(path.read_text()) or {}
    if not isinstance(raw, dict):
        raise ValueError(f"profile {path} must be a YAML mapping")
    if "name" not in raw or not str(raw.get("name", "")).strip():
        raise ValueError(f"profile {path} missing required 'name'")

    cfg: dict[str, Any] = dict(_DEFAULTS)
    cfg.update(raw)

    if cfg["vocab_mode"] not in _VOCAB_MODES:
        raise ValueError(f"profile {path}: vocab_mode must be one of {_VOCAB_MODES}, "
                         f"got {cfg['vocab_mode']!r}")
    if cfg["inject_evidence"] and not cfg["evidence_provider"]:
        raise ValueError(f"profile {path}: inject_evidence is true but evidence_provider is null")
    for key, typ in (("vocab", list), ("synonyms", dict), ("evidence_options", dict),
                     ("domain_guards", dict), ("models", list), ("join", dict)):
        if not isinstance(cfg[key], typ):
            raise ValueError(f"profile {path}: '{key}' must be a {typ.__name__}")

    # Resolve prompt_template against, in order: absolute → skill_root → the profile's own
    # directory. This lets a profile anywhere reuse a SHIPPED template ("templates/...",
    # relative to skill_root) OR ship its own template alongside it.
    if skill_root is None:
        skill_root = path.parent.parent  # shipped profiles live in <skill>/profiles/
    skill_root = Path(skill_root)

    tmpl = cfg.get("prompt_template")
    tmpl_path: Path | None = None
    if tmpl:
        raw_tmpl = Path(tmpl)
        if raw_tmpl.is_absolute():
            candidates = [raw_tmpl]
        else:
            candidates = [(skill_root / raw_tmpl), (path.parent / raw_tmpl)]
        tmpl_path = next((c.resolve() for c in candidates if c.exists()), None)
        if tmpl_path is None:
            tried = ", ".join(str(c) for c in candidates)
            raise FileNotFoundError(f"profile {path}: prompt_template {tmpl!r} not found (tried: {tried})")

    return Profile(
        name=str(cfg["name"]),
        inject_evidence=bool(cfg["inject_evidence"]),
        evidence_provider=cfg["evidence_provider"],
        evidence_options=dict(cfg["evidence_options"]),
        aux_table_key=cfg["aux_table_key"],
        prompt_template=tmpl_path,
        vocab_mode=str(cfg["vocab_mode"]),
        vocab=list(cfg["vocab"]),
        synonyms=dict(cfg["synonyms"]),
        fallback_label=str(cfg["fallback_label"]),
        novel_prefix=str(cfg["novel_prefix"]),
        clean_novel=bool(cfg["clean_novel"]),
        domain_guards=dict(cfg["domain_guards"]),
        models=list(cfg["models"]),
        consensus_model=cfg["consensus_model"],
        consensus_threshold=float(cfg["consensus_threshold"]),
        entropy_threshold=float(cfg["entropy_threshold"]),
        max_discussion_rounds=int(cfg["max_discussion_rounds"]),
        join=dict(cfg["join"]),
        source_path=path,
    )
