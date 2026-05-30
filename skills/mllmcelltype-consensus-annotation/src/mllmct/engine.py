"""AnnotationEngine — orchestration only (the generalized, de-projectified 15a).

One code path for cell-TYPE and cell-STATE; the active Profile decides the rest. Steps:

    load markers → (if inject_evidence) assemble evidence → fill + install prompt template
    → wrap call in token-capture + DEBUG-capture → interactive_consensus_annotation
    → recompute consensus metrics in Python (#1) → harmonize labels (#8/#10)
    → validate (gates the write, #9) → write labels.csv + trace tree (#6) + cost_summary (#7)

All transformation logic lives in ``mllmct.core`` / the EvidenceProvider; this module only
sequences it and owns the validation gate. Nothing here names a biology.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from mllmcelltype import interactive_consensus_annotation

from .core import capture as _capture
from .core import cost as _cost
from .core import debug_capture as _debug
from .core import prompt as _prompt
from .core.consensus import recompute_lens_metrics
from .core.evidence import load_markers_csv, load_provider
from .core.guards import make_reject_guard
from .core.harmonize import harmonize_label, is_novel_label
from .core.trace import TraceWriter, git_head
from .profile import Profile

log = logging.getLogger(__name__)

_OPEN_CLAUSE = ("If a cluster clearly does not fit any vocabulary term, you MAY propose a "
                "short novel name instead of forcing a fit.")
_CLOSED_CLAUSE = "Do NOT invent labels outside this vocabulary."


@dataclass
class AnnotationResult:
    lens: str
    labels: dict[str, str]            # cluster -> harmonized label
    out_dir: Path
    labels_csv: Path
    validation: dict
    status: str                       # PASS / FAIL


def _numeric_cluster_key(cid: str) -> tuple[int, str]:
    try:
        return (int(cid), "")
    except ValueError:
        return (10**9, cid)


def _build_evidence_block(evidence: dict[str, str]) -> str:
    return "\n".join(
        f"  cluster {cid}: {evidence.get(cid, '')}"
        for cid in sorted(evidence, key=_numeric_cluster_key)
    )


def _fill_template(template: str, *, lens: str, profile: Profile, evidence_block: str) -> str:
    """Fill the engine's ``<<...>>`` slots; leave ``{species}/{tissue}/{markers}`` for the library."""
    mode_clause = _OPEN_CLAUSE if profile.vocab_mode == "open" else _CLOSED_CLAUSE
    guardrails = str(profile.domain_guards.get("guardrails_text", "") or "")
    repl = {
        "<<lens>>": lens,
        "<<vocab>>": ", ".join(profile.vocab),
        "<<fallback_label>>": profile.fallback_label,
        "<<mode_clause>>": mode_clause,
        "<<guardrails>>": guardrails,
        "<<evidence_block>>": evidence_block or "(no per-cluster evidence)",
    }
    for k, v in repl.items():
        template = template.replace(k, v)
    return template


class AnnotationEngine:
    def __init__(self, profile: Profile, git_cwd: str | Path | None = None) -> None:
        self.profile = profile
        self.git_cwd = git_cwd

    # -- evidence + prompt ------------------------------------------------
    def _assemble_evidence(self, evidence_csv: str | None, aux_path: str | None) -> dict[str, str]:
        p = self.profile
        if not p.inject_evidence:
            return {}
        if not evidence_csv:
            raise ValueError("profile sets inject_evidence:true but no --evidence table was given")
        provider = load_provider(p.evidence_provider, options=p.evidence_options, aux_path=aux_path)
        evidence = provider.assemble(evidence_csv)
        log.info("evidence assembled for %d clusters via %s", len(evidence), p.evidence_provider)
        return evidence

    def _prepare_template(self, lens: str, evidence: dict[str, str]) -> str:
        if self.profile.prompt_template is None:
            raise ValueError("profile has no prompt_template")
        raw = Path(self.profile.prompt_template).read_text()
        return _fill_template(raw, lens=lens, profile=self.profile,
                              evidence_block=_build_evidence_block(evidence))

    # -- harmonization ---------------------------------------------------
    def _harmonize(self, consensus: dict[str, str]) -> dict[str, str]:
        p = self.profile
        guard_fn = make_reject_guard(p.reject_patterns) if p.reject_patterns else None
        out: dict[str, str] = {}
        for cid, raw in consensus.items():
            out[str(cid)] = harmonize_label(
                str(raw), p.vocab, p.synonyms,
                guard_fn=guard_fn, mode=p.vocab_mode,
                fallback_label=p.fallback_label, novel_prefix=p.novel_prefix,
                clean_novel=p.clean_novel,
            )
        return out

    # -- validation (gates the write, #9) --------------------------------
    def _validate(self, lens: str, consensus: dict[str, str], harmonized: dict[str, str],
                  marker_clusters: set[str]) -> dict:
        p = self.profile
        failures: list[str] = []
        warnings: list[str] = []

        missing = marker_clusters - set(harmonized)
        if missing:
            failures.append(f"clusters with no label: {sorted(missing)}")

        guard_fn = make_reject_guard(p.reject_patterns) if p.reject_patterns else None
        for cid, lab in harmonized.items():
            if guard_fn is not None and guard_fn(str(lab).lower()):
                failures.append(f"cluster {cid}: guard-rejected label leaked: {lab!r}")
            if is_novel_label(lab, p.novel_prefix):
                warnings.append(f"cluster {cid}: novel (open-vocab) label for review: {lab!r}")
            elif p.vocab_mode == "closed" and p.vocab and lab not in p.vocab and lab != p.fallback_label:
                failures.append(f"cluster {cid}: off-vocab label in closed mode: {lab!r}")

        return {
            "lens": lens,
            "n_clusters": len(marker_clusters),
            "n_labeled": len(harmonized),
            "failures": failures,
            "warnings": warnings,
            "status": "FAIL" if failures else "PASS",
        }

    # -- run --------------------------------------------------------------
    def run(
        self,
        *,
        markers_csv: str,
        evidence_csv: str | None,
        aux_path: str | None,
        species: str,
        tissue: str,
        models: list[str],
        api_keys: dict[str, str],
        lens: str,
        out_dir: str | Path,
        n_markers: int = 15,
        cache_dir: str | Path | None = None,
        consensus_model: str | None = None,
    ) -> AnnotationResult:
        p = self.profile
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        trace_root = out_dir / "trace"
        lens_trace = trace_root / lens
        lens_trace.mkdir(parents=True, exist_ok=True)
        if cache_dir is None:
            cache_dir = out_dir / "cache"
        cache_dir = Path(cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)

        # 1. markers + evidence
        markers = load_markers_csv(markers_csv, n_markers)
        marker_clusters = set(markers)
        evidence = self._assemble_evidence(evidence_csv, aux_path)

        # 2. prompt template (fill our slots, install the library global seam)
        template = self._prepare_template(lens, evidence)
        rendered = _prompt.render_prompt_preview(species, tissue, template, markers)
        log.info("lens=%s clusters=%d markers/cluster=%d prompt_chars=%d models=%s",
                 lens, len(markers), n_markers, len(rendered), models)

        # 3. trace: prompt + meta BEFORE the call
        tw = TraceWriter(lens, base_dir=lens_trace, git_cwd=self.git_cwd)
        tw.write_prompt(rendered)
        tw.write_meta(models=models, sources={"markers": str(markers_csv),
                                              "evidence": str(evidence_csv) if evidence_csv else "none",
                                              "aux": str(aux_path) if aux_path else "none"},
                      extra={"profile": p.name, "vocab_mode": p.vocab_mode,
                             "inject_evidence": p.inject_evidence})

        # 4. provider routing for token capture (model slug with '/' → OpenRouter)
        uses_openrouter = any("/" in m for m in models)
        capture = _capture.OpenRouterCapture() if uses_openrouter else _capture.GeminiTokenCapture()

        call_kwargs: dict[str, Any] = dict(
            marker_genes=markers, species=species, tissue=tissue, models=models,
            api_keys=api_keys,
            consensus_threshold=p.consensus_threshold,
            entropy_threshold=p.entropy_threshold,
            max_discussion_rounds=p.max_discussion_rounds,
            use_cache=True, cache_dir=str(cache_dir), verbose=True,
        )
        cm = consensus_model or p.consensus_model
        if cm:
            call_kwargs["consensus_model"] = cm

        # 5. the call — wrapped for token capture + DEBUG survival + template restore
        _debug.configure_llm_logging(lens, lens_trace)
        _prompt.install_prompt_template(template)
        try:
            with _debug.llm_debug_capture(lens, lens_trace), capture as cap:
                result = interactive_consensus_annotation(**call_kwargs)
            token_usage = cap.usage
        finally:
            _prompt.restore_prompt_template()

        consensus = {str(k): v for k, v in result.get("consensus", {}).items()}
        model_annotations = result.get("model_annotations", {})
        log.info("consensus clusters=%d controversial=%d", len(consensus),
                 len(result.get("controversial_clusters", [])))

        # 6. Python-recomputed metrics (#1) + harmonize (#8/#10)
        py_cp, py_entropy, py_majority = recompute_lens_metrics(model_annotations)
        harmonized = self._harmonize(consensus)

        # 7. trace: responses / discussion / tokens + cost AFTER
        tw.write_model_responses(model_annotations)
        tw.write_discussion(result.get("discussion_logs", {}))
        tw.write_tokens(token_usage)
        _cost.report_cost(token_usage, models)

        # 8. validation gates the write
        validation = self._validate(lens, consensus, harmonized, marker_clusters)
        (out_dir / "validation.json").write_text(json.dumps(validation, indent=2, default=str))
        for w in validation["warnings"]:
            log.warning(w)
        if validation["status"] == "FAIL":
            for f in validation["failures"]:
                log.error("VALIDATION FAIL: %s", f)
            raise SystemExit(1)

        # 9. write labels.csv (py_* are authoritative; llm_reported_* preserved separately)
        llm_cp = result.get("consensus_proportion", {})
        llm_ent = result.get("entropy", {})
        rows = []
        for cid in sorted(harmonized, key=_numeric_cluster_key):
            rows.append({
                "cluster": cid,
                "consensus_label": consensus.get(cid, ""),
                "harmonized_label": harmonized[cid],
                "py_consensus_proportion": py_cp.get(cid),
                "py_entropy": py_entropy.get(cid),
                "py_majority_label": py_majority.get(cid),
                "llm_reported_proportion": llm_cp.get(cid),
                "llm_reported_entropy": llm_ent.get(cid),
                "model_annotations": json.dumps(
                    {m: ann.get(cid) for m, ann in model_annotations.items() if cid in ann}),
            })
        labels_csv = out_dir / "labels.csv"
        pd.DataFrame(rows).to_csv(labels_csv, index=False)
        log.info("wrote %s (%d clusters)", labels_csv, len(rows))

        # 10. cost summary (recoverable from disk even on a cached re-run)
        _cost.write_cost_summary([lens], models, base_dir=trace_root, git_cwd=self.git_cwd)

        return AnnotationResult(
            lens=lens, labels=harmonized, out_dir=out_dir, labels_csv=labels_csv,
            validation=validation, status=validation["status"],
        )
