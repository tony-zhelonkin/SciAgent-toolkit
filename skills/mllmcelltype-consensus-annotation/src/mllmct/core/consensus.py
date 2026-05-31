"""Deterministic, Python-recomputed consensus metrics (remediation fix #1).

mLLMCelltype asks the LLM to *emit* consensus_proportion and Shannon entropy as text,
parses them verbatim, and overwrites them with the POST-DISCUSSION round's numbers.
Those are (a) a stochastic model's self-reported arithmetic and (b) describe agreement
AFTER the models were made to argue, masking the original split. We instead compute
both as exact functions of the per-model label multiset. The LLM's chosen ``consensus``
*label* is kept as authoritative; only the *metrics* are ours.

Ported from ``cellstate_llm.compute_consensus_metrics`` / ``recompute_lens_metrics`` —
already domain-agnostic.
"""

from __future__ import annotations

import math
from collections import Counter


def compute_consensus_metrics(
    per_model_labels: dict[str, str],
) -> tuple[str, float, float]:
    """Deterministically recompute ``(majority_label, consensus_proportion, entropy)``.

    SOURCE OF TRUTH: the per-model labels for ONE cluster
    (``result['model_annotations'][model][cluster]``).

        consensus_proportion = max_label_count / n_models
        entropy (Shannon)    = -Σ p_i · log2(p_i)   over the label distribution

    Ties for the majority resolve deterministically by (count desc, label asc). Returns
    ``('', nan, nan)`` only if no labels are present.
    """
    labels = [str(v) for v in per_model_labels.values() if v is not None and str(v).strip()]
    n = len(labels)
    if n == 0:
        return ("", float("nan"), float("nan"))

    counts = Counter(labels)
    majority_label, majority_count = sorted(
        counts.items(), key=lambda kv: (-kv[1], kv[0])
    )[0]

    consensus_proportion = majority_count / n
    entropy = -sum((c / n) * math.log2(c / n) for c in counts.values() if c > 0)
    entropy = abs(entropy)  # guard against -0.0 from log2
    return (majority_label, consensus_proportion, entropy)


def recompute_lens_metrics(
    model_annotations: dict[str, dict[str, str]],
) -> tuple[dict[str, float], dict[str, float], dict[str, str]]:
    """Vectorize ``compute_consensus_metrics`` over all clusters of one annotation run.

    Args:
        model_annotations: ``{provider:model -> {cluster_id -> label}}`` (the
            ``model_annotations`` key of an mllmcelltype result / saved JSON).

    Returns:
        ``(consensus_proportion, entropy, majority)`` dicts keyed by cluster id.
        ``majority`` is the deterministic majority label (informational; the LLM's
        ``consensus`` label is kept as authoritative per fix #1).
    """
    cluster_ids: set[str] = set()
    for per_cluster in model_annotations.values():
        cluster_ids.update(str(k) for k in per_cluster.keys())

    cp: dict[str, float] = {}
    ent: dict[str, float] = {}
    maj: dict[str, str] = {}
    for cid in cluster_ids:
        per_model = {
            model: str(per_cluster[cid])
            for model, per_cluster in model_annotations.items()
            if cid in per_cluster
        }
        m_label, m_cp, m_ent = compute_consensus_metrics(per_model)
        cp[cid] = m_cp
        ent[cid] = m_ent
        maj[cid] = m_label
    return cp, ent, maj
