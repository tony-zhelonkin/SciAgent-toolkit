# mLLMCelltype — Python API Reference

Public exports (`from mllmcelltype import ...`):
`annotate_clusters`, `interactive_consensus_annotation`, `get_model_response`,
`get_provider`, `get_supported_providers`, `parse_marker_genes`,
`format_discussion_report`, `setup_logging`, `get_current_log_file`,
`clear_cache`, `get_cache_stats`.

---

> ## ⚠️ What this skill OVERRIDES in the library (read first)
>
> This page documents the **upstream `mllmcelltype` API**. `mllmct` deliberately does **not**
> trust two things the library emits, and recovers a third it discards:
>
> - **Consensus metrics are recomputed in Python.** The library's `consensus_proportion` and
>   `entropy` (see the return-dict table below) are **post-discussion** numbers — a stochastic
>   model's self-reported arithmetic describing agreement *after* the models were made to argue.
>   `mllmct` recomputes both deterministically from the per-model label multiset
>   (`core.consensus`) and writes them as the **authoritative** `py_consensus_proportion` /
>   `py_entropy` columns in `labels.csv`. The library's numbers are kept only as
>   `llm_reported_proportion` / `llm_reported_entropy` and **never override**. The LLM's chosen
>   `consensus` *label* is still authoritative; only the *metrics* are ours.
> - **Token usage + USD cost are captured by wrapping the provider call.** The library reads
>   only the response text and **discards** `usage_metadata` (Gemini) / the `usage`+`cost` block
>   (OpenRouter). `mllmct` wraps the provider call to recover them (and to force
>   `temperature=0, seed=0`). So "the library returns no token usage" is expected — `mllmct`
>   captures it out-of-band.
>
> See **`monkeypatch-internals.md`** for the four version-sensitive seams that make this work,
> and **`cell-state-annotation.md`** for the `py_*` vs `llm_reported_*` columns and the
> profile/evidence model.

---

## `annotate_clusters(...) -> dict[str, str]`

Single-model annotation. Returns `{cluster_id: label}` with **no** uncertainty metrics.

```python
annotate_clusters(
    marker_genes,            # dict[str, list[str]]  OR  pandas DataFrame with 'cluster' & 'gene' cols
    species,                 # "human" | "mouse" | ...
    provider="openai",       # "openai" | "anthropic" | "gemini" | "qwen" | "deepseek" |
                             #   "zhipu" | "stepfun" | "minimax" | "grok" | "openrouter"
    model=None,              # e.g. "gpt-5.5", "claude-sonnet-4-6"; None -> provider default
    api_key=None,            # falls back to <PROVIDER>_API_KEY env var
    tissue=None,             # tissue context, improves accuracy
    additional_context=None, # free-text hints (disease, protocol, expected lineages)
    prompt_template=None,    # override the built-in prompt
    use_cache=True, cache_dir=None,
    log_dir=None, log_level="INFO",
    base_urls=None,          # str or {provider: url} for proxies / self-hosted gateways
)
```

## `interactive_consensus_annotation(...) -> dict[str, Any]`

Multi-model consensus + iterative discussion. **Primary entry point.**

```python
interactive_consensus_annotation(
    marker_genes,                # dict[str, list[str]]
    species,
    models=None,                 # list of str names and/or {"provider":..., "model":...} dicts
    api_keys=None,               # {provider: key}; else env vars
    tissue=None,
    additional_context=None,
    consensus_threshold=0.7,     # agreement proportion below this => controversial
    entropy_threshold=1.0,       # Shannon entropy above this => controversial
    max_discussion_rounds=3,     # deliberation rounds for controversial clusters
    use_cache=True, cache_dir=None,
    verbose=False,
    consensus_model=None,        # optional dedicated model to arbitrate the final consensus
    base_urls=None,
    clusters_to_analyze=None,    # restrict to a subset of cluster IDs
    force_rerun=False,           # ignore cache for this run
)
```

### Return dict keys

| Key | Type | Meaning |
|---|---|---|
| `consensus` | `{cluster: label}` | Final cell-type label per cluster |
| `consensus_proportion` | `{cluster: float 0–1}` | Fraction of models agreeing on the final label. **Post-discussion / untrusted** — `mllmct` keeps this only as `llm_reported_proportion` and uses its own `py_consensus_proportion` (see top-of-page note). |
| `entropy` | `{cluster: float}` | Shannon entropy of model votes (higher = more disagreement). **Post-discussion / untrusted** — kept as `llm_reported_entropy`; `mllmct` uses `py_entropy` instead. |
| `controversial_clusters` | list | Clusters that triggered discussion |
| `resolved` | dict | Post-discussion labels for controversial clusters |
| `model_annotations` | `{model: {cluster: label}}` | Per-model raw predictions |
| `discussion_logs` | list | Structured deliberation transcript |

Map into AnnData:

```python
adata.obs["cell_type"]            = adata.obs["leiden"].astype(str).map(res["consensus"])
adata.obs["consensus_proportion"] = adata.obs["leiden"].astype(str).map(res["consensus_proportion"])
adata.obs["entropy"]              = adata.obs["leiden"].astype(str).map(res["entropy"])
```

## OpenRouter / mixed-model and free tiers

```python
res = interactive_consensus_annotation(
    marker_genes=marker_genes, species="human", tissue="blood",
    models=[
        {"provider": "openrouter", "model": "meta-llama/llama-4-maverick:free"},
        {"provider": "openrouter", "model": "deepseek/deepseek-v4-pro:free"},
        "claude-sonnet-4-6",            # mix native providers with OpenRouter
    ],
    consensus_threshold=0.7, max_discussion_rounds=2,
)
```

OpenRouter free tier (`:free` suffix) needs no credits but has rate limits
(~50 req/day, 20 req/min). Format: `'provider/model-name'`.

## Helpers

- `get_supported_providers()` → list of valid provider strings.
- `parse_marker_genes(...)` → normalize a DataFrame / file into the dict form.
- `format_discussion_report(res)` → human-readable Markdown of the deliberation.
- `clear_cache()` / `get_cache_stats()` → manage the response cache.
- `setup_logging(log_dir=..., level=...)` → only needed for a custom log location; logging auto-initializes otherwise.
