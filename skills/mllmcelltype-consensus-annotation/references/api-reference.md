# mLLMCelltype — Python API Reference

Public exports (`from mllmcelltype import ...`):
`annotate_clusters`, `interactive_consensus_annotation`, `get_model_response`,
`get_provider`, `get_supported_providers`, `parse_marker_genes`,
`format_discussion_report`, `setup_logging`, `get_current_log_file`,
`clear_cache`, `get_cache_stats`.

---

> ## ⚠️ What this skill OVERRIDES in the library (read first)
>
> This page documents the **upstream `mllmcelltype` API**. `mllmct` distrusts the library's
> consensus metrics and recovers token/cost the consensus entrypoint leaves unreachable:
>
> - **Consensus metrics are recomputed in Python.** The returned `consensus_proportion` /
>   `entropy` are **post-discussion**, self-reported numbers. `mllmct` recomputes both from the
>   per-model label multiset (`core.consensus`) as the authoritative `py_consensus_proportion` /
>   `py_entropy` in `labels.csv`, keeping the library's as `llm_reported_*` (never override). The
>   chosen `consensus` *label* stays authoritative; only the *metrics* are ours.
> - **Token usage + USD cost are captured by wrapping the provider call.** 2.0.7 reads usage
>   natively into a `usage_sink`, but that sink is **not** threaded through
>   `interactive_consensus_annotation` — so via the entrypoint `mllmct` uses, tokens/cost stay
>   unreachable. The wrap recovers them and forces `temperature=0, seed=0` (the library exposes no
>   public determinism knob). Custom prompts, by contrast, are threaded natively (`prompt_template=`).
>
> See **`monkeypatch-internals.md`** for the version-sensitive seams, and
> **`cell-state-annotation.md`** for the `py_*` vs `llm_reported_*` columns and profile/evidence model.

---

## `annotate_clusters(...) -> dict[str, str]`

Single-model annotation. Returns `{cluster_id: label}` with **no** uncertainty metrics.

```python
annotate_clusters(
    marker_genes,            # dict[str, list[str]]  OR  pandas DataFrame with 'cluster' & 'gene' cols
    species,                 # "human" | "mouse" | ...
    provider="openai",       # "openai" | "anthropic" | "gemini" | "qwen" | "deepseek" |
                             #   "zhipu" | "stepfun" | "minimax" | "grok" | "openrouter"
    model=None,              # e.g. "gpt-5", "claude-sonnet-4.6"; None -> provider default
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
    prompt_template=None,        # custom prompt; validated + threaded to create_prompt (native since 2.0.7)
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
| `consensus_proportion` | `{cluster: float 0–1}` | Fraction agreeing on the final label. Post-discussion — kept as `llm_reported_proportion` (see note). |
| `entropy` | `{cluster: float}` | Shannon entropy of votes (higher = more disagreement). Post-discussion — kept as `llm_reported_entropy` (see note). |
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
        "claude-sonnet-4.6",            # mix native providers with OpenRouter
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
