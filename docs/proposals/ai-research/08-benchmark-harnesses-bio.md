# 08 — Bio benchmark harnesses for ADR-006 tier-2

## BixBench

**Scope:** 205 questions derived from 60 real-world published Jupyter notebooks; 53 expert-curated analytical scenarios, 296 open-answer questions. From FutureHouse + ScienceMachine.

**Scoring:** majority vote for MCQ; LLM-as-judge for open-ended responses; exact-match for MCQ. Tests three capabilities: exploring biological datasets, multi-step computational analysis, contextual interpretation.

**Public availability:** Yes — dataset on Hugging Face, code at [github.com/Future-House/BixBench](https://github.com/Future-House/BixBench). License Apache-2.0.

**Integration via sciagent:**
- External agents integrate by implementing a `custom_rollout` function in `generate_trajectories.py` matching BixBench's trajectory format.
- CLI: `python bixbench/generate_trajectories.py --config_file <yaml>` then `python bixbench/postprocessing.py --config_file <yaml>`.
- YAML configs control model, prompts, batch size.
- **Cost: significant** — documentation explicitly warns "5 replicas across GPT-4o and Claude-3.5-Sonnet" yields "significant API costs" and 24–48 hours of wallclock. For a single agent + single config, scale roughly linearly (5–10 hours, a few hundred dollars API).
- Reported state-of-the-art is low — original BixBench paper reports ~17% best accuracy on open-answer, ~chance on MCQ for GPT-4o and Claude 3.5 Sonnet. Headroom is enormous; even mediocre sciagent performance will look credible against baselines.

**Verdict for sciagent:** the obvious tier-2 anchor. Standalone integration cost is real (1–2 days to wire up the trajectory format, plus API spend per run). High signal for a methods paper.

Sources:
- [BixBench arxiv 2503.00096](https://arxiv.org/html/2503.00096v1)
- [Future-House/BixBench repo](https://github.com/Future-House/BixBench)
- [BixBench announcement — FutureHouse](https://www.futurehouse.org/research-announcements/bixbench)
- [Future-House/data-analysis-crow (Finch agent framework)](https://github.com/Future-House/data-analysis-crow)

## OpenProblems.bio batch integration

**Scope:** atlas-level single-cell integration evaluated with 14 scIB metrics — graph connectivity, HVG conservation, isolated label F1, kBET, plus 10 others. Tracks ERA's headline single-cell result.

**Scoring:** weighted aggregate of metrics across batch-removal and bio-conservation. Public leaderboard.

**Public availability:** Yes. Task code at [github.com/openproblems-bio/task_batch_integration](https://github.com/openproblems-bio/task_batch_integration). Leaderboard at [openproblems.bio/results/batch_integration_feature](https://openproblems.bio/results/batch_integration_feature/). Built on Nextflow / Viash. scIB Python package at [scib.readthedocs.io](https://scib.readthedocs.io/).

**Integration via sciagent:**
- Method implementations follow OpenProblems' adapter pattern (Nextflow module wrapping a Python or R script).
- A sciagent run wraps "given an AnnData with `batch` and `label` obs columns, produce an integrated representation," then scIB scores it.
- Cost: moderate. The benchmark is designed for unattended runs. Nextflow + Viash is a learning curve if you haven't used them, but it's the standard pattern.
- The benchmark is exactly what ERA's single-cell claim is measured on — head-to-head comparison is natural.

**Verdict:** second tier-2 anchor specifically because it gives a direct ERA comparison. The infrastructure cost is the Nextflow/Viash dance — pencil in 3–5 days to wire up properly.

Sources:
- [OpenProblems batch integration task](https://github.com/openproblems-bio/task_batch_integration)
- [scIB benchmarking — Nature Methods 2022](https://www.nature.com/articles/s41592-021-01336-8)
- [scib reproducibility website](https://theislab.github.io/scib-reproducibility/)
- [scib python package docs](https://scib.readthedocs.io/)

## BioReason (VEP-Coding, VEP-Non-SNV)

**Scope:** variant effect prediction. VEP-Coding: 50,083 entries from ClinVar/gnomAD. VEP-Non-SNV: 36,088 ClinVar indels.

**Scoring:** F1. BioReason itself reports 80.0% F1 on VEP-Coding and 89.91% on VEP-Non-SNV, vs. max 49.19% / 66.51% for fDNA baselines and 39.58% / 76.21% for fLLM baselines.

**Public availability:** Yes (benchmarks in BioReason paper). Less mature harness than BixBench / OpenProblems — closer to "here is a labeled dataset, score predictions."

**Integration:** straightforward — variant + context → predicted effect label. No agent loop needed; this is more of a static benchmark than an agent task. Cheaper than BixBench (no multi-turn agent overhead).

**Verdict:** Optional tier-2. Useful for showing sciagent's literature-grounded interpretation works, but doesn't exercise the orchestrator-skill graph deeply. Could be deferred.

Sources:
- [BioReason: Incentivizing Multimodal Biological Reasoning — themoonlight review](https://www.themoonlight.io/en/review/bioreason-incentivizing-multimodal-biological-reasoning-within-a-dna-llm-model)
- [Variant Effect Predictors — varianteffect.org](https://www.varianteffect.org/veps/)

## GIFT-Eval

**Scope:** time-series forecasting, 23 datasets, 144K time series, 177M data points, **7 domains: Web/CloudOps, Finance, Energy, Traffic, Weather, Economics, Sensor Data.**

**Bio relevance: none.** No biology, epidemiology, healthcare, genomics, or single-cell datasets. The spec's framing of GIFT-Eval as a sciagent tier-2 benchmark is a category mistake — it's a general time-series forecasting benchmark, not a bio benchmark. ERA evaluated against it because ERA is domain-agnostic; sciagent is bio-focused.

**Public availability:** Yes, [github.com/SalesforceAIResearch/gift-eval](https://github.com/SalesforceAIResearch/gift-eval), HuggingFace dataset. Minutes to set up. Uses GluonTS evaluation patterns.

**Verdict for sciagent: drop from tier-2.** Worth running only if sciagent ships a forecasting-overlay role. For the bio-anchored MVP, this is noise.

Sources:
- [GIFT-Eval arxiv 2410.10393](https://arxiv.org/pdf/2410.10393)
- [SalesforceAIResearch/gift-eval](https://github.com/SalesforceAIResearch/gift-eval)
- [GIFT-Eval HuggingFace](https://huggingface.co/datasets/Salesforce/GiftEval)

## PaperQA2 LitQA2

**Scope:** literature QA on full-text scientific papers. From FutureHouse. PaperQA2 reports precision 85.2% / accuracy 66% on LitQA2; superhuman vs. PhD/postdoc biologists.

**Public availability:** Yes — [Future-House/paper-qa](https://github.com/Future-House/paper-qa), [Future-House/litqa](https://github.com/Future-House/litqa). MIT-licensed.

**Integration:** moderate. PaperQA2 is itself a pipeline; LitQA2 is the labeled QA dataset. To benchmark a sciagent role, the role plays the QA role, PaperQA2 (or its evaluator) scores. Pencil in 1–2 days.

**Verdict:** good fit for *literature-search* skills specifically. If sciagent ships a `lit-research` role or skill (which the current roles do not — `bioinf-librarian` agent is the closest), wire it in. Otherwise low priority.

Sources:
- [PaperQA2 / WikiCrow announcement — FutureHouse](https://www.futurehouse.org/research-announcements/wikicrow)
- [Future-House/paper-qa repo](https://github.com/Future-House/paper-qa)
- [LitQA Benchmark](https://github.com/Future-House/litqa)

## Order-of-magnitude integration cost summary

| benchmark | wire-up days | per-run API $ | per-run wallclock | bio fit | priority |
|---|---|---|---|---|---|
| BixBench | 1–2 | $200–800 | 5–48 h | exact | tier-2 anchor |
| OpenProblems batch-integration | 3–5 | <$50 (mostly GPU) | hours | exact | tier-2 anchor (ERA head-to-head) |
| BioReason VEP | 0.5–1 | $20–100 | <1 h | exact (variant-level) | optional |
| GIFT-Eval | 0.5 | <$50 | <1 h | **none** | drop |
| PaperQA2 LitQA2 | 1–2 | $50–200 | hours | exact (literature) | optional |

## Recommendation

For ADR-006 tier-2 MVP, ship BixBench + OpenProblems batch-integration. Defer BioReason and PaperQA2 to a v2. Remove GIFT-Eval from the tier-2 list.
