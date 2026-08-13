# 06 — ERA: "An AI system to help scientists write expert-level empirical software"

> **Outcome note (superseded).** Superseded by the three-verb toolkit
> architecture documented in [`docs/architecture.md`](../../architecture.md).
> Retained as a dated research record; the original body below is unchanged.

Grounds ADRs 004 (`scorable` overlay) and 009 (uv-pinned execution).

## The paper exists; the spec's attribution is mostly correct

- **Title:** *An AI system to help scientists write expert-level empirical software*
- **Authors:** Aygün, E.; Belyaeva, A.; Comanici, G. et al. (~40 authors, Google DeepMind + collaborators; full list includes Cory Y. McLean, Subhashini Venugopalan, Michael P. Brenner, James Manyika, etc.)
- **Venue:** *Nature*, May 19 2026
- **DOI:** 10.1038/s41586-026-10658-6
- **Preprint:** [arXiv:2509.06503](https://arxiv.org/abs/2509.06503) (also referenced as 2025-09-09 — preprint was up before journal acceptance)
- **Paper link:** [nature.com/articles/s41586-026-10658-6](https://www.nature.com/articles/s41586-026-10658-6) (paywalled; redirects to login)
- **DeepMind blog:** [Empirical Research Assistance (ERA): From Nature publication to catalyzing Computational Discovery](https://research.google/blog/empirical-research-assistance-era-from-nature-publication-to-catalyzing-computational-discovery/)

Spec attribution as "Aygün et al., Nature 2026" is correct on first author and year. It is part of Google DeepMind's [Gemini for Science](https://blog.google/innovation-and-ai/technology/research/gemini-for-science-io-2026/) initiative announced at Google I/O May 2026 — context the spec omits but worth knowing for grilling.

## Code repository — exists and is open source

[github.com/google-research/era](https://github.com/google-research/era) — Apache-2.0.

Structure:
- `implementation/` — core algorithm + examples (e.g., `playground_s3e1.py` for regression)
- `era_applications/` — real-world scientific applications (single-cell, COVID forecasting, etc.)
- `docs/` — documentation
- 90.6% Jupyter Notebook, 9.4% Python

Dependencies: pandas, numpy, scikit-learn, google-generativeai. Python 3.10+.

The link `github.com/google-research/era/tree/main/era_applications` referenced in the DeepMind blog resolves correctly. The spec's mention of `github.com/google-research/era` in skill 2.1 (`tree-search-controller`) is verified.

## Methodology — partial confirmation

The spec describes "ERA-style tree search" with "mutator + scorer + selector" prompts and "UCB1 selection." Repository inspection clarifies:

- **Tree search algorithm:** **Flat UCB Tree Search (FUTS)**, in `implementation/futs.py`. Uses Flat UCB (PUCT-style) — not classical UCB1. The spec saying "UCB1" is approximately right but technically loose. The constant `c` in `score + c·sqrt(ln(N_total)/N_node)` (skill 2.1) is from UCB1; FUTS uses a similar but PUCT-derived formula.
- **Mutator:** `generate_fn(problem, past_solution) → new_solution`. Prompts an LLM. Stateless per call. Matches spec.
- **Scorer:** `execute_fn(code) → float`. "Runs in a sandboxed environment against the problem's metric." Matches spec.
- **Selector:** present in the FUTS implementation; not separately documented in the README excerpt.

The spec's three-agent decomposition (mutator/scorer/selector) is accurate to the architecture, though the names map loosely.

## Empirical claims — partially confirmed

| spec claim | verification |
|---|---|
| Beats hand-written methods on OpenProblems.bio (batch integration) | Paper claims "40 novel methods for single-cell data analysis that outperformed the top human-developed methods on a public leaderboard." OpenProblems is the implied leaderboard but not named in the abstract. |
| Beats CDC ensemble for COVID forecasting | Confirmed: "14 models that outperformed the CDC ensemble and all other individual models for forecasting COVID-19 hospitalizations." |
| GIFT-Eval improvement | Not in abstract. Paper claims time-series forecasting results across multiple domains. GIFT-Eval is not named in publicly accessible material. |
| Other domains | Paper claims: genomics, epidemiology, geospatial analysis, neural activity prediction, time-series forecasting, numerical analysis. |

The OpenProblems.bio result is plausibly the headline single-cell story. The "GIFT-Eval" mention in the spec is **unverifiable from publicly accessible material** — may be in the paper body, may not. Worth confirming before any downstream claim depends on it.

## Sandbox / uv discussion

The repo README mentions "sandboxed environment" for `execute_fn` but does not explicitly reference `uv`. ADR-009's claim that the Antigravity Skills paper / ERA paper "reports significant reproducibility loss" and "adopted `uv` to pin environments per skill execution" is **not confirmed in publicly accessible ERA material**. The full paper text (gated) may discuss it; the abstract, blog post, and repo README do not. This is a spec claim that needs source-paper verification.

Two possibilities:
- (a) The claim is in the paper body, just not in the abstract/blog/README.
- (b) The claim conflates uv usage in some other Google source (e.g., Gemini for Science skills) with ERA.

ADR-009 is consequential (uv as a hard dependency) — its empirical justification deserves a real citation, not a possibly-misremembered one.

## Per-search token / compute cost

Not reported in abstract or blog. The paper (78 pages, 31 figures per ResearchGate) likely contains compute budgets but they are not in publicly accessible summaries. Skill 2.1's recommendation of "100-1000 candidate evaluations" and "budgets <50 nodes — overhead dominates" are reasonable but not sourced from ERA.

## Bottom line for ADRs 004 and 009

- **ADR-004 is safe.** ERA is real, the code is open, the architecture is roughly as described. The `scorable` overlay can borrow patterns directly from `implementation/futs.py`.
- **ADR-009 is exposed.** The empirical justification for `uv` as a hard dependency may not survive direct fact-checking against the ERA paper. The personal-take section in ADR-009 already flags this. Re-read the paper body before committing to "uv-pinned execution" as a Day-1 architecture requirement.
- **Naming nit:** "Flat UCB Tree Search (FUTS)" is the actual algorithm. The skill body should say so rather than "UCB1" — this is a citation precision question, not a design question.

Sources:
- [arxiv 2509.06503 — *An AI system to help scientists write expert-level empirical software*](https://arxiv.org/abs/2509.06503)
- [google-research/era repository (Apache-2.0)](https://github.com/google-research/era)
- [DeepMind blog — ERA from Nature publication](https://research.google/blog/empirical-research-assistance-era-from-nature-publication-to-catalyzing-computational-discovery/)
- [Gemini for Science launch — Google blog](https://blog.google/innovation-and-ai/technology/research/gemini-for-science-io-2026/)
- [Tech Times coverage — ERA beats CDC](https://www.techtimes.com/articles/316901/20260520/gemini-science-launches-peer-reviewed-benchmarks-era-beat-cdc-forecasting-model.htm)
- [Harvard SEAS news coverage](https://seas.harvard.edu/news/ai-system-automates-coding-scientific-research)
- [Nature paper landing page (paywalled)](https://www.nature.com/articles/s41586-026-10658-6)
