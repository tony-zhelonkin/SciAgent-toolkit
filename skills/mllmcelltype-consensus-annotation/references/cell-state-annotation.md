# Cell-STATE annotation — how it differs, and how to adapt it without touching code

> Notes for when I want functional-STATE calls (not identity) and need to point the tool at a
> *new* dataset. The whole adaptation is YAML + (maybe) one small EvidenceProvider class — no
> edits to `engine.py` or `core/`. If you just need plain cell-TYPE labels, use
> `profiles/celltype.yaml` and skip this.

## TYPE vs STATE in one breath

Both modes run the **same** code path (`engine.AnnotationEngine.run`); the active **profile**
decides everything else.

| | cell-TYPE | cell-STATE |
|---|---|---|
| What the LLM is asked | "what IS this cluster?" (identity/lineage) | "what is this cluster DOING?" (functional state) — identity is FIXED |
| Prompt input | marker genes only | markers **+ injected per-cluster evidence** |
| `inject_evidence` | `false` | `true` (runs an `EvidenceProvider`) |
| Vocab | usually open / empty (free-text passthrough) | a curated `vocab` + `synonyms`, still `open` so novel states survive |
| Guards | none | `reject_patterns` so a lineage name can never be emitted as a state |
| Template | `templates/celltype_prompt.txt` | `templates/cellstate_prompt.txt` (adds vocab/guardrails/evidence slots) |

Switching modes = passing a different `--profile`. Nothing else.

---

## The profile knobs (the whole adaptation surface)

A profile is a YAML file loaded by `profile.load_profile()` and merged over built-in defaults
(`profile.py`). The cell-state-relevant knobs:

| Knob | Meaning |
|---|---|
| `inject_evidence: true` | engine instantiates the `EvidenceProvider` and renders an evidence block; `false` ⇒ markers-only |
| `evidence_provider` | `"module.path:ClassName"` resolved by `core.evidence.load_provider` (must `import`-resolve in the locked venv) |
| `evidence_options` | a free-form dict handed to the provider constructor — **this replaces importing any project `config.py`** |
| `aux_table_key` | which CLI input (`--aux`) feeds the provider's `aux_path` (e.g. a program-decode map) |
| `prompt_template` | path to the template with the `<<...>>` slots (resolved absolute → skill_root → profile dir) |
| `vocab` | preferred/canonical label list (rendered into `<<vocab>>`) |
| `synonyms` | `{substring: canonical}`, longest-key-wins substring match |
| `vocab_mode` | `open` → unmapped kept as `Novel:*`; `closed` → snapped to `fallback_label` |
| `fallback_label` | the "I can't tell" bucket (default `Ambiguous_LowSignal`) |
| `novel_prefix` / `clean_novel` | how a kept-novel label is tagged + normalized (`Novel:DNA_Damage_Response`); set `novel_prefix:""`+`clean_novel:false` for verbatim passthrough (cell-type) |
| `domain_guards.reject_patterns` | regexes matched against the **lowercased** label → forces `fallback_label` (so a fixed lineage term can't leak in as a state) |
| `domain_guards.restricted_states` | `{state: [allowed_group,...]}` → a warning (not a fail) when a state lands outside its allow-list |
| `domain_guards.guardrails_text` | free text dropped into the `<<guardrails>>` prompt slot |
| `join` | consumed ONLY by `mllmct reconcile` (see two-axis section); `annotate` ignores it |

> Note the **shipped `cellstate.yaml` is a GENERIC neutral example** — textbook functional
> states + abstract axis names, wired to the synthetic test fixtures. It is a template, not a
> finished config. Your dataset's real biology (`vocab`, `synonyms`, `axis_order`,
> `reject_patterns`, the guardrails text) belongs in **your own analysis-repo profile**, copied
> from this one — not committed back into the shared skill.

### The template's `<<...>>` slots

`templates/cellstate_prompt.txt` uses two layers of placeholders:
- `<<lens>>`, `<<vocab>>`, `<<fallback_label>>`, `<<mode_clause>>`, `<<guardrails>>`,
  `<<evidence_block>>` — filled by the **engine** (`engine._fill_template`) from the profile;
- `{species}`, `{tissue}`, `{markers}` — left for the **library's** `create_prompt` to fill.

`<<mode_clause>>` is auto-set from `vocab_mode` (an "you MAY propose a novel name" clause for
open, a "do NOT invent labels" clause for closed). Use `mllmct preview-prompt …` to see the
exact bytes the model will receive — no API call.

---

## Writing your own EvidenceProvider

The ABC lives in `core/evidence.py`:

```python
class EvidenceProvider(ABC):
    def __init__(self, options: dict | None = None, aux_path: str | None = None): ...
    @abstractmethod
    def assemble(self, evidence_csv: str) -> dict[str, str]:
        """Return {cluster_id: evidence_string}. Empty dict ⇒ markers-only."""
```

Contract:
- `options` is exactly `profile.evidence_options` (your knobs — read from here, never from a
  global `config.py`, so the provider stays unit-testable from fixtures).
- `aux_path` is the optional `--aux` table (a secondary file, e.g. a program-decode map).
- `assemble(evidence_csv)` reads your per-cluster evidence table and returns
  `{cluster_id: evidence_string}` with **string** cluster ids (matched to `adata.obs`
  categories). Return `{}` to fall back to markers-only.
- The engine sorts evidence numerically and renders one `  cluster N: <string>` line per
  cluster into `<<evidence_block>>`.

Point a profile at it with `evidence_provider: "your.module:YourProvider"`. The class must
`import`-resolve inside the locked `.venv` (so either ship it under `src/mllmct/plugins/` or
make your module importable on the venv's path) and subclass `EvidenceProvider`.

---

## The reference `PanelEvidenceProvider` — four mechanics

`plugins/panel_evidence.py` is the worked, **domain-neutral** example (no real biology). It
demonstrates the four assembly mechanics that made the original cell-state run trustworthy.
Every name/threshold comes from `options`:

1. **Program decode.** An opaque program id (`P7`) is decoded via the `aux_path` program-map to
   `program=P7[Category:gene1,gene2,…]`, so the LLM never sees a bare id. Unmapped → `(unmapped)`;
   a program in a `drop_categories` category is omitted, or rendered `(artifact)` when
   `flag_artifact_not_drop:true`.
2. **Signature de-saturation.** Globally-dominant `always_on_signatures` are pulled OUT of the
   ranked top-`signature_top_k` and surfaced as compact `name+` flags, so a constant anchor
   can't dominate the evidence.
3. **Binned-axis panel with companion pairing.** Only NON-`baseline_bins` axes are shown, capped
   to `max_axes` ranked by `|z|`, in the stable `axis_order`. A `companion_axis` is
   **force-rendered** whenever its `companion_trigger_axis` fires (the "never show X without its
   companion" rule), even if it would otherwise be capped or baseline.
4. **Category + flag tokens.** A dominant `category` token, plus a thresholded boolean
   `flag_metric` rendered with `flag_labels`.

**`evidence_csv` columns it expects** (all column names configurable via `options`):
`cluster`, `dominant_category`, `dominant_program`, `signatures_top` (`;`-joined),
`companion_bin`, `flag_metric`, and per-axis `bin_<axis>` / `z_<axis>` pairs.
Optional `aux_path` program-map columns: `program`, `category`, `top_genes` (`;`-joined).

> The original (real) provider behind this is the Th1/Th17 cNMF + AUCell assembler in
> `02_Analysis/helpers/cellstate_llm.py`. `PanelEvidenceProvider` is the de-projectified
> generalization of it: cNMF program decode → mechanic 1, AUCell de-saturation → mechanic 2,
> the binned-z temperature/axis panel → mechanic 3.

---

## Authoritative metrics: `py_*` vs `llm_reported_*`

`labels.csv` carries **both** sets of uncertainty numbers:
- `py_consensus_proportion` / `py_entropy` / `py_majority_label` — **authoritative**,
  recomputed in Python (`core.consensus`) directly from the per-model label multiset.
- `llm_reported_proportion` / `llm_reported_entropy` — the library's own emitted numbers,
  **preserved but never trusted**: they describe agreement AFTER the discussion round and are a
  stochastic model's self-reported arithmetic. See `api-reference.md` and
  `monkeypatch-internals.md` for why.

Always read the `py_*` columns. `backfill-metrics` can recompute them offline from
`model_annotations` if you ever need to.

---

## Open-vocab `Novel:*` are HYPOTHESES, not findings

In `vocab_mode: open`, a label that doesn't map to `vocab`/`synonyms` (and isn't guard-rejected)
is kept as `Novel:<normalized>` rather than being forced into the nearest vocab term — this is
deliberately discovery-safe. But a single-run novel discovery is **sampling-sensitive**: even
with `temperature=0`+`seed`, the panel composition and discussion path can surface or drop a
novel state run-to-run. Treat every `Novel:*` as a **hypothesis needing orthogonal validation**
(an independent marker check, a held-out scoring axis, a re-run with a different panel) before
it becomes a claim. The validation gate flags novels as **warnings**, not failures, precisely
so you remember to look at them.

---

## Two-axis reconcile (`mllmct reconcile`)

**When to use it.** Only when each cell carries **two orthogonal annotations** — e.g. the same
clusters annotated under two different lenses/embeddings (call them axis A and axis B) — and
collapsing them to one consensus would throw away biology. It is a deterministic, no-LLM,
no-API post-processing join (`reconcile.py`), gated by `profile.join.enabled`.

**The `JoinConfig` fields** (from `profile.join`, via `JoinConfig.from_profile_join`):

| Field | Meaning |
|---|---|
| `distinct_labels` | axis-A labels that are "distinct enough" to KEEP both axes as a compound `A · B` |
| `order` | canonical priority for naming a two-axis pair (guarantees `{A,B}` renders ONE string, fixing the old `A · B`/`B · A` duplication) |
| `uninformative_labels` | labels carrying no state info → treated as missing |
| `conf_floor` | a per-axis confidence below this → that axis is treated as uninformative |
| `fallback_label` | used when both axes are uninformative |
| `sep` | the join separator (default ` · `, U+00B7) |

The per-cell join rule (`join_cell`): both uninformative → fallback; exactly one informative →
that one; agree → the shared label; disagree & A distinct → canonical `A · B`; disagree & A
generic → defer to the finer B label.

**The per-cell axes table is YOUR job to build.** `reconcile` consumes a CSV with columns
`axis_a_label, axis_a_conf, axis_b_label, axis_b_conf` (one row per cell). Producing that table
from your AnnData — mapping each cell's two cluster assignments to their respective `labels.csv`
states + confidences — is dataset-specific and lives in your analysis repo, **not** in this
skill. `reconcile` only does the deterministic join + crosstab once you hand it that table.
