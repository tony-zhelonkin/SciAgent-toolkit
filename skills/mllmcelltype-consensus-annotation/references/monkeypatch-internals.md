# Monkeypatch internals — the version-sensitive seams

> Deep notes for when I (or whoever inherits this) need to understand *why* `mllmct`
> reaches into `mllmcelltype`'s guts, where exactly it binds, and what silently rots if a
> dependency moves. If you just want to *run* annotation, this file is not for you — read
> `cell-state-annotation.md` or the SKILL.md instead.

`mllmcelltype==2.0.7` does its job (multi-LLM consensus over marker lists) but throws away
three things I care about: token usage, sampling determinism, and the raw model responses.
None are exposed through a public API, so I recover them by wrapping internal seams *without
editing the vendored package*. Every wrapper binds to an **attribute name + call signature**,
not a stable public contract — which is why the target packages are pinned `==` and gated by
`checks/smoke_check_versions.py`.

A fourth concern — injecting a custom (cell-STATE) prompt template — *used* to need a
monkeypatch (seam (i), the `DEFAULT_PROMPT_TEMPLATE` global swap) but is **native since 2.0.7**:
`interactive_consensus_annotation` now accepts `prompt_template=` and threads it to
`create_prompt`. The engine passes the filled template through directly (`engine.py`), so seam
(i) is retired; the smoke-check now asserts the native parameter instead of the global.

The wrappers live in:
- `core/capture.py` — Gemini `generate_content` wrap + OpenRouter `requests.post` wrap
- `core/debug_capture.py` — logger DEBUG survival
- `core/prompt.py` — **no longer a monkeypatch**; renders the byte-faithful prompt preview via
  the library's own `create_prompt(prompt_template=…)`

The three surviving wrappers were lifted from the original real-world implementation in
`02_Analysis/helpers/cellstate_llm.py` / `cellstate_obs.py`. Line numbers below are verified
against the locked `.venv` copy of `mllmcelltype 2.0.7`; treat them as "true as of 2.0.7" and
re-verify on any bump.

---

## Version-sensitivity table

| # | Seam (module.attr) | Verified location (2.0.7) | Pin that guards it | If the seam moves… |
|---|---|---|---|---|
| i | ~~`mllmcelltype.prompts.DEFAULT_PROMPT_TEMPLATE` global swap~~ **RETIRED — native since 2.0.7** | native param `consensus.py:2598` (`prompt_template=`), validated `:2651`, threaded `:2676` | `mllmcelltype==2.0.7` | custom template silently ignored → models get the stock cell-TYPE prompt; smoke-check asserts the param exists |
| ii | `google.genai.models.Models.generate_content`; `GenerateContentResponse.usage_metadata`; `GenerateContentConfig.{temperature,seed}` | hardcoded config `providers/gemini.py:128`; native usage now at `gemini.py:34`/`:136` (unsurfaced via the consensus entrypoint) | `google-genai==2.6.0` (+`pydantic==2.13.3`) | tokens lost (cost=$0 everywhere) and/or determinism lost (re-runs drift); no error |
| iii | `mllmcelltype.providers.openrouter.requests.post`; the chat-completions body/parser that omits `temperature` and drops `usage`/`cost` | `providers/openrouter.py:5` (`import requests`), `:61` (`post_func=requests.post`); body `providers/common.py:270` (`data=json.dumps`, `request_json=False` default `:448`); no temperature sent (`common.py:89,102-103`) | `requests==2.33.1` | OpenRouter tokens + native USD lost; `temperature=0`/`seed` not injected → non-deterministic; no error |
| iv | `mllmcelltype.logger.setup_logging` "just update level" fast path; logger name `"llmcelltype"` | fast path `logger.py:75-78`; logger name `logger.py:9` | `mllmcelltype==2.0.7` | raw-response DEBUG lines never land (empty `llmcelltype_debug.log`); no error |

The unifying failure mode: **every one fails *silently***. The run completes, labels come out,
nothing throws — you just quietly lose tokens, reproducibility, your prompt, or your debug trace.
That's why the gate is structural-assertion-on-lock, not a runtime try/except.

---

## (i) Prompt-template injection — RETIRED (native since 2.0.7)

**History.** In `2.0.5`, `interactive_consensus_annotation` did not thread a `prompt_template`
argument; `create_prompt` fell back to the module global
`mllmcelltype.prompts.DEFAULT_PROMPT_TEMPLATE`, read *at call time*. The skill injected its
cell-STATE prompt by swapping that global (`install_prompt_template` / `restore_prompt_template`
in `core/prompt.py`, wrapped in a try/finally around the call).

**What 2.0.7 gives natively.** `interactive_consensus_annotation(prompt_template=…)`
(`consensus.py:2598`) validates the template (`validate_prompt_template`, `:2651`) and threads
it to `create_prompt` (`:2676`). The engine now sets `call_kwargs["prompt_template"] = template`
and calls straight through — no global swap, no restore. `core/prompt.py` keeps only
`render_prompt_preview`, which passes `prompt_template=` to the library's own `create_prompt`
for the byte-faithful trace.

**One behavioral note from the bump.** `create_prompt` (and `validate_prompt_template`) now run
`_get_prompt_template_fields`, which rejects any `{placeholder}` outside
`SUPPORTED_PROMPT_PLACEHOLDERS = (species, tissue, markers, context)` (`prompts.py:72`). The
engine fills its own `<<…>>` slots *before* the library sees the template, leaving only
`{species}/{tissue}/{markers}`, so the shipped templates validate — but any literal `{`/`}`
injected via evidence or guardrail text would now raise `Invalid prompt_template format` (it did
not in 2.0.5). This applies to the native and the (removed) global-swap paths equally; it is a
property of the 2.0.7 bump, caught by `test_end_to_end.py`.

**Guarded by:** `mllmcelltype==2.0.7`. Smoke-check asserts `create_prompt` **and**
`interactive_consensus_annotation` each accept a `prompt_template` parameter.

---

## (ii) Gemini `generate_content` wrap

**What the library does/discards.** The Gemini provider (`providers/gemini.py`) calls
`client.models.generate_content(...)` with a hardcoded
`config=types.GenerateContentConfig(temperature=0.7, max_output_tokens=4096)` (`gemini.py:128`)
— **no seed**, so re-runs drift. 2.0.7 added native usage extraction (`extract_gemini_usage`,
`gemini.py:34`; `capture_usage(response, usage_sink, …)`, `gemini.py:136`), but `usage_sink` is
**not** threaded through `interactive_consensus_annotation`, so a caller on the consensus
entrypoint still cannot retrieve tokens without the wrap. Determinism is still hardcoded to
`temperature=0.7` with no public knob.

**The seam we bind to.** `GeminiTokenCapture` (`core/capture.py`) monkeypatches the class
method `google.genai.models.Models.generate_content`. The wrapper does two things:

1. **Determinism (merge, don't clobber).** `force_deterministic_config` rewrites the incoming
   `config` to `temperature=0` and `seed=0`, *preserving every other field* (notably
   `max_output_tokens=4096`). It handles three shapes:
   - a `GenerateContentConfig` **pydantic BaseModel** → `model_copy(update=overrides)`
     (pydantic v1 fallback: `copy(update=...)`); `seed` is only added if the model class
     actually declares a `seed` field (read off `cls.model_fields`, on the class — instance
     `.model_fields` access is deprecated in pydantic ≥2.11);
   - a plain `dict` (`GenerateContentConfigDict`) → copy + key override (checked **first**,
     because dict also has `.copy`);
   - `None` / anything unrecognised → returned unchanged (never break the call).
   The config travels as a kwarg in 2.0.7, but the wrapper also handles a positional
   `config` (signature `generate_content(model, contents, config)`).
2. **Token capture.** After the real call, it reads `response.usage_metadata` and accumulates
   `prompt_token_count` / `candidates_token_count` / `total_token_count` per model into
   `self.usage = {model: {prompt, output, total, calls}}`.

**`seed` is best-effort.** Per Google's own field docstring, `seed` makes the model
*best-effort* deterministic, not bitwise-guaranteed. `temperature=0` is the real lever.

**Silent breakage.** If `Models.generate_content` is renamed/moved, or `usage_metadata` is
renamed (e.g. `usage`), or `GenerateContentConfig` drops `seed`/`temperature`, or pydantic
drops `model_copy(update=...)` → tokens read as zero (cost=$0) and/or determinism silently
reverts to the library's `temperature=0.7`. No exception — `getattr(response, "usage_metadata",
None)` just returns `None`, and `force_deterministic_config` returns the config as-is.

**Guarded by:** `google-genai==2.6.0` (the SDK class/field shapes) **and** `pydantic==2.13.3`
(`model_copy(update=...)` + `model_fields` on the class). Smoke-check asserts
`Models.generate_content` is callable, and that `GenerateContentConfig` has `temperature`+`seed`
and `GenerateContentResponseUsageMetadata` has the three token-count fields.

---

## (iii) OpenRouter `requests.post` wrap

**What the library does/discards.** The OpenRouter provider posts via plain HTTP. At
`providers/openrouter.py:61` it passes `post_func=requests.post` into the shared
`call_openai_compatible_api` helper. That helper (`providers/common.py`) sends the body as
`data=json.dumps(body)` — **not** `json=` — because `request_json=False` is the default
(`common.py:448`, branch at `common.py:263-270`). `build_chat_completions_body` (`common.py:85`)
sends `temperature` only when it is not `None`, and the default is `None`
(`common.py:89,102-103`) — so no temperature and no seed go out, each vendor samples at its own
default. 2.0.7 will read a `cost`/`usage` block if present (`common.py:39-41`) and forwards a
`usage_sink`, but — as with Gemini — `usage_sink` is not surfaced through
`interactive_consensus_annotation`, and the provider never asks OpenRouter for cost
(`usage:{include:true}`), so the wrap still injects both.

**The seam we bind to.** `OpenRouterCapture` (`core/capture.py`) monkeypatches
`mllmcelltype.providers.openrouter.requests.post` — which *is* the global `requests.post`,
because the provider evaluates `requests.post` from its module-level `import requests` at call
time. The wrapper:

1. **Injects determinism into the outgoing body.** It inspects `kwargs["data"]` (the
   `json.dumps(body)` string), `json.loads` it, and if it looks like a chat-completions body
   (`"messages" in body`), sets `body["temperature"] = 0`, `body.setdefault("seed", 0)`, and
   `body.setdefault("usage", {"include": True})` (the last asks OpenRouter to return cost),
   then re-`json.dumps` it back into `kwargs["data"]`. It also handles a `json=` kwarg shape
   defensively, in case the library ever flips `request_json`.
2. **Reads usage + native USD off the returned response.** After the real `post`, it re-parses
   `response.json()` (cheap; the body is already buffered, so the library's own `.json()` still
   works) and accumulates `{model: {prompt, output, total, calls, cost_usd}}`, preferring
   `usage.cost`, falling back to top-level `cost`.

The wrapper **never raises** — any failure to inject or read degrades to a normal call.
`self.usage` is the Gemini shape **plus** a native `cost_usd` that `report_cost` prefers when
`> 0` (so the rate table in `core/cost.py` mostly matters only for the direct-Gemini path).

**Silent breakage.** If the provider switches to `request_json=True` (body as `json=`), or the
parser starts surfacing usage itself, or `mllmcelltype.providers.openrouter` stops importing
`requests` as a module-level name, or the response schema's `usage`/`cost` keys move → tokens +
native cost are lost and `temperature=0` is not injected. The `smoke_check` "`…openrouter.requests
IS requests`" assertion is the tripwire for the import-shape half; the body-shape half is only
caught by `test_capture.py`'s fake-`post` assertions.

**Guarded by:** `requests==2.33.1` (the `requests.post` identity the patch swaps) and
`mllmcelltype==2.0.7` (the provider's import shape + body/parse contract). Smoke-check asserts
`mllmcelltype.providers.openrouter.requests is requests`.

---

## (iv) Logger DEBUG survival

**What the library does/discards.** `mllmcelltype` logs raw model responses **only at DEBUG**.
But `interactive_consensus_annotation` calls `annotate_clusters` once **per model**, and each
of those calls `setup_logging(...)` against the same default log dir. With an already-initialized
logger on the same dir, `setup_logging` takes the **"just update level" fast path**
(`logger.py:75-78`): it walks the existing handlers and calls `setLevel(INFO)` on them. That
downgrades *any* DEBUG `FileHandler` I attached before the first provider call — so by the time
the raw-response DEBUG line fires, my handler is at INFO and **zero** raw lines are persisted.
(For a *different* dir it takes a different branch that removes file handlers entirely — also
fatal to my handler, but the fast path is the one that bites in practice since the dir doesn't
change mid-run.) The logger name is hardcoded `"llmcelltype"` (`logger.py:9`) — note the single
`l`, not `llmcelltype`/`mllmcelltype`.

**The seam we bind to.** `core/debug_capture.py`:
- `configure_llm_logging(lens, log_dir)` attaches a DEBUG `FileHandler` (tagged
  `mllmct_debug::<lens>`) to `logging.getLogger("llmcelltype")` and forces logger+handler to
  DEBUG. Idempotent per (lens, file) — it dedupes by `baseFilename`.
- `llm_debug_capture(lens, log_dir)` is the ctx manager that makes it *survive*. It wraps
  `mllmcelltype.logger.setup_logging` (and every module that imported the name *by value* —
  `mllmcelltype.annotate` does `from .logger import setup_logging`, and the package re-exports
  it). After each library `setup_logging` runs, the wrapper **re-asserts**: re-adds the handler
  if dropped and forces DEBUG back on. Because the library re-inits once per model, the
  re-assertion fires once per model too. Originals are restored on `__exit__`.

**Silent breakage.** If the fast path stops calling `setLevel` on existing handlers, or
`setup_logging` is renamed, or the logger name changes, or `annotate` stops importing the name
by value (so my by-value patch list misses the real callsite) → DEBUG survives anyway *or*
(worse) my re-assert no longer fires and `llmcelltype_debug.log` is empty. Either way: no error.
This is the only seam with **no dedicated smoke-check structural assertion beyond
"`setup_logging` is callable"** — the per-model re-init behaviour is behavioural, not
structural, so trust `test_*` + a real dry-run here.

**Guarded by:** `mllmcelltype==2.0.7`. Smoke-check asserts `mllmcelltype.logger.setup_logging`
is callable (a weak guard — see above).

---

## Regeneration procedure (run after any `uv lock`)

The pins are the contract. **Never** regenerate the lock without re-running the structural gate:

```bash
SKILL=skills/mllmcelltype-consensus-annotation
uv lock --project "$SKILL"                                  # regenerate lock
uv run --project "$SKILL" python "$SKILL/checks/smoke_check_versions.py"
# equivalently:  uv run --project "$SKILL" mllmct check-env
```

`smoke_check_versions.py` is the tripwire. It asserts, **offline, no API key, no AnnData**:

1. the four exact versions in `mllmct._version.PINNED`
   (`mllmcelltype 2.0.7 / google-genai 2.6.0 / pydantic 2.13.3 / requests 2.33.1`); and
2. that every structural seam still *exists*: `create_prompt` and
   `interactive_consensus_annotation` each accept a `prompt_template` parameter (the native seam
   that replaced the retired global swap), `logger.setup_logging` callable,
   `providers.openrouter.requests IS requests`, `Models.generate_content` callable, and
   `GenerateContentConfig`/`usage_metadata` carry the expected fields.

**If a seam check FAILs:** do **not** trust the wrappers. The run would still "succeed" and lose
tokens/determinism/prompt/DEBUG with no error. Reconcile **this file** against the new pinned
source first — re-verify each file:line, fix the wrapper in `core/`, re-run `pytest tests/`,
*then* update `mllmct._version.PINNED` and the `pyproject.toml` pin. Only after the smoke-check
goes green is the new lock safe to commit.

### The four pins and what each guards

| Pin | Guards |
|---|---|
| `mllmcelltype==2.0.7` | seam (iv) the per-model `setup_logging` fast-path behaviour + logger name; the OpenRouter provider's `data=json.dumps` body shape + no-temperature body in seam (iii); and the native `prompt_template=` parameter that retired seam (i) |
| `google-genai==2.6.0` | seam (ii): `Models.generate_content` location/signature, `usage_metadata` field names, `GenerateContentConfig.{temperature,seed}` |
| `pydantic==2.13.3` | seam (ii) determinism path: `model_copy(update=...)` semantics + `model_fields` on the class |
| `requests==2.33.1` | seam (iii): the `requests.post` callable identity the OpenRouter patch swaps in place |

`pandas` / `pyyaml` / `python-dotenv` are floors (`>=`), not pins — nothing monkeypatches them.
