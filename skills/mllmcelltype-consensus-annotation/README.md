# mllmct

Version-locked CLI for **mLLMCelltype** multi-LLM consensus annotation of scRNA-seq
clusters — **cell-type** *and* **cell-state** — with token/cost capture, forced
determinism (`temp=0, seed=0`), **Python-recomputed** consensus metrics, and a full
reproducibility trace. The same tool does both modes; you switch by passing a
different YAML **profile**.

This package is the *implementation*. If you just want to **use** the tool, read
[`SKILL.md`](SKILL.md) — it is the thin interface contract (what it is, its arguments,
how to run it). You do not need to read the source to use it.

## Bootstrap (one time)

```bash
SKILL=skills/mllmcelltype-consensus-annotation     # adjust to your path
uv sync --project "$SKILL"                          # build the locked sandbox (.venv) from uv.lock
uv run --project "$SKILL" python "$SKILL/checks/smoke_check_versions.py"   # gate: pins + patch seams intact
uv run --project "$SKILL" mllmct selftest           # offline core-logic smoke (no API)
```

`uv` (>=0.11) is the only prerequisite. The sandbox is isolated — it does **not**
inherit the system or `/opt/venvs/base` site-packages.

## Run

```bash
# stable launcher (no venv activation needed):
"$SKILL/bin/mllmct" annotate --profile "$SKILL/profiles/celltype.yaml" \
  --markers markers.csv --species human --tissue "peripheral blood" \
  --models 'google/gemini-2.5-flash' --lens leiden --out runs/demo --env-file .env
```

See `mllmct <subcmd> --help` and [`SKILL.md`](SKILL.md) for the full surface.

## Test

```bash
uv run --project "$SKILL" pytest "$SKILL/tests" -q
```

Tests are **offline** (no API, no AnnData) and run against the pinned `mllmcelltype`
so the version-sensitive monkeypatches are exercised. Fixtures under `tests/fixtures/`
are **synthetic and de-identified** — see `tests/fixtures/README.md`.

## How it works / extend it

Internals are documented in `references/` (not needed to use the tool):
`monkeypatch-internals.md`, `cell-state-annotation.md`, `packaging-template.md`.

## Reproducibility note

The four `==` pins (`mllmcelltype`, `google-genai`, `pydantic`, `requests`) are
load-bearing: the tool wraps version-sensitive internals of these packages. Only
regenerate the lock with `uv lock` followed by `checks/smoke_check_versions.py`.
