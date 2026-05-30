# Packaging template — "locked tool + CLI + logging + tests"

> Notes on the *pattern* this skill embodies, so the next SciAgent-toolkit skill that needs to
> wrap a fragile third-party library (or just ship real, tested executable code) can copy it
> instead of reinventing it. `mllmct` is the reference implementation; this doc is the recipe.

The shape is: **a self-contained `uv` project living inside the skill directory, exposing one
console-script CLI through a symlink-safe `bin/` wrapper, with a committed `uv.lock`, an offline
test suite, and a smoke-check that gates lock regeneration.** The SKILL.md stays thin; all the
weight is hidden in `core/` + `references/`.

---

## Why a `uv` project inside a skill

- **Isolated sandbox.** `uv sync` builds `<skill>/.venv` from the lock and — critically —
  does **not** inherit the system or `/opt/venvs/base` site-packages (`[tool.uv]` note in
  `pyproject.toml`). So the skill's dependency graph can't be poisoned by, or poison, the
  host environment. The agent's main Python stays untouched.
- **Committed `uv.lock` = exact, reproducible graph.** Every transitive is pinned. Anyone who
  `uv sync`s gets byte-identical deps, which is the whole point when the tool wraps
  version-sensitive internals.
- **Exact pins for the fragile deps, floors for the rest.** The deps whose *internals* are
  monkeypatched are pinned `==` (here: `mllmcelltype`, `google-genai`, `pydantic`, `requests`
  — see `monkeypatch-internals.md`); everything the skill merely *uses* through a public API
  (`pandas`, `pyyaml`, `python-dotenv`) is a `>=` floor. Pin only what you reach into.

The rejected alternative — a flat `scripts/` package reached via `PYTHONPATH` — can't give you a
real locked sandbox or a clean console-script, so it loses the reproducibility guarantee.

---

## Console-script + symlink-safe `bin/` wrapper

Two layers:

1. **Console script** declared in `pyproject.toml`:
   ```toml
   [project.scripts]
   mllmct = "mllmct.cli:main"
   ```
   `uv run --project <skill> mllmct …` then works from anywhere.

2. **A `bin/<tool>` bash wrapper** so a human or agent can call it *without* knowing about uv or
   activating a venv — and so it survives being **symlinked** into `.claude/skills/` by
   `sciagent activate`. The wrapper dereferences its own `BASH_SOURCE` through any symlink chain
   before resolving the skill dir, then `exec uv run --project "$SKILL_DIR" <tool> "$@"`:
   ```bash
   src="${BASH_SOURCE[0]}"
   while [ -L "$src" ]; do
     dir="$(cd -P "$(dirname "$src")" && pwd)"; src="$(readlink "$src")"
     [[ "$src" != /* ]] && src="$dir/$src"
   done
   SKILL_DIR="$(cd -P "$(dirname "$src")/.." && pwd)"
   ```
   It also fails *loudly and actionably* when `uv` is missing (exit 127) or the lock is absent
   (exit 1, with the `uv sync` command to run). Naively using `$(dirname "$0")` would break the
   moment the wrapper is symlinked — that's the bug this guards against.

---

## Thin SKILL.md over a deep, tested module

The governing philosophy: **SKILL.md is the namespace interface; the internals are hidden.**

- SKILL.md says *what the tool is, its arguments, and how to invoke it* — nothing about how the
  monkeypatches work. A human or agent reads it, calls `bin/<tool> …`, and **outsources
  execution**: the implementation never enters their context window.
- The real logic lives in `src/<pkg>/core/` (domain-agnostic, reusable) plus a thin
  `cli.py`/`engine.py` that only sequences it. Deep "how it works / how to extend" material
  goes in `references/*.md`, read **only** by someone who wants internals.
- The README is the *human bootstrap* doc (uv sync, run, test, extend); SKILL.md is the
  *agent-facing* doc. They're different audiences — keep them separate.

The payoff: the tool is a black box at the call site, but a glass box for whoever maintains it.

---

## The smoke-check that gates lock regeneration

`checks/smoke_check_versions.py` is the tripwire that makes the `==` pins safe to live with. It
asserts, **offline**: (a) the exact pinned versions, and (b) that every structural seam the
monkeypatches bind to still *exists*. The rule is: **never `uv lock` without re-running it.**

```bash
uv lock --project "$SKILL"
uv run --project "$SKILL" python "$SKILL/checks/smoke_check_versions.py"   # must exit 0
```

A FAIL means a dependency bump silently moved a seam — the tool would still "work" but lose its
guarantees with no error. The version-of-truth is `<pkg>._version.PINNED` (kept in sync with the
`pyproject.toml` pins), so the check has one place to read the expected versions from. This is
the generic move: **if you pin internals, also ship a structural assertion that the internals
are where you think they are, and wire it into lock-regen + a `check-env` subcommand.**

---

## Offline tests + synthetic fixtures + a `selftest` subcommand

Three layers of self-verification, all **offline (no API, no AnnData, no network)**:

1. **A `pytest` suite** (`tests/`) run *inside the locked venv* so the version-sensitive imports
   resolve against the pinned library:
   ```bash
   uv run --project "$SKILL" --extra test pytest "$SKILL/tests" -q
   ```
   The end-to-end test uses a **function-level monkeypatch** (`FakeConsensus`) that patches the
   library entry point *as imported by the engine module* — not at the package root — so a
   `from lib import name` local binding can't escape the patch. Assertions key on the **minimum**
   set of dict keys the engine actually consumes, so a benign upstream schema addition doesn't
   break them.
2. **Synthetic, de-identified fixtures** (`tests/fixtures/`): hand-authored CSV/JSON chosen to
   hit every render/branch, carrying **zero** real biology — no identities, conditions, donors,
   barcodes, accessions, or project paths. A `fixtures/README.md` documents exactly what was
   kept out and why. (Public gene *symbols* are fine — they're not PII and are load-bearing for
   prompt rendering.)
3. **A `selftest` subcommand** (`<tool> selftest`) that runs a dependency-free invariant check
   *without* pytest, so it works immediately after `uv sync` (which installs runtime deps only,
   not the `test` extra). It exercises the core invariants (here: consensus metrics, harmonize,
   guards, determinism rewrite, token-capture non-clobber, cost, reconcile branches) and exits
   0/1 — a true post-bootstrap health check.

---

## Skip-gracefully toolkit wiring

The skill ships its own `tests/run_skill_tests.sh`, and the toolkit's `tests/run-all.sh` does a
single discovery sweep over `../skills/*/tests/run_skill_tests.sh`. The shim **skips gracefully**
(prints `SKIP …`, exits 0) when `uv` or the lock are absent, so the toolkit's bash CI never
blocks on an un-bootstrapped machine:

```bash
command -v uv >/dev/null 2>&1 || { echo "SKIP [tool]: uv not on PATH"; exit 0; }
[ -f "$SKILL/uv.lock" ] || { echo "SKIP [tool]: uv.lock absent"; exit 0; }
env -u VIRTUAL_ENV uv run --project "$SKILL" python "$SKILL/checks/smoke_check_versions.py"
env -u VIRTUAL_ENV uv run --project "$SKILL" --extra test pytest "$SKILL/tests" -q
```

(`env -u VIRTUAL_ENV` drops any inherited venv so uv targets *this* skill's sandbox cleanly.)
The toolkit side already exists — `run-all.sh` loops the sweep and folds results into its
pass/fail summary — so a new skill only has to drop in its own `run_skill_tests.sh`.

---

## To copy this for a new skill — checklist

1. **Scaffold the `uv` project** under `skills/<name>/`: `pyproject.toml`
   (`requires-python`, deps, `[project.scripts]`, `[tool.hatch.build.targets.wheel]
   packages=["src/<pkg>"]`, `[tool.uv]`), `src/<pkg>/{__init__,_version,cli}.py`.
2. **Pin only the fragile deps `==`**, floors for the rest. Put the expected versions in
   `<pkg>._version.PINNED`.
3. **Write `checks/smoke_check_versions.py`** asserting the pins AND each structural seam your
   code reaches into. Run `uv lock`, then verify the smoke-check exits 0.
4. **Add the `bin/<tool>` symlink-safe wrapper** (dereference `BASH_SOURCE`, exec
   `uv run --project`); fail loudly if uv/lock missing.
5. **Keep `core/` domain-agnostic**; put any project-specific piece behind an interface (here:
   the `EvidenceProvider` ABC + a reference plugin), parameterized by config you pass in — never
   by importing a project `config.py`.
6. **Offline tests + synthetic fixtures + `fixtures/README.md`**; an end-to-end test that
   patches the external entry point as imported by *your* module.
7. **A dependency-free `selftest` subcommand** + a `check-env` that delegates to the smoke-check.
8. **Drop in `tests/run_skill_tests.sh`** that skips gracefully and runs smoke-check + pytest in
   the locked sandbox. The toolkit's `run-all.sh` discovers it automatically.
9. **Thin SKILL.md** (agent-facing interface) + a **README** (human bootstrap/run/test/extend) +
   `references/*.md` for the deep internals. Validate with
   `python skills/skill-creator/scripts/quick_validate.py skills/<name>/`.
10. **Document the lock-regen ritual** prominently: `uv lock` is only ever followed by the
    smoke-check before the new lock is trusted/committed.
