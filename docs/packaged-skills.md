# Packaged skills — design spec

A toolkit-level tier for skills that ship real, tested, version-locked executable code rather than
prose. The reference implementation is `skills/mllmcelltype-consensus-annotation/` — and the
governing principle is **the code is the spec**: there is no `_TEMPLATE_PACKAGED/`, you copy the
reference skill. This doc defines the tier; that skill *is* the worked example.

---

## 1. The idea in one breath

A docs-only skill is a SKILL.md an agent *reads and interprets* — so the agent re-derives the
procedure each time, non-deterministically. A **packaged skill** is a thin SKILL.md *interface*
over a deep, version-locked, tested module that an agent or human **outsources execution to** by
calling one CLI. The implementation never enters the caller's context window; the behaviour is the
same on every run. Packaging trades flexibility for **determinism + reuse**: the skill becomes a
black box at the call site and a glass box for whoever maintains it.

## 2. Package it, or stay docs-only?

Most skills should stay flat SKILL.md files — probably no need to over-build. But there`s value in 
**Packaging a skill when it ships executable logic that must run identically every time**, specifically when any of these hold:

- it wraps a **version-sensitive / monkeypatched** dependency (a loose bump would silently break a
  guarantee with no error);
- it wants a **tested CLI** with an enforceable output contract, not a recipe an agent re-types;
- determinism, captured cost/trace, or reproducibility is itself part of the deliverable.

If none apply — it's guidance, conventions, a decision tree, a thin wrapper over a stable public API
— a flat SKILL.md is correct. Packaging has carrying cost (a lock to maintain, a smoke-check to
keep honest); pay it only when a skill earns it.

## 3. The contract

A packaged skill MUST satisfy all of the following. The reference column points at the file in
`skills/mllmcelltype-consensus-annotation/` that demonstrates each — read those, don't re-derive them.

| # | Requirement | Reference |
|---|---|---|
| 1 | **Self-contained `uv` project** inside the skill dir: `pyproject.toml` + a **committed `uv.lock`**. `[tool.uv]` keeps the sandbox from inheriting system / base site-packages. | `pyproject.toml` |
| 2 | **EXACT `==` pins for the fragile deps** (the ones whose internals you reach into / monkeypatch), **`>=` floors for the rest** (deps you use through a stable public API). Pin only what you reach into. | `pyproject.toml` |
| 3 | **Console script** (`[project.scripts]`) + a **symlink-safe `bin/<tool>` launcher** that dereferences its own `BASH_SOURCE` before resolving the skill dir, so it survives being symlinked into `.claude/skills/` by `sciagent activate`. Fails loudly when `uv`/lock are absent. | `bin/mllmct` |
| 4 | **A smoke-check that GATES lock regeneration**: offline, asserts the exact pins AND that every structural seam the code binds to still exists. Rule: never `uv lock` without re-running it. Surface it via a `check-env` subcommand too. | `checks/smoke_check_versions.py` |
| 5 | **Offline tests + synthetic de-identified fixtures + a `selftest` subcommand.** Tests run inside the locked sandbox (so version-sensitive imports resolve against the pins). Fixtures carry zero real biology/identities. `selftest` is dependency-free (works right after `uv sync`, no `test` extra). | `tests/`, `tests/fixtures/`, `selftest` cmd |
| 6 | **Thin SKILL.md (interface only)**; internals live in `src/` + `references/`. SKILL.md is agent-facing (what/args/how-to-call); README is the human bootstrap doc; deep internals are `references/*.md`, read only by a maintainer. | `SKILL.md`, `src/`, `references/` |
| 7 | **STRICTLY self-contained**: no imports from other skills or from the host toolkit/project (no `config.py`). Project-specific bits go behind an interface (the reference uses an `EvidenceProvider` ABC + a plugin) parameterized by config you pass in. This keeps the package **extractable**. | `src/mllmct/` |
| 8 | **`tests/run_skill_tests.sh`** that skips gracefully (exit 0) when uv/lock are absent; the toolkit's `tests/run-all.sh` discovers it automatically via its `skills/*/tests/run_skill_tests.sh` sweep. | `tests/run_skill_tests.sh` |
| 9 | **A one-line tier note at the top of the SKILL.md body** — "Packaged skill: execution is outsourced to the pinned CLI below." There is no frontmatter marker; the presence of `pyproject.toml` + a lockfile + `src/` is the actual, checkable signal. | `SKILL.md` body |

## 4. Distribution stance

Opinionated, in order:

- **Monorepo now.** Packaged skills live in `skills/<name>/` like every other skill. No premature
  splitting.
- **Keep them extractable** (contract #7). A packaged skill should be liftable to its own repo with a
  `git mv` and zero import rewiring.
- **Promote to PyPI / its own repo ONLY when a skill earns it**: reused across **≥3 unrelated
  projects**, OR an external consumer wants it standalone, OR its **release cadence diverges** from the
  toolkit's. Until then, the carrying cost of a separate repo isn't justified.
- **NEVER nested git submodules.** The toolkit is *already* a submodule of analysis projects; nesting a
  packaged skill as a sub-submodule means a three-level pointer-bump every change — don't.
- **Separate interface from implementation.** When a skill is promoted, the **SKILL.md interface
  contract stays discoverable** in the toolkit (so activation + routing still see it); the
  **implementation package can live anywhere** (PyPI, its own repo). The launcher resolves the tool;
  the toolkit need not vendor the code.

## 5. Copy-this checklist (new packaged skill)

The reference skill IS the template. To stamp out a new one:

1. **Copy the reference skill dir** `skills/mllmcelltype-consensus-annotation/` to `skills/<name>/`.
2. **Strip the domain bits**: the `EvidenceProvider` plugin, the `profiles/`, the fixtures, and any
   mllmct-specific `references/*.md`. Keep the *shape* — `pyproject.toml`, `bin/`, `checks/`,
   `tests/`, `src/<pkg>/{__init__,_version,cli}.py`, `tests/run_skill_tests.sh`.
3. **Rename** the package (`src/<pkg>/`), the `[project.name]`, the `[project.scripts]` console-script,
   the `bin/<tool>` launcher, and `[tool.hatch.build.targets.wheel] packages`.
4. **Re-pin deps**: `==` for the fragile ones, `>=` floors for the rest; put the expected versions in
   `<pkg>._version.PINNED`. Run `uv lock`.
5. **Rewrite `checks/smoke_check_versions.py`** to assert YOUR pins + YOUR structural seams, then run
   it (`uv run … checks/smoke_check_versions.py` must exit 0). Lock + smoke-check is one ritual.
6. **Write your tests + synthetic fixtures + `fixtures/README.md`**; an end-to-end test that patches
   the external entry point *as imported by your module*. Wire a `selftest` + `check-env` subcommand.
7. **Write a thin SKILL.md** (interface) opening with the packaged-skill tier note, a README
   (human bootstrap), and `references/*.md` for internals. Validate with
   `python skills/skill-creator/scripts/quick_validate.py skills/<name>/`.

The reference skill's files are the spec for each step above — when in doubt, diff against them.
