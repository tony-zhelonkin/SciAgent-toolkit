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

---

## 6. `compatibility:` — declaring what a skill needs from outside itself

Applies to **every** skill, packaged or docs-only. Most skills need nothing outside their own
directory and must **omit the key entirely**; declaring a dependency a skill does not have is as
wrong as omitting one it does.

Add it only when the skill would misbehave — not merely read a little oddly — outside a SciAgent
analysis repo. Place it after `license:`. It is a **quoted** YAML scalar (the value contains
`: `, which is not a legal plain scalar) of at most **500 characters**
(`skills/skill-creator/scripts/quick_validate.py:86-92`).

```yaml
---
name: figure-style
description: ...
license: MIT
compatibility: "sciagent-toolkit: figure-style; sciagent-scaffold: 02_analysis/helpers/figure_style.R, 02_analysis/config/analysis_config.yaml, 03_results/"
---
```

**Grammar** — clauses joined by `"; "`, items within a clause by `", "`:

```
compatibility := clause ( "; " clause )*
clause        := flavour ": " item ( ", " item )*
flavour       := "sciagent-scaffold" | "sciagent-toolkit" | "sibling-skill" | "external-module"
item          := non-empty token; no ";", no ",", no "#";
                 no leading or trailing whitespace
```

**The four flavours** — the set is closed; anything else is a hard validate failure, because a
misspelled flavour reads correct and enforces nothing:

| Flavour | Item is | Checked by `sciagent validate` |
|---|---|---|
| `sciagent-scaffold` | a **repo-root-relative path** into the analysis-repo layout (`02_analysis/config/analysis_config.yaml`, `03_results/objects/`). A trailing `/` means "directory". | Shape only — must be relative, no `..`. The consumer project is absent at validate time and is deliberately never stat()ed: `activate` runs `validate --quiet` as a pre-flight, so a project-state check here could hard-block activation. |
| `sciagent-toolkit` | a **bare directory name** under the toolkit's `lib/`. Only `figure-style` and `interactive-style` exist — those are the two `symlink_create_helper_lib` mounts into `02_analysis/helpers/`. | The directory must exist in this checkout. |
| `sibling-skill` | a **bare directory name** under `skills/`. Use when the skill sources another skill's `scripts/` or requires its `references/`. | The skill directory must exist in this checkout. |
| `external-module` | free-form: a submodule or package the toolkit does **not** provide — `RNAseq-toolkit`, `TE-RNAseq-toolkit`, `pathway-explorer`. | Nothing. Unverifiable by construction; that is the point of the flavour. |

Do **not** file a `01_modules/<other-toolkit>` dependency under `sciagent-scaffold`: SciAgent does
not ship it, and the declaration would be a lie. That is what `external-module` is for.

**Do not use `metadata.requires:`.** That key belonged to a removed taxonomy resolver;
`skills/skill-creator/scripts/add_requires_field.py` is an inert tombstone that refuses to run. See
`skills/README.md`.

Enforcement lives in `lib/sciagent/validate.sh` (`_validate_compatibility`); every rule above is
mutation-tested in `tests/test_validate_compatibility.sh`, and the shipped declarations are pinned
by `tests/test_skill_compatibility_declared.sh`. The audit that produced the current set is
`docs/proposals/2026-08-11-offline-distribution/50_ADRs.md` ADR-D6.
