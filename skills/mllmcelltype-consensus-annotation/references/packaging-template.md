# Packaging notes — this skill as reference implementation

The "thin interface over a deep, version-locked, tested module" philosophy is now a toolkit-level
spec: **[`../../../docs/packaged-skills.md`](../../../docs/packaged-skills.md)**. That's where the
contract, the package-vs-docs-only decision, the distribution stance, and the copy-this checklist
live. **This skill is its reference implementation** — to build a new packaged skill, copy this dir
and follow the checklist there.

What follows is only the handful of things specific to *this* skill that the toolkit spec doesn't
cover. For the seam/pin internals see [`monkeypatch-internals.md`](monkeypatch-internals.md); for
the cell-state plugin see [`cell-state-annotation.md`](cell-state-annotation.md).

---

## Skill-specific decisions

- **The four `==` pins** — `mllmcelltype`, `google-genai`, `pydantic`, `requests` — are pinned exact
  because mllmct monkeypatches their *internals* (attribute names + call signatures, not public
  APIs). Everything else (`pandas`, `pyyaml`, `python-dotenv`) is a `>=` floor. Which is which, and
  why, is in [`monkeypatch-internals.md`](monkeypatch-internals.md).
- **`checks/smoke_check_versions.py` asserts these four seams specifically**: the prompts template
  global, the consensus entrypoint + logger, the OpenRouter `requests.post` identity, and the
  google-genai config/usage fields. The version-of-truth is `mllmct._version.PINNED`. Never `uv lock`
  without re-running the smoke-check.
- **The project-specific seam is the `EvidenceProvider` ABC** (+ a reference plugin) and the YAML
  `profiles/` — that's how cell-state biology is injected *without* importing any analysis-repo
  `config.py`, keeping the package extractable. Details in
  [`cell-state-annotation.md`](cell-state-annotation.md).
- **End-to-end test** patches the consensus entry point *as imported by the engine module* (a
  function-level `FakeConsensus`), not at the package root — so a `from lib import name` local
  binding can't escape the patch. Assertions key on the minimum dict keys the engine consumes.
