# 09 — Validator tooling patterns for ADR-007

ADR-007 specifies `sciagent doctor` with five check classes:
1. SKILL.md YAML/markdown parse
2. `metadata.sciagent.*` namespace conformance (ADR-001)
3. `requires:` graph: resolution, cycles, orchestrator/atomic invariants (ADR-002, 003)
4. agentskills.io reference validator delegation (skills-ref)
5. Symlink hygiene + manifest consistency

This file maps each check to existing tooling.

## YAML frontmatter validation

Three real choices:

- **[Yamale](https://github.com/23andMe/Yamale)** — schema-in-YAML for YAML. Pythonic, lightweight. The schema language itself is YAML, which makes it easy to read but harder to integrate with IDE tooling. Good fit if sciagent stays Python-side for validation.
- **JSON Schema** (via `jsonschema` Python package or [check-jsonschema CLI](https://check-jsonschema.readthedocs.io/)) — industry standard, broad tooling support (IDEs, CI), can validate YAML by converting to JSON first. `check-jsonschema` is a single-binary CLI; works against `.yaml` files natively.
- **Pydantic + ruamel.yaml** — programmatic; gives you typed Python objects post-validation. Heavier dependency but worth it if `lib/sciagent/roles.sh` ever gains a Python sibling.

**Recommendation:** JSON Schema via `check-jsonschema`. Reasons:
- IDE editors (VS Code, JetBrains) already render JSON Schema validation against YAML in-editor — your skill authors get autocomplete and red squiggles for free.
- CI integration is one-line in any CI system.
- The schema itself is portable — if anyone wants to reuse sciagent's SKILL.md extensions in another tool, they consume the schema.
- Pure-Python via `jsonschema` if `check-jsonschema` is not preferred.

The schema would express ADR-001 (metadata.sciagent.* namespace) and ADR-003 (scope ∈ {orchestrator, atomic}) as a JSON Schema document checked into `schemas/sciagent-skill.schema.json`. ADR-002's `requires:` graph is *structural* (graph constraints over multiple files) and cannot be expressed in JSON Schema — that's a separate check.

## DAG / cycle detection for `requires:` graph (ADR-002)

Trivial. Python `graphlib.TopologicalSorter` (stdlib since Python 3.9) does this in 10 lines:

```python
import graphlib
ts = graphlib.TopologicalSorter()
for skill, deps in requires_graph.items():
    ts.add(skill, *deps)
try:
    list(ts.static_order())  # raises CycleError on cycle
except graphlib.CycleError as e:
    print(f"Cycle detected: {e.args[1]}")
```

If sciagent stays in bash (matches the existing `lib/sciagent/*.sh` style), a topological sort in bash is ~30 lines. The existing repo already includes `lib/sciagent/skill_deps.sh` per the `bin/sciagent` source — depth-first dep resolution is present. Cycle detection can be added there without adding Python.

Resolution check (every name in `requires:` corresponds to a directory in `skills/`) is a `for skill in requires; do [[ -d skills/$skill ]] || error; done`. Trivial bash.

Orchestrator/atomic invariants (ADR-003: orchestrators must have non-empty requires, atomics must have empty requires) are similarly trivial per-file checks.

## Symlink hygiene

Patterns from real package managers:

- **Homebrew / Nix** detect broken/dangling symlinks via `find -L <dir> -type l -! -exec test -e {} \;` (POSIX) or `readlink -e` (returns nothing for dangling). Both are bash-native.
- **pnpm** maintains a `.pnpm-store/` directory and validates symlinks reference items in it.

For sciagent:
- `.sciagent/manifest.json` already records which symlinks sciagent created (per architecture doc §9). `doctor` walks that list, calls `readlink -e` on each, reports dangling or missing entries.
- Walking `.claude/skills/`, `.claude/agents/`, `.claude/commands/`, `.agents/skills/`, etc. for symlinks not in the manifest catches "orphaned" symlinks from prior aborted activations.

This is ~50 lines of bash. No new dependencies.

## Manifest consistency

Cross-check `.sciagent/manifest.json` against actual filesystem state:
- Every symlink in manifest exists and points where it claims to.
- No "extra" symlinks in the dirs sciagent owns that aren't in the manifest.
- Stack in manifest matches stack in AGENTS.md managed block (block hash check).

The architecture doc already specifies the manifest format. Validation is a straightforward diff between two sets.

## Delegation to agentskills.io skills-ref

`uvx --from git+https://github.com/agentskills/agentskills#subdirectory=skills-ref skills-ref validate <skill-dir>` — runs the upstream validator on each skill.

Per file 05, skills-ref **does not** check dependency graphs or cross-skill invariants. Its scope is single-skill structural validation (frontmatter fields, naming, token budgets). So:

- `sciagent doctor` runs skills-ref per skill as the first-pass linter.
- ADR-002/003 invariants are sciagent-owned implementations.

Subtle issue: skills-ref may flag the structured-value-under-`metadata.sciagent.*` pattern from ADR-001 if it strictly enforces the spec's "string → string" type. Test on a sample skill before assuming smooth delegation. See file 05's caveat.

## Lockfile pattern (for future ADR, not 007)

`requires: [{name: harmonypy, version: ">=0.1.0"}]` (ADR-002) implies version resolution. Current state of the art for skill graphs:

- **Cargo (Rust)** — `Cargo.lock` is a flat list of resolved package@version pairs.
- **pnpm** — `pnpm-lock.yaml` records the full resolved DAG including transitive deps and their resolved versions.
- **uv** — `uv.lock` (file 07) is a similar flat resolved list with hash-pinning.

For sciagent: a `.sciagent/skills.lock` would record the resolved version of every skill in the active stack's transitive closure. Generated on `activate`, used by `doctor` to detect post-activate drift. Out of scope for ADR-007 but a natural ADR-007.5.

## Build vs. wrap upstream

The recommendation is hybrid:
- **Wrap skills-ref** for single-skill structural validation.
- **Build sciagent-native** for everything cross-skill (graph, scope invariants, symlinks, manifest).

Approximate effort: ~300 LOC bash + ~50 LOC bash for the schema check (or ~150 LOC Python if going JSON Schema route). Matches the spec's own estimate.

## One small thing the spec misses

`--fix` mode for `sciagent doctor` (auto-repair) is dangerous for graph invariants. Re-creating a dangling symlink is safe. "Fixing" a `requires:` cycle is not. ADR-007 should narrow `--fix` to symlinks + manifest entries only — never edits SKILL.md frontmatter. Worth grilling.

Sources:
- [Yamale — 23andMe](https://github.com/23andMe/Yamale)
- [JSON Schema Everywhere — YAML validation](https://json-schema-everywhere.github.io/yaml)
- [check-jsonschema CLI](https://check-jsonschema.readthedocs.io/)
- [Python graphlib — TopologicalSorter](https://docs.python.org/3/library/graphlib.html)
- [agentskills/skills-ref](https://github.com/agentskills/agentskills/tree/main/skills-ref)
- [Cargo book — Cargo.lock](https://doc.rust-lang.org/cargo/guide/cargo-toml-vs-cargo-lock.html)
- [pnpm lockfile structure](https://pnpm.io/git#lockfiles)
