# 05 — agentskills.io spec, metadata semantics, and the reference validator

Grounds ADRs 001 (`metadata.sciagent.*` namespace), 002 (`requires` graph), 003 (orchestrator/atomic scope), 007 (`sciagent doctor`).

## Where the spec lives

The canonical spec is at [agentskills.io/specification](https://agentskills.io/specification). Spec source is at [github.com/agentskills/agentskills](https://github.com/agentskills/agentskills) with a [skills-ref/](https://github.com/agentskills/agentskills/tree/main/skills-ref) reference implementation.

Adopters: Anthropic (Claude Code), Pi coding agent, Codex CLI. The format is convergent across major harnesses; there is no real competitor.

## Frontmatter schema (verbatim from spec)

| field | required | constraints |
|---|---|---|
| `name` | yes | ≤64 chars, `^[a-z0-9]+(-[a-z0-9]+)*$`, must match parent directory name |
| `description` | yes | ≤1024 chars, non-empty, should describe what + when |
| `license` | no | license name or reference to a bundled license file |
| `compatibility` | no | ≤500 chars, environment requirements |
| `metadata` | no | **map from string keys to string values; clients can store additional properties** |
| `allowed-tools` | no | space-separated, experimental |

Body: Markdown after the frontmatter. No format restrictions. Recommended sections: step-by-step instructions, IO examples, edge cases.

**Progressive disclosure model** (load behavior, not metadata constraint):
1. `name` + `description` (~100 tokens): loaded at startup for every skill.
2. Full `SKILL.md` body (<5000 tokens recommended): loaded when skill activates.
3. `scripts/`, `references/`, `assets/`: loaded on demand.

Spec recommends keeping `SKILL.md` under 500 lines and pushing detail into `references/`.

## The "reasonably unique key names" quote — verified

ADR-001's context paraphrases the spec's metadata-extensibility guidance. Exact text from [agentskills.io/specification](https://agentskills.io/specification):

> The optional `metadata` field:
> - A map from string keys to string values
> - Clients can use this to store additional properties not defined by the Agent Skills spec
> - **We recommend making your key names reasonably unique to avoid accidental conflicts**

ADR-001's claim is correctly attributed. The exact mitigation (namespace under `metadata.sciagent.*`) is one valid response. An alternative — using a prefix like `sciagent-requires`, `sciagent-scope` — would also satisfy the spec's "reasonably unique" guidance. The namespacing approach is more readable and is what virtually every spec-extension convention does (think `x-*` HTTP headers, `xmlns:`, etc.).

**Caveat the spec does not address:** `metadata` is typed as `string → string`. The current SciAgent-toolkit skills use **structured values** (lists, nested maps) under `metadata.requires`, `metadata.complementary-skills`, `metadata.contraindications`. This is a spec violation, mild but real — the spec defines a flat string map. ADR-001 inherits this violation by moving the structured keys under `metadata.sciagent.*`. The reference validator (skills-ref) may flag this.

Two options if it does:
- (a) Encode structured values as JSON strings (`requires: "[\"foo\", \"bar\"]"`). Ugly but spec-conformant.
- (b) Treat the spec as guidance and accept structured YAML under `metadata.sciagent.*`. Trust that skills-ref will not strictly enforce types under the extension namespace.

Worth confirming with skills-ref on a sample skill before committing to ADR-001's migration. **This is a concrete grilling question for ADR-001 that the spec does not address.**

## Versioning conventions

The spec does **not** define `version`, `last-reviewed`, `requires`, `complementary-skills`, `contraindications`, `scope`, `tier`, or any other field beyond the six in the table above. Anything else is consumer-defined metadata.

Existing sciagent skills (e.g., `coresh-signature-search/SKILL.md`) place `version`, `last-reviewed`, `category`, `tier`, `tags`, `upstream-docs`, `skill-author` directly under `metadata.*` and structured fields (`requires`, `complementary-skills`, `contraindications`) also at `metadata.*`. ADR-001 proposes splitting these — namespace-relevant ones under `metadata.sciagent.*`, generic ones (`skill-author`, `version`, `last-reviewed`, `upstream-docs`, `tags`, `category`, `license`) stay at `metadata.*`. This split mirrors the de-facto convention emerging in [openclaw/skills/agentskills-io](https://github.com/openclaw/skills/blob/main/skills/killerapp/agentskills-io/SKILL.md), though no formal convention exists.

## skills-ref reference validator

CLI commands (three total):

```
skills-ref validate <path>          # check structure + frontmatter
skills-ref read-properties <path>   # output metadata as JSON
skills-ref to-prompt <path>         # generate XML system-prompt block for an agent
```

Installation:

```
# permanent
uv tool install git+https://github.com/agentskills/agentskills#subdirectory=skills-ref
# one-shot
uvx --from git+https://github.com/agentskills/agentskills#subdirectory=skills-ref skills-ref validate ./skill
```

Scope: **single-skill structural validation only**. Checks frontmatter fields, directory naming, token budgets. Does **not**:
- Validate dependency graphs across skills
- Detect `requires:` cycles
- Cross-check that referenced skills exist
- Validate orchestrator-vs-atomic invariants
- Validate symlink hygiene

**This is a load-bearing finding for ADR-007.** `sciagent doctor` cannot delegate dependency-graph validation to skills-ref. ADR-007 must implement ADR-002/003 invariants natively. skills-ref is a useful first-pass linter for individual SKILL.md files — call it per skill — but everything cross-skill is sciagent's responsibility. See file 09 for validator-tooling patterns.

## Two Anthropic reference SKILL.md files for comparison

[anthropics/skills/skill-creator/SKILL.md](https://github.com/anthropics/skills/blob/main/skills/skill-creator/SKILL.md) — frontmatter is minimal: `name`, `description`. No `metadata` block at all. Body is the whole skill.

[anthropics/skills](https://github.com/anthropics/skills) document skills (pdf, docx, pptx, xlsx) follow the same pattern: minimal frontmatter, body-heavy.

Compare to sciagent's current `coresh-signature-search/SKILL.md`: ~30 lines of frontmatter with structured `metadata`. sciagent's metadata-rich style is non-standard — neither violating nor following any precedent. ADR-001's namespacing is a way to make this defensible: anything sciagent-specific is clearly tagged as such, anything generic uses the same flat keys Anthropic uses.

Sources:
- [agentskills.io specification](https://agentskills.io/specification)
- [agentskills.io specification — markdown](https://agentskills.io/specification.md)
- [agentskills/agentskills repo](https://github.com/agentskills/agentskills)
- [skills-ref reference library](https://github.com/agentskills/agentskills/tree/main/skills-ref)
- [SKILL.md spec deep dive — DeepWiki](https://deepwiki.com/agentskills/agentskills/2.2-skill.md-specification)
- [SKILL.md spec field reference — agensi.io](https://www.agensi.io/learn/skill-md-format-reference)
- [anthropics/skills repo](https://github.com/anthropics/skills)
- [skill-creator SKILL.md](https://github.com/anthropics/skills/blob/main/skills/skill-creator/SKILL.md)
