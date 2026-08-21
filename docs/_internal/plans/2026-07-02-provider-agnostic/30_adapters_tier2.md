# Tier 2 — Medium: the per-harness adapter layer (shape B, built)

**Effort:** ~1–2 weeks. **Risk:** medium (touches the mutating core). **Depends on:** Tier 1.
**Delivers:** `activate` (project-level) fans a role into *any* requested harness via a dispatch table.
This is the "harness-agnostic core + thin per-harness adapters" recommendation from
`ai-research/11`, finally implemented.

---

## T2.1 — Extract the Claude backend behind a contract

Refactor the Claude-specific materialization currently inline in `activate.sh` + `claude_settings.sh`
into `lib/sciagent/harness/claude.sh`, implementing a fixed contract every adapter satisfies:

```
harness_<name>_ensure_context     # AGENTS.md shim / native; global + project
harness_<name>_ensure_skills      # symlink/point the role's skills into the harness's skill dir
harness_<name>_ensure_agents      # sub-agents (Claude/opencode/pi have them; codex/agy vary)
harness_<name>_ensure_commands    # slash commands (Claude/opencode/pi)
harness_<name>_ensure_settings    # power-user defaults in the harness's config schema
harness_<name>_ensure_hooks       # guardrail hooks (layer c) in the harness's hook surface
harness_<name>_teardown           # remove only what this adapter wrote (manifest-tracked)
harness_<name>_capabilities       # declare which of the above it supports (for honest skipping)
```

**Invariant:** extracting Claude must be behavior-preserving — the existing
`test_claude_settings_lifecycle.sh`, `test_activate_*`, `test_inject_*` suites must pass unchanged.
This is the load-bearing refactor; do it first and prove it green before adding new adapters.

`common.sh` holds the shared, already-neutral work: render CRAFT/ROLES managed block into `AGENTS.md`,
symlink `.agents/skills/`, write the manifest. Adapters call into `common.sh`.

## T2.2 — `activate --harness`

```
sciagent activate <base> [overlay] [--harness claude,pi,opencode,...|all|auto]
```

- `auto` (default) = every harness detected in the project/user env.
- Loops the requested adapters; each honors its `capabilities` and logs honestly what it skipped
  (e.g. "codex: no slash-command surface — role commands not materialized").
- Manifest gains a per-harness section so `deactivate`/`eject` tear down each cleanly.

## T2.3 — The four new adapters (capability-honest, thin)

Grounded in `research/02_cli_capability_matrix.md`. Each is ~100–200 LOC.

### codex (`~/.codex/`)
- **context:** `AGENTS.md` native (walks root→cwd). Nothing to do beyond ensuring `AGENTS.md` exists.
- **skills:** `~/.codex/skills/` (agentskills.io) — point at the role's skills, or set pi/codex to read
  `.agents/skills` where supported.
- **settings:** `config.toml` — `model`, `model_reasoning_effort`, `[mcp_servers.*]`, and per-project
  trust `[projects."<abs>"] trust_level="trusted"`. Merge via a TOML-aware step (not `jq`); tolerate
  absence of a TOML tool by copy-if-absent only. Profiles (`~/.codex/<name>.config.toml`) can encode a
  "science" profile.
- **agents/commands/hooks:** limited — codex has no first-class sub-agent spawn or slash-command surface
  comparable to Claude. `capabilities` declares these unsupported; role delegation to codex stays via
  the existing `skills/delegate-cli` (headless `codex exec -m gpt-5.5 --skip-git-repo-check`).

### agy (`~/.gemini/antigravity-cli/`)
- **context:** `AGENTS.md` **and** `GEMINI.md` native; agy even writes rules back into AGENTS.md
  (Global vs Workspace Customizations Root). Ensure AGENTS.md; skip GEMINI.md unless legacy.
- **skills:** `.agents/skills/<name>/SKILL.md` (workspace) + `~/.gemini/antigravity-cli/skills/`
  (global) — both paths are literally in the binary. **This is a direct win for the `.agents/` mirror.**
- **settings:** MCP via `mcp_config.json` (`serverUrl` for remote). **No `--model` flag** (auto-select)
  — do not try to pin a model.
- **agents:** agy has an RPC `define_subagent`/`invoke_subagent` layer — future depth, not Tier 2.
- **hooks:** `PreToolUse`/json-hooks exist in the binary; wire the no-ephemeral/caption guards if the
  hook schema is stable. `agy inspect` is the debugging oracle for "what got loaded."

### opencode (`~/.config/opencode/`, project `opencode.json`)
- **context:** `AGENTS.md` native; extra files via `instructions[]` globs.
- **skills:** `.opencode/skills/` (+ singular tolerated).
- **agents/commands:** **strongest fit after Claude** — first-class `agent` (primary/subagent) in
  `opencode.json` or `agents/*.md` (frontmatter `mode/model/tools`), invoked via the **Task tool** or
  `@mention`. Map SciAgent sub-agents → opencode subagents directly. Slash `command` supported.
- **settings:** `opencode.json` (`$schema` opencode.ai/config.json), `permission` allow/ask/deny globs,
  `mcp`. Merge JSON via `jq` (project > custom > global).
- **hooks:** `.opencode/plugins/` add tools/hooks.

### pi (`~/.pi/`, project `.pi/`)
- **context:** `AGENTS.md` **or** `CLAUDE.md` native; system prompt override via `.pi/SYSTEM.md` /
  append via `APPEND_SYSTEM.md`.
- **skills:** `.agents/skills/` native (+ reads `~/.claude/skills`, `~/.codex/skills` if listed in
  `.pi/settings.json:skills[]`). **The `.agents/` dual-track finally has a first-class consumer.**
- **settings:** `.pi/settings.json` / `~/.pi/agent/settings.json` (deep-merged), resource arrays
  `extensions[]`, `skills[]`, `prompts[]`, trust `defaultProjectTrust`.
- **agents/commands/hooks:** this is the seam to Tier 3 — the adapter registers the SciAgent pi
  extension (`"pi":{"extensions":[...]}`), which is where roles become sub-agents. Tier 2's pi adapter
  can stop at context+skills+settings; Tier 3 adds the extension.

## T2.4 — Capability matrix drives honest degradation

The toolkit must never pretend. Each adapter's `capabilities` gates what `activate` attempts and what
the status/validate output reports. Example honest line:
`"role 'science-architect' → codex: context+skills OK; sub-agents/commands unsupported (delegated via delegate-cli)."`

## T2.5 — Guardrail layer (c) across harnesses

Layer (c) — "unskippable" — is the hardest to make portable, because each harness's hook surface
differs (Claude settings.json hooks; opencode plugins; pi `pi.on("tool_call")`; agy PreToolUse; codex
execpolicy rules). Strategy:

- **Portable floor:** `sciagent validate --check ...` is already harness-agnostic (it inspects the repo,
  not the harness) and is mirrored in CI via `check-provenance`. This is the guaranteed guardrail for
  every harness — it runs regardless of which agent produced the artifacts.
- **Per-harness ceiling:** where a hook surface exists and is stable, the adapter wires the live
  no-ephemeral / caption-sweep guards. Where it doesn't (codex), rely on the validate floor.

Record the "portable floor + per-harness ceiling" split as ADR-P5.

---

## Tier-2 acceptance checks

- [ ] Claude backend extracted to `harness/claude.sh`; all existing tests green, no behavior change.
- [ ] `activate <role> --harness opencode` produces working opencode subagents + AGENTS.md + skills.
- [ ] `activate <role> --harness all` on a machine with 3 harnesses installed materializes all three;
      each manifest section tears down cleanly on `deactivate`.
- [ ] Every adapter reports capabilities honestly; unsupported surfaces are logged, not faked.
- [ ] `validate --check` floor runs identically regardless of harness.
