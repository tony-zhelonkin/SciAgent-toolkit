# Prior-art map (what's already decided/proposed/deferred)

So this plan builds on the corpus instead of duplicating it. Status labels are strict.

## Philosophy (DECIDED, load-bearing)
- "The scaffold IS the interface" (`2026-06-04_mindpalace…:18`, `2026-06-23…/00_ARCHITECTURE.md:16`).
- Three-layer model: state once (CRAFT managed block) → make doable (skill/helper/command) → make
  unskippable (validate check + hook). "hooks/validate > goodwill" (`00_INDEX.md:95`).
- Ethos belongs in an always-on managed block, NOT a skill (`2026-06-04…:68`).
- Coding-style files: `craft.yaml`, `docs/guidelines/code_style.md`, `skills/{figure-style,
  reasoning-trace,scrna-pipeline-conventions,architecture-first-dev}`, `validate --check`, Claude hooks.

## Harness-agnosticism (DECIDED)
- AGENTS.md canonical; CLAUDE.md/GEMINI.md are `@AGENTS.md` shims (v3.0.0).
- Standing preferences: **harness-agnostic-preference** + **minimal-but-no-less** (`kickoff.md:25,138`).
- Provider-agnostic source-dir naming locked (`system-prompts/`, `AGENTS.md:47`).
- Recommended deployment = **shape B: harness-agnostic core + thin per-harness adapters (~150 LOC)**;
  don't depend on any single pi extension (`ai-research/11:62–66`). **Recommended, never built.**
- Two-tier skill scope `concept`/`implementation` (durability axis).
- AGENTS.md must NOT be LLM-authored (ETH finding, `00_ARCHITECTURE.md:214`).
- **v3.0.0 deleted harness installers** ("roles/agents/skills only"). Do not resurrect.

## pi extension (PROPOSED, `sciagent-extension-design-spec.md`, awaiting decision)
- Convergent backbone: `package.json` `"pi":{"extensions":[…]}`, `pi.registerTool` (discriminated-union
  actions), optional `pi.registerCommand`, optional `session_shutdown` hook (`ai-research/11:38–45`).
- **[now outdated]** claimed pi has no Task/sub-agent primitive → simulated by child `pi -p`. The
  installed package's `examples/extensions/subagent/` *is* that primitive, shipped. See
  `02_cli_capability_matrix.md`.
- Roles routing = base+overlay, depth-2, last-wins.

## Trace / reproducibility (DEFERRED as recorder; research authoritative)
- Both harnesses write JSONL: Claude `~/.claude/projects/<enc>/<uuid>.jsonl` (undocumented); pi
  `~/.pi/agent/sessions/` (documented `version:3`) + separate observability stream.
- Cross-harness standard = **sciagent-native JSONL on disk + optional OpenInference exporter**
  (`ai-research/03:67–129`). Reject OTel-GenAI/LangSmith/Langfuse as primary.
- Recorder ships iff a 2nd consumer (ADR-004 optimizer or ADR-006 bench) ships (`ai-research/11:96`).
- The "missing ADR": promote **ADR-011 — sciagent trace/result substrate** (`ai-research/10:220`).
- What shipped instead (v3.2.0): CRAFT no-ephemeral + `reasoning-trace` skill + PreToolUse hook +
  `validate --check provenance` (markdown traces, deliberately lightweight).

## Resolved/deferred ADRs (`kickoff.md §9`, `deffered.md`)
001 metadata namespace (DECIDED) · 002 `requires` graph (DECIDED) · 003 two-tier scope (DECIDED) ·
004 ERA optimizer (DEFERRED) · 005 trace recorder+distiller (DEFERRED) · 006 benchmark (DEFERRED) ·
007 doctor→thin `validate` (DISSOLVED) · 008 lab-loop (DROPPED) · 009 uv pinning (DEFERRED) ·
010 cloud-lab MCP (DEFERRED). Latent: stack-depth-2 "rubber-banding"; empirical-claims audit.

## The gap this plan fills
Prior corpus = in-repo context injection + one deferred trace/optimizer thread. **Not covered:**
(1) installer/provisioning orchestration for codex/agy/pi/opencode (only `delegate-cli` *documents*
codex/agy invocation; opencode unmentioned; shape-B adapters never built); (2) user-level settings
provisioning (only just-landed project-level, Claude-only); (3) session-persistence fix. This plan's
novel surface = the provisioning + cross-provider user-settings + session-persistence layer *above* the
existing injection core, realizing shape B for the actually-used harnesses.
