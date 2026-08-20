The recommendation is **example-first structure plus a scio-owned validator, with no new lifecycle hook**. The filesystem teaches the pattern; `scio lint` defines the measurable contract. Harness hooks may improve timing, but they cannot be the cross-harness authority.

## 1. Position on the example-versus-instruction puzzle

I agree with the hypothesis in a narrower form:

> Structure propagates when the tree contains an authentic, useful instance whose form can be copied.

A `.gitkeep`, placeholder form, or generic README is documentation about a pattern. A current `session.md`, a real reasoning note, or a figure beside its source table demonstrates the pattern.

The survey supports this strongly:

- All five STING children retained empty `handoffs/`; actual continuity records went to `sessions/`.
- Populated result stages repeatedly reproduce `{figures,tables}/`, adjacent captions, and same-stem provenance.
- The largest internal trees show that an unconstrained container teaches “put miscellaneous private material here,” producing caches, logs, environments, binaries, and a model checkpoint.

One qualification: the survey cannot establish that the tree alone caused compliance. Current project instructions explicitly route handoffs to `sessions/`, and results layout is reinforced by CRAFT, lint, and the caption hook. The evidence supports authentic examples as the strongest teacher, with concise instructions resolving ambiguity.

What follows:

- Do not mass-scaffold the ten projects lacking `docs/_internal/`.
- On the next substantive work in each project, create the first genuine record:
  - `_project/session.md` for pre-stage or cross-stage work.
  - `<stage-stem>/session.md` for stage work.
  - `<stage-stem>/reasoning/<topic>.md` when durable rationale exists.
- Create directories only with that first real file. No `.gitkeep`.
- Update `session.md` in place. Git history carries chronology.
- Treat existing memory opportunistically when the project becomes active; avoid a fleet-wide relocation exercise.

Mirroring `02_analysis/helpers/` is a category error. Helpers are reusable machinery and may serve several stages. File their rationale with the consuming stage, or under `_project/reasoning/` when the decision genuinely spans stages.

## 2. The mechanism

I recommend three neutral surfaces:

1. **A populated filesystem example.** The first real entry establishes the local pattern.
2. **One compact CRAFT router.** It says when and where to write, without restating the full grammar.
3. **A scio lint predicate.** This is the authoritative structural enforcement.

The work events are:

- **First substantive stage work:** create the stage’s first real memory file.
- **Handoff or material change of direction:** update that stage’s `session.md`.
- **Before declaring work complete or committing:** run the project’s pinned `scio lint --check internal-memory --strict`.

The proposed lint check should:

- Treat an absent `docs/_internal/` as a clean no-op.
- When present, allow `_project/` and stage stems corresponding to `02_analysis/stages/NN_<stem>.*`.
- Require each created stage directory to contain a nonempty `session.md` or `reasoning/*.md`.
- Flag `.gitkeep`, empty category directories, dated session piles, and legacy top-level containers such as `handoffs/`, `sessions/`, `reports/`, and `research/`.
- Flag caches, virtual environments, checkpoints, bytecode, logs, parquet, and other non-memory payloads.
- Warn normally and exit nonzero under `--strict`, matching current lint semantics.
- Exempt nested-repository control files if that topology is retained.

It cannot determine whether reasoning is scientifically sound, whether every important thought was recorded, or whether a note is candid and current. Those are semantic judgments.

### Surface ranking for this job

| Rank | Surface | Why | Cannot enforce |
|---:|---|---|---|
| 1 | Populated filesystem | Universal reach and almost no mechanism cost | Cannot veto, detect omissions, or judge retention |
| 2 | `scio lint` | Harness-blind, versioned with the project, deterministic | Acts only when invoked; cannot see vanished `/tmp` work |
| 3 | Git `pre-commit` calling lint | Harness-neutral and capable of blocking | Late, clone-local, bypassable, misses uncommitted work and complicates the nested repo |
| 4 | Claude/Codex lifecycle hooks | Immediate advisory or veto capability | Harness-specific configuration, trust, incomplete tool coverage, and no future-harness parity |

I recommend **no new hook now**, including no Git hook. A pre-commit hook is a reasonable escalation if evidence later shows that agents routinely skip lint: trigger `pre-commit`, run only the shared lint predicate, block on strict findings. `commit-msg` has no role because this design requires no maintained intent metadata.

The existing Claude hooks concern reproducible results and captions. They can remain orthogonal; the internal-memory design adds no third hook.

## 3. Harness parity

Codex can technically reproduce Claude’s hook pattern, but using that capability would create unnecessary project-level duplication.

I verified locally that:

- `codex-cli 0.147.0` reports `hooks` as a stable, enabled feature.
- It discovers hooks in `~/.codex/{hooks.json,config.toml}` and `<repo>/.codex/{hooks.json,config.toml}`.
- Its events include `PreToolUse`, `PostToolUse`, `Stop`, `SessionStart`, `SessionEnd`, compaction, prompt, permission, and subagent events.
- `PreToolUse` can observe and block Bash, `apply_patch`, MCP calls, and most local function tools. `Stop` can force another continuation.
- Blocking uses structured JSON or exit code 2.
- Non-managed hooks are trusted by definition hash; changed hooks are skipped until reviewed.
- Tool coverage explicitly has exceptions, so it is a guardrail rather than a complete boundary. These details agree with the current [official OpenAI Hooks documentation](https://developers.openai.com/codex/hooks).

I also checked the local installation:

- `~/.codex/config.toml` declares no hooks.
- No user hook file or installed plugin hook was present.
- The 24-project survey found no project `.codex/`.

I did not execute a live Codex hook because doing so would require adding configuration, contrary to “modify nothing.” I therefore did not empirically verify the runtime payload or trust UI. A future harness’s event surface is inherently unknown.

Codex gets parity through the neutral mechanism: it reads the same `AGENTS.md`, sees the same populated tree, and runs the same pinned `scio lint`. No project `.codex/` is needed. A future harness receives those same three surfaces.

## 4. Ownership verdict

**Named seam:** scio owns the project’s path grammar, examples, CRAFT router, and lint predicates; dev-env owns user-global habits and optional harness conveniences.

The test I applied is:

> If it must follow the repo across humans and machines, it belongs to scio; if it must follow this human across unrelated repos, it belongs to dev-env.

Therefore, “write reasoning to `docs/_internal/<stage>/…`” belongs to scio. “Leave a durable handoff before ending substantive work, using the project’s declared convention” is a dev-env-level habit.

The staged dev-env global template currently hardcodes the old `plans/`, `reasoning/`, and `sessions/` routes. Those project paths should be removed before adoption so the global layer cannot drift from a repo’s pinned scio version.

This work should not force ADR-D5. Neutral examples and lint require no harness adapter. The ADR record itself contains a demolition note saying the adapter premise was superseded, while the handoff calls D5 unimplemented; that status discrepancy can be resolved separately.

## 5. What this makes unnecessary

This design lets the toolkit remove:

- `.gitkeep` files and empty internal category scaffolds.
- Default `handoffs/`, `sessions/`, `plans/`, `reports/`, and `research/` namespaces.
- The current lint finding that treats a missing `docs/_internal/` as an error and recommends `scio link`, even though `link` does not create it.
- The current CRAFT routes to global `docs/_internal/reasoning/` and `sessions/`.
- The root `_scratch/` claim, which names a location the fleet does not use.
- Project-level `.codex/` configuration and a Codex structural-hook adapter.
- A new Claude structural hook.
- Git commit-message metadata, approval fields, and manually maintained status knobs.
- A memory mirror for `02_analysis/helpers/`.

The result hooks may remain for their existing result-specific jobs. Their references to a generic `docs/_internal/` reasoning trace should eventually use the stage route.

## 6. The first move

Delete the “`docs/_internal/` missing — run: scio link” predicate from [`lib/scio/lint.sh`](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/lib/scio/lint.sh:601), and add a test proving that a project with `docs/` but no `_internal/` passes `docs-layout`.

That single change:

- Stops pressuring ten existing projects toward empty scaffolding.
- Corrects a remediation command that cannot create the directory.
- Preserves every useful lint check once real internal content exists.
- Improves the toolkit even if the rest of this recommendation changes.

I made no modifications. The toolkit worktree remained clean; dev-env’s pre-existing modification to `scripts/aibox.sh` was left untouched.