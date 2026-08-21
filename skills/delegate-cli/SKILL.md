---
name: delegate-cli
description: >-
  Delegate or consult through local codex and agy CLIs. Use when Anton asks
  either tool to implement, review, consult, or research, or when another skill
  needs model-tiered worker fan-out. Covers prompt construction, capability
  probing, safe launches, concurrency, artifacts, and review.
license: MIT
---

# Delegate through codex and agy

Use these peer CLIs for bounded work that benefits from an independent implementer, reviewer, or
researcher. Prefer codex for implementation and code review; prefer agy for web research and very
large-context synthesis. Honor a user's named tool.

Keep fan-out bounded to two or three independent units. Give every writer an exclusive file or
directory scope. Use a patch-only deliverable when exclusive paths are unavailable. Treat network
access, workspace writes, and unsandboxed execution as separate authorizations.

## Prepare the prompt

The worker does not share the current chat. Create a file-backed prompt that states:

- the task and exact deliverable;
- the absolute working directory and scoped files;
- permission boundaries, including forbidden git writes, installs, network use, and out-of-scope
  changes;
- acceptance commands and expected results;
- expected output paths for every material artifact.

Run each acceptance command before delegating it. Scope checks to the files the worker may change.
Keep runtime prompts and logs temporary by default. Record durable decisions, deviations, and outcomes
in the project's existing memory structure when they carry lasting value.

Three rules apply to every implementer:

1. Put exact numeric targets in the prompt and require the worker to report a deviation instead of
   editing an artifact to match the stated number.
2. Forbid weakening or deleting tests. Require the worker to report any blocking test by name.
3. Arrange independent review for consequential work, then re-run the checks and inspect the touched
   artifacts yourself.

Choose models by task risk and observed availability. Flash-class workers benefit from smaller units
and verification between hops. Pro- and GPT-5.x-class workers can hold broader multi-part units.

## Probe codex capabilities

Resolve asset paths relative to this `SKILL.md`. Before each codex launch, run:

```bash
"$skill_dir/assets/probe.sh" \
  --workdir "$workdir" --model "$model" --sandbox "$sandbox"
```

Add the probe's `--web` option when the task requires web access. Add `--file-read-check` before work
that depends on inspecting local files, especially in a Linux devcontainer. Stop on a nonzero result
and report the emitted facts and exit status. The probe may report that codex is absent; report that
condition without substituting another tool or installing software.

The host investigation on 2026-08-21 used `codex-cli-exec 0.147.0`. Treat that result as dated evidence
and rely on the probe for the active capability surface.

Some Linux devcontainers prevent codex file inspection because their container policy blocks the
bwrap namespace operation. See `scbio-docker/docs/ai-integration.md` for the seccomp cause and container
remediation. A bypass removes the sandbox and requires explicit user authorization; prompt rules do
not recreate sandbox enforcement.

## Launch codex

Use the launcher after a successful probe:

```bash
"$skill_dir/assets/launch.sh" \
  --unit "$unit" --workdir "$workdir" --prompt "$prompt" \
  --model "$model" --sandbox "$sandbox"
```

The asset owns codex flag spelling, file-backed stdin delivery, final-message capture, stream logging,
run directories, status records, byte counts, and locks. Pass optional launcher arguments according to
the task:

- `--rules FILE` prepends reusable rules to the task prompt.
- `--effort LEVEL` selects reasoning effort.
- `--web` requests the capability-gated web mechanism.
- `--expect PATH` names each artifact that must exist after the run.
- `--bg` returns after recording the PID and run paths.
- `--parallel-ok` affirms that a different active unit has a disjoint write scope.
- `--bypass` is reserved for an explicitly authorized unsandboxed run.

Use a stable, specific unit name. A repeated active unit indicates overlapping work and must be
resolved before relaunching. Background completion is determined from the launcher's PID and status
artifacts. Use those records instead of process-name matching.

## Run agy headlessly

The agy surface below was verified against agy 1.1.10. Its prompt flag consumes the next token, the
last prompt flag wins, and the first bare argument ends flag parsing. Supply the prompt as the value
and place every other flag afterward:

```bash
agy -p "$(<"$prompt")" \
  --model "$agy_model" \
  --print-timeout 25m \
  > "$output" 2>&1
```

Observe these parser and runtime constraints:

- A bare prompt starts interactive mode and fails headlessly on `/dev/tty`.
- `--model` accepts the display label or the slug printed by `agy models`. An invalid value exits
  before the run. Interactive sessions may rewrite the configured default.
- Omit `--effort`; agy model names already encode their effort tier and reject that option.
- `--print-timeout` requires a Go duration with a unit, such as `25m` or `300s`.
- Print mode buffers output until completion. Capture stdout with shell redirection. Asking agy to
  write the output path in the prompt switches it into agentic file-write mode.
- Headless tool calls require the appropriate permission authorization. With
  `--dangerously-skip-permissions`, confirm the log records that permissions were auto-approved.
- Read `~/.gemini/antigravity-cli/cli.log` after the run to confirm the selected model; reading it
  before the run reports the preceding session.
- A denied shell command can end a headless run without a useful answer. Grant only the command
  prefixes the scoped task needs.

For agy fan-out, use unique prompt and output paths and independent sessions. Capture each background
PID at launch and wait on that PID or on the expected output artifact.

## Review the result

Use the launcher's reported paths and byte counts to locate the codex result. A final message is a
worker report; acceptance depends on repository state and rendered artifacts.

For every delegated change:

1. Inspect every touched file and expected artifact directly from disk, including files in ignored
   directories.
2. Re-run the scoped tests and acceptance commands yourself. Report any stated count that differs
   from the observed count.
3. Check git state in every reachable repository. Treat each submodule as its own repository and
   inspect its HEAD, reflog, and worktree separately when git writes were forbidden.
4. Open rendered outputs when text length, row count, labels, or layout changed; text-only checks do
   not expose visual collisions.
5. Relay a concise result with the worker/model used, changed artifacts, verification evidence, and
   unresolved deviations.

## See also

- `reasoning-trace`
- `architecture-first-dev`
- `decision-gate-notebook`
- `interactive-breakpoint-explorer`
