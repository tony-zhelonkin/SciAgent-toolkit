## Recommendation

I agree with both identified gaps.

The trace verifies that the session failure was primarily delivery: the corrected assets were absent from the pinned checkout. It also verifies that the currently shipped design remains incomplete in two respects:

- `launch.sh --bg` creates observable files but no harness-owned completion signal.
- `status.tsv` appears only after completion, so it cannot answer an in-flight progress question.

Disk observability and orchestrator wake-up are separate contracts. The revision needs both.

## 1. Status-read contract

Add a separate `assets/status.sh`. A read-only query deserves its own entry point: it requires neither a prompt nor the launch validations, and it can also observe a run started by another caller. A `launch.sh --status` mode would overload the mutating launcher.

The ordinary command should be:

```bash
"$skill_dir/assets/status.sh" --unit "$unit"
```

Also support an exact run path returned by the launcher:

```bash
"$skill_dir/assets/status.sh" --file "$status_file"
```

`--unit` reads the latest run through:

```text
${TMPDIR:-/tmp}/scio-delegate/<unit>/current/status.tsv
```

`current` should be an atomically replaced relative symlink to the run directory. Each run retains its own `status.tsv`, preserving history without duplicating status content.

### Output

Print one TSV summary followed by artifact records. Quote string values using Bash `printf %q`; use decimal integers or `-` for numeric fields.

```text
SUMMARY	unit=u05_embed_scvi	run_id=20260821T191116Z-106635	state=running	launcher_pid=106620	child_pid=106635	elapsed_s=600	heartbeat_age_s=1	activity_age_s=2	stream_bytes=425956	final_bytes=0	expected_present=0/2	codex_exit=-	result_exit=-
ARTIFACT	role=prompt	bytes=8274	mtime_epoch=1787358620	baseline_bytes=8274	baseline_mtime=1787358620	path=/tmp/.../prompt.md
ARTIFACT	role=stream	bytes=425956	mtime_epoch=1787359218	baseline_bytes=0	baseline_mtime=-	path=/tmp/.../stream.log
ARTIFACT	role=final	bytes=0	mtime_epoch=-	baseline_bytes=0	baseline_mtime=-	path=/tmp/.../final.md
ARTIFACT	role=expect	bytes=0	mtime_epoch=-	baseline_bytes=0	baseline_mtime=-	path=/work/.../result.h5ad
```

This answers “how far along” using the quantities the trace verified were useful: stream bytes, artifact bytes, mtimes, and activity age. It must not manufacture a phase, percentage, or semantic completion estimate.

### `status.sh` exit codes

For a single read:

- `0`: a complete, valid snapshot was read, regardless of run state.
- `2`: invalid arguments.
- `3`: no current status exists—unit never started in this runtime root.
- `4`: status exists but is unreadable or fails its schema/integrity checks.

The worker’s exit belongs in `codex_exit` and `result_exit`; it should not be conflated with whether the status command successfully read the record.

For harness shadowing, add:

```bash
"$skill_dir/assets/status.sh" --file "$status_file" --wait
```

`--wait` polls the same file at a fixed internal cadence, prints the terminal snapshot, and exits with `result_exit`. It adds no owner-tunable timeout or staleness knob.

## 2. What the launcher writes

Keep the existing per-run files:

```text
prompt.md
final.md
stream.log
status.tsv
pid
```

Add only the per-unit `current` symlink. Do not create separate heartbeat, byte-count, progress, or completion files.

`status.tsv` should be an atomically replaced snapshot, written through a temporary file in the same directory followed by `mv`. It should contain:

- schema version;
- unit and run ID;
- `state=starting|running|succeeded|failed|interrupted`;
- launcher and child PIDs;
- start, sample/heartbeat, last-activity, and terminal epochs;
- real Codex exit code;
- final launcher result code;
- prompt, stream, final, and expected-artifact paths;
- current bytes and mtimes;
- initial bytes and mtimes for expected artifacts.

The lifecycle should be:

1. Write `starting` before spawning Codex.
2. Write `running` immediately after recording the child PID.
3. Refresh one atomic snapshot every fixed two seconds while the child runs.
4. Advance `last_activity_epoch` when the stream, final, or an expected artifact changes in size or mtime.
5. `wait` for the actual child and capture its real status.
6. Write a terminal snapshot:
   - `succeeded`, result `0`;
   - `failed`, actual Codex status;
   - `failed`, result `30` for absent final output;
   - `failed`, result `31` for absent expected artifacts;
   - `interrupted`, with the signal-derived result when the launcher can record it.

A `SIGKILL`, machine loss, or filesystem failure may leave the last state as `running`. That limitation should remain visible through the frozen heartbeat timestamp.

The status file must not contain raw reasoning tails, prompt contents, process-name matches, appended heartbeat history, inferred percentages, or a success code derived from the existence of a wrapper banner. The trace verified that raw tails were inert and that the old wrapper’s printed exit code was decorative.

## 3. Backgrounding prescription and `--bg`

Remove `--bg` in its current form. During migration, recognizing it and failing with an actionable message is safer than silently continuing:

```text
launch.sh: --bg has been removed; run the foreground launcher through the caller's background-task mechanism
```

Prescribe callers explicitly:

- **Owner through `!`:** invoke the foreground `launch.sh` command with no `&` or `nohup`. The trace verified that timeout promotion gave these commands harness task IDs.
- **Orchestrator launching the work:** submit the foreground `launch.sh` command using the harness’s native background-task facility. The harness owns the launcher process and produces the completion notification.
- **Orchestrator observing an owner/external launch without its own task ID:** submit `status.sh --file … --wait` as a harness background task. This shadows the run and gives the harness something it can notify on.
- **Interactive shell outside an orchestrator:** foreground execution remains the default. A human who uses shell job control owns that terminal’s observation lifecycle.

The skill should state the governing rule plainly: keep `launch.sh` attached; background the attached launcher or a disk-status watcher through the harness. Disk state supports polling, while the harness task supplies wake-up.

The trace verifies the need for this distinction: every long dead-time launch was detached without a harness-owned watcher, while task-owned launches and watchers remained under two minutes.

## 4. Hung versus slow

A Bash predicate cannot distinguish “hung” from “slow” without a policy threshold or semantic progress instrumentation.

Report:

- heartbeat age;
- activity age;
- stream bytes and mtime;
- expected-artifact bytes and mtimes;
- elapsed time;
- recorded state.

Interpretation remains judgment:

- Fresh heartbeat plus old activity means the supervisor is alive and observable; Codex may be computing, blocked, or hung.
- Stale heartbeat means the status writer stopped updating; suspension, host sleep, launcher death, and filesystem trouble remain possible.
- A terminal state is definitive because the launcher reaped the child.

I therefore disagree with the trace’s requirement 4 only in wording: last activity makes hung and slow more assessable; it does not mechanically distinguish them.

## 5. Line-level `SKILL.md` changes

In the current [SKILL.md](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/skills/delegate-cli/SKILL.md), I recommend:

- Delete lines **61–62**. The dated `0.147.0` fact is superseded by runtime probing.
- Delete lines **64–65**. The probe owns detection of the file-read failure; the linked platform document owns the incident explanation.
- Retain the authorization judgment from lines **66–67**, condensed into the launch permissions guidance.
- Delete lines **79–87**. The launcher’s `--help` and implementation should own mechanics and option spelling. Line 87’s `--bg` description is specifically harmful.
- Retain and rewrite lines **88–89** as judgment:
  - `--parallel-ok` affirms both dependency independence and disjoint writes.
  - `--bypass` requires explicit authorization.
- Delete lines **91–93**. Lock refusal, PID/status interpretation, and avoidance of process-name matching should be enforced by `launch.sh` and `status.sh`.
- Replace line **129** with the concrete `status.sh` read and the returned final/stream paths.
- Delete lines **145–150** to make room for the observability prescription without growing the skill.

Rewrite lines **17–19** so “independent” means:

1. neither unit consumes an output the other may create or modify; and
2. their write scopes are disjoint.

Keep the agy section at lines **97–125** and the delegation verification checklist at **127–143**, as requested. The agy parser hazards remain mechanically unowned because no agy launcher asset exists.

## 6. Delivery

Nothing in an older pinned copy can discover this update under the stated constraints.

The failing checkout had neither the assets nor future lint logic. It has no external reference against which “stale” can be defined. A version embedded in the old skill would only demonstrate internal consistency with that same old checkout.

New lint can enforce that the current toolkit contains executable, `bash -n`-clean `probe.sh`, `launch.sh`, and `status.sh`, and can forbid the obsolete `--bg` prescription. It cannot retroactively notify a session using an older pin.

Delivery is the fleet sweep’s job: update the vendored pin and remount/reload the skill. A local version check would add ceremony without closing the verified failure mode.

## 7. Concurrency

The current lock shape is sound:

- same active unit: refuse;
- different active unit: refuse unless `--parallel-ok`.

The duplicate `u1_build` run verifies the value of the first rule. The `u08_annot_prep` violation shows that “disjoint writes” is insufficient for the second.

Redefine `--parallel-ok` as an affirmation that:

- there is no dependency edge in either direction; and
- write scopes are disjoint.

Do not add `--after`, a dependency manifest, or another scheduling DSL. Such a mechanism would require hand-maintained intent that the launcher cannot derive or verify. When a dependency exists or is uncertain, omit `--parallel-ok`; the existing live-unit refusal serializes the launches. The project runbook may document the reason for the order, while the lock mechanically prevents unaffirmed overlap.

## Disagreements and evidence boundary

I agree with the brief’s conclusions and with the trace’s primary diagnosis.

Two qualifications matter:

- The trace’s “delivery gap, decisively” is verified for the observed session. The currently shipped assets still have the two residual design gaps identified here.
- A PID binding is useful for ownership, signals, and locks. The status reader should derive its answer from the launcher-written snapshot rather than `kill -0`; PID reuse and abrupt launcher death prevent a PID probe from being durable truth.

Everything about the session timings, task notifications, false process matches, inert `tail`, duplicate run, and dependency violation is verified by the [forensic trace](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/docs/_internal/research/2026-08-21-delegation-monitoring/TRACE_jr-mc-session.md). The proposed status schema, fixed heartbeat, caller rules, and `--parallel-ok` semantics are design recommendations inferred from that evidence.