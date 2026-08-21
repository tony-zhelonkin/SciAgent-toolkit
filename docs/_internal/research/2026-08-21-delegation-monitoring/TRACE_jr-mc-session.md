# TRACE — delegation monitoring in the JR-MC session

Forensic reconstruction of how codex work units were launched and observed in the
JR-MC session on 2026-08-21. Read-only investigation; nothing under
`/workspaces` or `/data2` was written and no git command that writes was run.

## Provenance of the evidence

| Source | How it was read | Status |
|---|---|---|
| Transcript `f03e8c64-1aa8-459e-98a5-46d96fd52a7d.jsonl`, 5,229,845 B / 1,621 lines | `docker cp` to `/tmp/jrmc-trace/t.jsonl`, then `jq` field extraction | verified by reading |
| `/workspaces/JR-MC/.devcontainer/scripts/codex_unit.sh`, 1,845 B | `cat` in container | verified by reading |
| `/workspaces/JR-MC/logs/codex/` mtimes and sizes | `ls -la --time-style=full-iso` | verified by reading |
| `/workspaces/JR-MC/01_modules/SciAgent-toolkit` at `e33b635` | `git log -1`, `ls`, `sed` | verified by reading |
| Background task output files under `/tmp/claude-788715489/.../tasks/` | `ls -l`, `cat` | verified by reading |
| Live process table | `ps -ef`, `ps -o etime=` | verified by running |

**Times.** Transcript timestamps are UTC. Container filesystem mtimes are CDT
(UTC−5). Every time in this document is UTC; log mtimes have been converted by
adding five hours. Completion time for a unit is taken as the mtime of
`logs/codex/<unit>.final.md`, which is the moment codex wrote its final message.

**Two runs were in flight during this investigation** and were left alone:
`u08_annot_prep` (pid 104212, launched 19:16:56Z, exited 19:32:30Z during the
investigation) and `u05b_pilot_fix` (pid 106635, launched 19:23:51Z, still
running at 19:33:53Z). No signal was sent to either.

---

## 1. Timeline of delegation attempts

Twenty-one launches. Six were typed by the owner behind `!`; fifteen were issued
by the orchestrator. "First true observation" is the first transcript event after
completion in which the orchestrator learned the run's real state — a check that
returned the truth, not a probe that returned a false positive.

| # | Time (Z) | Unit | Command | Detach mechanism | Completed (Z) | First true observation (Z) | Dead time | What woke the orchestrator |
|---|---|---|---|---|---|---|---|---|
| 1 | 05:03:20 | u5_references | `codex_unit.sh u5_references --search &` | owner `!` + `&` | died on bad flag | 05:03:48 | ~0 | own `sleep 20` probe |
| 2 | 05:05:26 | u1_build | `codex_unit.sh u1_build` | owner `!`, **foreground**; harness promoted at 120 s → task `bs3gpurlh` | 05:15:04 | 05:15:41 | 0.6 min | own artifact-wait watcher `b47cnd64g` |
| 3 | 05:08:34 | u5_references | `codex_unit.sh u5_references --web &` | owner `!` + `&` | 05:38:00 | 05:39:50 | 1.8 min | own turn |
| 4 | 05:10:01 | u1_build (duplicate) | `codex_unit.sh u1_build` | owner `!`, foreground; promoted → `b2lwx21kq` | killed 05:12:49 | 05:12:29 | n/a | owner `status?` at 05:12:02 |
| 5 | 05:17:33 | u1b_symbol_fix | `codex_unit.sh u1b_symbol_fix` | owner `!`, foreground; promoted → `bj4zvpzmh` | 05:19:58 | 05:21:40 | 1.7 min | task notification 05:21:31 |
| 6 | 05:24:02 | u1c_stage_assertion | `codex_unit.sh u1c_stage_assertion` | owner `!`, foreground; promoted → `bngzbq4qf` | 05:26:42 | 05:26:50 | 0.1 min | owner message 05:26:02 |
| 7 | 05:28:07 | u2_droplet | `nohup codex_unit.sh u2_droplet > /dev/null 2>&1 &` | orchestrator `nohup`, no watcher | 05:44:36 | 11:34:26 | **349.8 min** | owner message 11:33:22 |
| 8 | 11:35:03 | u3_qc | `nohup … u3_qc &` + `bg` watcher on `/proc/34443` | `nohup` + watcher `bjb9zf6uj` | 12:02:17 | 12:02:33 | 0.3 min | task notification 12:02:23 |
| 9 | 12:04:32 | u3b_scrublet_threshold | `nohup … u3b_scrublet_threshold &` | `nohup`, **no watcher** | 12:20:23 | 16:00:03 | **219.7 min** | owner message 15:59:13 |
| 10 | 16:02:33 | u3c_qc_removal | `nohup … u3c_qc_removal &` | `nohup`, no watcher | 16:34:43 | 16:40:43 | 6.0 min | owner message 16:39:32 |
| 11 | 16:44:08 | u3d_qc_figures | `nohup … u3d_qc_figures &` | `nohup`, no watcher | killed 16:46:52 | 16:46:52 | n/a | owner message 16:45:53 |
| 12 | 16:47:55 | consult_gpt55 | `nohup codex exec -m gpt-5.5 … &` | `nohup` + pid watcher `bx3ohteug` | 16:53:54 | 16:54:11 | 0.3 min | task notification 16:54:06 |
| 13 | 16:49:50 | consult_doublet | `nohup codex exec … &` | `nohup` + group watcher `bloee0e6l` | 16:54:20 | 16:55:47 | 1.5 min | own turn |
| 14 | 16:50:30 | consult_ambient | `nohup codex exec … &` | `nohup` (same watcher) | 16:55:28 | 16:55:47 | 0.3 min | own turn |
| 15 | 17:04:56 | u3d_qc_figures (relaunch) | `nohup … u3d_qc_figures &` | `nohup`, no watcher | 17:46:27 | 18:03:36 | **17.2 min** | owner message 18:03:13 |
| 16 | 18:05:05 | consult_plan | `nohup codex exec … &` | `nohup`, no watcher | 18:08:10 | 18:32:26 | **24.3 min** | owner `status?` 18:32:20 |
| 17 | 18:34:00 | u8_restructure | `nohup … u8_restructure &` | `nohup`, no watcher | 19:01:41 | 19:10:12 | **8.5 min** | owner nag 19:09:18 |
| 18 | 19:11:16 | u05_embed_scvi | `codex_unit.sh u05_embed_scvi 2>&1 \| tail -5` | Bash `run_in_background: true` → `bapxzb86w` | 19:22:04 | 19:22:08 | 0.07 min | task notification 19:22:04 |
| 19 | 19:11:23 | u08_annot_prep | `codex_unit.sh u08_annot_prep 2>&1 \| tail -5` | `run_in_background: true` → `buxipoqgp` | 19:15:54 | 19:15:58 | 0.07 min | task notification 19:15:54 |
| 20 | 19:16:56 | u08_annot_prep (relaunch) | `… u08_annot_prep 2>&1 \| tail -5` | `run_in_background: true` → `bo7tjmo6s` | 19:32:30 | after the session window | — | — |
| 21 | 19:23:51 | u05b_pilot_fix | `… u05b_pilot_fix 2>&1 \| tail -4` | `run_in_background: true` → `b265q0woe` | still running | — | — | — |

Three detach mechanisms appear, and they are not equivalent:

- **Owner `!` in the foreground** (#2, #4, #5, #6). The harness promoted these at
  the 120-second mark. The tool result is explicit about the mechanism:
  `Command did not complete within its 120s timeout and was moved to the
  background (ID: bs3gpurlh).` Promotion produces a task id, which is what
  eventually generated the notification at 05:21:31.
- **`nohup … &` inside a synchronous Bash call** (#7 and eleven more). The Bash
  call returns in seconds; the codex process is reparented and the harness has no
  record of it. There is no task id and therefore no completion notification.
  Eleven of the fifteen orchestrator launches used this.
- **`run_in_background: true`** (#18–#21). Verified from the tool input:
  `jq` on the `tool_use` block shows `run_in_background` set to `true`, and the
  harness answered `Command running in background with ID: bapxzb86w. … To check
  interim output, use Read on that file path.`

### The shape of the orchestrator's day

Every idle gap in the session is bounded on the right by a user message. The
orchestrator never resumed on its own except through a background-task
notification.

```
GAP 152.7 min : 02:29:57 (own turn) -> 05:02:36 (owner)
GAP 352.1 min : 05:41:13 (own turn) -> 11:33:22 (owner)
GAP  17.8 min : 11:44:33 (own turn) -> 12:02:23 (task notification, u3 watcher)
GAP 234.3 min : 12:04:57 (own turn) -> 15:59:13 (owner)
GAP  36.6 min : 16:02:57 (own turn) -> 16:39:32 (owner)
GAP  21.0 min : 17:24:09 (own turn) -> 17:45:08 (owner "status?")
GAP  17.6 min : 17:45:38 (own turn) -> 18:03:13 (owner)
GAP  26.7 min : 18:05:37 (own turn) -> 18:32:20 (owner "status?")
GAP  34.9 min : 18:34:26 (own turn) -> 19:09:18 (owner nag)
```

Only one of the nine gaps ends in a notification, and that is the gap covering
`u3_qc`, the single unit for which the orchestrator armed a `/proc` watcher.

---

## 2. The dead-time measurements

### Worst case by wall clock: u2_droplet, 349.8 minutes

Launched 05:28:07Z with `nohup … &` and no watcher. `logs/codex/u2_droplet.final.md`
carries mtime 05:44:36Z. The orchestrator's last turn of the night ended at
05:41:13Z, three and a half minutes before the run finished. Nothing observed it
until 11:34:26Z, when — one minute after the owner's `ok, makes sense proceed` at
11:33:22Z — the orchestrator ran `ps -p 24787 … && echo "u2 RUNNING …" || echo "u2 exited"`
and got `u2 exited`. Dead time from completion to observation: 5 h 49 min 50 s.

This window overlaps the owner's likely absence overnight, so it overstates
inattention as a *behaviour*. It is nonetheless the true measure of how long a
finished run sat unread.

### Worst case with the owner demonstrably present: u3b_scrublet_threshold, 219.7 minutes

Launched 12:04:32Z at the tail of a Bash call that also wrote the unit prompt:

```
… UNIT
nohup .devcontainer/scripts/codex_unit.sh u3b_scrublet_threshold > /dev/null 2>&1 &
sleep 5; pgrep -f 'codex exec -m gpt-5.6-sol' | head -1
```

The orchestrator's turn ended at 12:04:57Z. `u3b_scrublet_threshold.final.md`
carries mtime 12:20:23Z — sixteen minutes later. The next event in the transcript
is the owner at 15:59:13Z. Between those two points the orchestrator issued zero
tool calls.

Its first statement after the owner's message asserted the opposite of the truth:

> `2026-08-21T15:59:36.907Z` — "Three things: `docs/_internal` as its own repo,
> then config knobs, then the removal unit. Setting up the repo first (**u3b is
> mid-run**, so I'll hold config edits until it exits):"

It then probed, and the probe agreed with it — falsely:

> `2026-08-21T15:59:55.557Z` — `pgrep -f 'codex exec -m gpt-5.6-sol' >/dev/null && echo "u3b RUNNING" || echo "u3b done"; ls -la logs/codex/u3b*.final.md`
> → `u3b RUNNING` … `-rw-r--r-- 1 devuser devgroup 4935 Aug 21 07:20 logs/codex/u3b_scrublet_threshold.final.md`

The `pgrep` pattern appears inside the polling shell's own command line, so it
matched itself. The `ls` on the same line already carried the answer — a
completed final report, timestamped 07:20 local. Eight seconds later a
`ps -eo pid,ppid,etime,args | grep -E 'codex exec|codex_unit' | grep -v grep`
returned nothing, and the orchestrator corrected itself:

> `2026-08-21T16:00:13.890Z` — "u3b finished (**the earlier check raced**)."

It did not race. It self-matched, which is a determinate failure of that probe
and is documented verbatim in the delegate-cli SKILL.md the session had loaded
(see §5). Dead time from completion to true observation: 3 h 39 min 40 s.

### The nag, measured backwards

The owner's message anchoring the complaint:

> `2026-08-21T19:09:18.592Z` — "I think you need to be fanning out codex as
> background jobs that you can actually monitor for, otherwise we are lost ,
> please proceed moving onward"

Counting back from that timestamp:

| Interval | Duration | What was happening |
|---|---|---|
| 18:34:26 → 19:09:18 | 34 min 52 s | orchestrator silent; `u8_restructure` in flight until 19:01:41, then finished and unread for 7 min 37 s |
| 18:32:54 → 18:34:26 | 1 min 32 s | orchestrator working; launched u8 |
| 18:05:37 → 18:32:20 | 26 min 43 s | orchestrator silent; `consult_plan` in flight until 18:08:10, then finished and unread for 24 min 10 s |

In the 95 minutes before the nag, the orchestrator was silent for 61.6 of them,
two detached runs completed, and neither completion was observed until the owner
spoke. The first thing it did after the nag was concede the mechanism:

> `2026-08-21T19:09:47.192Z` — "You're right — `nohup` detaches and I lose the
> completion signal. Switching to monitored background jobs."

and then, at 19:10:12.968Z, discover that the run it had launched 36 minutes
earlier was already over: "u8 already exited. Checking whether it completed or
died:".

### All measured dead times, sorted

| Unit | Dead time (completion → true observation) |
|---|---|
| u2_droplet | 349.8 min |
| u3b_scrublet_threshold | 219.7 min |
| consult_plan | 24.3 min |
| u3d_qc_figures (relaunch) | 17.2 min |
| u8_restructure | 8.5 min |
| u3c_qc_removal | 6.0 min |
| u5_references | 1.8 min |
| u1b_symbol_fix | 1.7 min |
| consult_doublet | 1.5 min |
| u1_build | 0.6 min |
| u3_qc | 0.3 min |
| consult_gpt55 | 0.3 min |
| consult_ambient | 0.3 min |
| u1c_stage_assertion | 0.1 min |
| u05_embed_scvi | 0.07 min |
| u08_annot_prep | 0.07 min |

The split is clean and it is entirely mechanical. Every dead time above six
minutes belongs to a `nohup … &` launch with no watcher. Every dead time under
two minutes belongs to a run the harness knew about — either promoted from the
owner's foreground `!`, launched with `run_in_background: true`, or shadowed by a
watcher loop that the harness *did* know about.

**Time to act on the result.** For every unit, the orchestrator acted within one
tool call of its first true observation: it read `final.md`, ran the guardrail
grep for git writes, and re-derived the numbers itself. Observation-to-action
latency is negligible throughout — 4 s for u05, 4 s for u08, 30 s for u3b. The
loss is entirely in launch-to-observation, never in observation-to-action.

---

## 3. What the orchestrator did while a run was in flight

Of 222 Bash calls in the session, 51 (23 %) contained a process-liveness probe
(`pgrep`, `ps aux`, `ps -eo`, `ps -o pid`, `ps -p`, `kill -0`, `/proc/`) and 55
touched `logs/codex`. Twenty-two calls combined a process probe with a log or
`final.md` check — the recurring "is it alive and has it written anything" pair.

**Never polled the harness.** Tool-call counts for the whole session: `Bash` 222,
`Read` 5, `Skill` 1, `AskUserQuestion` 1. There is no `BashOutput` call and no
`TaskOutput` call, at all. The five `Read` calls were all PNG figures; not one
read a `tasks/*.output` file, even though the harness told it to in the tool
result for every background launch: "To check interim output, use Read on that
file path."

**Behaviour split by in-flight window:**

- *Productive concurrency.* While `u2_droplet` ran (05:28–05:44) the orchestrator
  installed `scikit-image`, appended a lane-comparison addendum to the u3 prompt,
  and drafted the u4/u4b/u4c/u6/u6b unit prompts — a dozen tool calls of real
  work. Same pattern while `u3_qc` ran (11:35–12:02): it installed and configured
  `loupeR`, wrote the u7 prompt, and armed the `/proc` watcher.
- *Zero-turn windows.* `u3b` (12:04–12:20), `u3c` (16:02–16:34), `consult_plan`
  (18:05–18:08) and `u8_restructure` (18:34–19:01) each ran to completion with
  **zero** tool calls from the orchestrator. In each case its turn had ended
  seconds after the launch and nothing brought it back.
- *Owner-driven interleaving.* While `u3d` ran the second time (17:04–17:46) the
  orchestrator did substantial unrelated work — reading the `consensus` package
  source, correcting the u6b branch design, writing a 252-line `STATE.md`. Every
  one of those turns began with an owner message (17:11:14, 17:22:22, 17:45:08).
- *Repeated identical checks.* The `ps -eo comm|grep -c '^codex$'` idiom appears
  at 16:54:11, 16:55:47, 17:22:33, 17:45:17, 18:03:36, 18:32:26 and 19:16:19 —
  seven identical liveness counts, each the opening move of a turn the owner had
  triggered.

**Polling was structurally unavailable in-turn.** At 05:08:39 the orchestrator
tried `sleep 60; … tail -15 …` and a guardrail hook refused it:

> `<tool_use_error>Blocked: sleep 60 followed by: … To wait for a condition, use
> Monitor with an until-loop (e.g. `until <check>; do sleep 2; done`). To wait
> for a command you started, use run_in_background: true. Do not chain shorter
> sleeps to work around this block.</tool_use_error>`

It complied and used the recommended pattern five times — `b27c8jk84`
(`until [ -s logs/codex/u1_build.final.md ]`), `b47cnd64g`, `bjb9zf6uj`
(`while [ -d /proc/34443 ]`), `bx3ohteug` (`while ps -p 69248`), `bloee0e6l`.
All five worked. It then stopped using them after 16:50:53 and launched six more
`nohup` runs bare. That abandonment is the proximate cause of the long dead
times, and nothing in the transcript explains it.

---

## 4. Whether `| tail -5` produced observable output before completion

It produced nothing, and this is measured rather than reasoned.

The two most recent launches were still in flight during this investigation.
Their harness output files:

```
2026-08-21T19:33:53Z
-rw-r--r-- 1 devuser devgroup   0 Aug 21 14:23  b265q0woe.output   # u05b_pilot_fix, 10 min elapsed
-rw-r--r-- 1 devuser devgroup 199 Aug 21 14:32  bo7tjmo6s.output   # u08_annot_prep, exited 19:32:30Z
```

`u05b_pilot_fix` had been running for exactly ten minutes (`ps -o etime=` reported
`10:00`) with **0 bytes** in the file the orchestrator would have to read, while
its own `logs/codex/u05b_pilot_fix.run.log` had grown to 425,956 bytes. The
information existed on disk; the pipe the orchestrator chose could not carry it.

`u08_annot_prep`'s relaunch was observed at 0 bytes and became 199 bytes at the
instant it exited — 15 min 34 s of a run with zero observable output, then the
whole of it at once. Its stream log at that moment stood at 1,368,614 bytes.

The 199 bytes are the wrapper's own banner and nothing else:

```
==> unit u05_embed_scvi  model gpt-5.6-sol  effort high
    final : logs/codex/u05_embed_scvi.final.md
    stream: logs/codex/u05_embed_scvi.run.log
==> u05_embed_scvi exited 0

[exited with code 0]
```

Three lines the orchestrator wrote itself, plus an exit line that can only ever
say `0` (see §5), plus the harness's own footer. `tail -5` selected exactly the
four lines that carry no information, and buffered them until they were moot.

**What made the harness say "Running in the background."** For #18–#21 it was an
explicit parameter, not a timeout promotion: `jq` on the `tool_use` input shows
`run_in_background` set to `true`, and the harness replied `Command running in
background with ID: …`. For the owner's `!` runs the same phrase came from
promotion, with different wording that names the cause: `Command did not complete
within its 120s timeout and was moved to the background`.

**Could the output be read afterwards?** Yes. The path was in the tool result and
again in the completion notification. The orchestrator never read it, going
straight to `logs/codex/<unit>.final.md` instead — which was the right call, since
the `.output` file held only the banner.

**What actually worked here.** Despite `| tail -5` being inert, #18 and #19 have
the two lowest dead times in the whole session (4 seconds each). The value of
`run_in_background: true` was never the output; it was the task id, and the
completion notification that id makes possible. The orchestrator's own summary
credited the wrong half:

> `2026-08-21T19:11:37.196Z` — "Two units running as **monitored** background
> jobs — I'll be notified on completion rather than having to poll"

The first clause is inaccurate — nothing was monitored — and the second is exactly
right.

---

## 5. Overlap

### Confirmed: two `u1_build` runs, 2 min 48 s

The owner typed `.devcontainer/scripts/codex_unit.sh u1_build` at 05:05:26.674Z
and again at 05:10:01.946Z. The orchestrator noticed at 05:12:27 and read the
process table:

> `2026-08-21T05:12:43.935Z` — "Two `u1_build` runs are racing — pids 14745
> (7 min) and 16233 (2.5 min), both writing the same files and both appending to
> the same log. Killing the younger duplicate:"

The kill went out at 05:12:44.439Z (`kill -TERM -- "-16228"`, the process group)
and by 05:12:49.445Z only pid 14745 remained. **Overlap: 05:10:01.9 → 05:12:49.4,
i.e. 2 min 47.5 s** — the earlier "~2.5 minutes" estimate is confirmed, and 2m48s
is the tight upper bound.

The collision is structural, independent of what codex chose to write. Both
invocations of the wrapper redirect to the same two paths:

```
    -o "${LOGS}/${UNIT}.final.md" … - > "${LOGS}/${UNIT}.run.log" 2>&1
```

`>` truncates, so the second launch destroyed the first's stream log at
05:10:01. Both also ran stage `00_build` writing the same 258 MB
`00_data/processed/00_raw.h5ad`. One detail of the orchestrator's own account is
wrong: they were not "appending" to the same log, they were truncating it.

### Confirmed: a dependency violation launched as parallelism

`u05_embed_scvi` (19:11:16) and `u08_annot_prep` (19:11:23) were launched seven
seconds apart with an explicit independence claim:

> `2026-08-21T19:11:37.196Z` — "These are genuinely independent: u05 needs
> `03_qc.h5ad`, u08 needs `03_qc.h5ad` plus stage 02's `votes_azimuth.csv`.
> Neither writes a path the other reads."

u08 came back blocked 4.5 minutes later, and the first of its four findings
refutes the claim:

> `2026-08-21T19:17:11.910Z` — "1. **It required `05_embed.h5ad`** — a dependency
> my renumbering carried over from the old u6 … Fixed to depend on `03_qc.h5ad`,
> which makes it genuinely parallel to u05 rather than falsely serial."

The prompt the orchestrator itself had written made u08 wait on u05's output. The
overlap wasted a full unit run.

### Concurrency inventory

| Window (Z) | Concurrent runs |
|---|---|
| 05:10:01–05:12:49 | u1_build ×2 — **same unit, same paths** |
| 05:17:33–05:19:58 | u5_references + u1b_symbol_fix |
| 05:24:02–05:26:42 | u5_references + u1c_stage_assertion |
| 05:28:07–05:38:00 | u5_references + u2_droplet |
| 16:49:50–16:53:54 | consult_gpt55 + consult_doublet + consult_ambient (all read-only by prompt) |
| 19:11:23–19:15:54 | u05_embed_scvi + u08_annot_prep — **false independence** |
| 19:23:51– | u08_annot_prep (relaunch) + u05b_pilot_fix — in flight, not evaluated |

**Could not determine: whether the u5/u1b, u5/u1c and u5/u2 pairs wrote the same
files.** Intersecting the repo paths *mentioned* in each run log gives overlaps
(`02_analysis/config/analysis_config.yaml` appears in all four), but the run logs
record no distinction between a read and a write and carry no timestamps, so a
mention cannot be promoted to a concurrent write. The orchestrator asserted
disjointness at 05:17:17 — "u5 is still going in the background; these don't
overlap" — and nothing in the evidence contradicts or confirms it.

### A separate collision the wrapper causes on its own

`u08_annot_prep.final.md` carries mtime 19:15:54Z (the first run) while
`u08_annot_prep.run.log` carries 19:32:30Z (the relaunch). The wrapper's paths
are keyed on unit name alone, so the relaunch at 19:16:56 overwrote the blocked
run's stream log. The evidence for *why* u08 blocked the first time no longer
exists on disk.

---

## 6. `codex_unit.sh` capability audit

1,845 bytes, tracked in git (`git ls-files --error-unmatch` succeeds). It reads
`docs/_internal/codex/00_house_rules.md` plus `docs/_internal/codex/<unit>.md`,
pipes the concatenation to `codex exec` with `-m gpt-5.6-sol`,
`model_reasoning_effort=high`, `--dangerously-bypass-approvals-and-sandbox`,
`-C $ROOT`, `--skip-git-repo-check`, and redirects.

**What it gives the orchestrator for observing a run:**

| Question | Provided? |
|---|---|
| pid file | **No.** Nothing is written anywhere. |
| status file | **No.** |
| exit code on disk | **No.** |
| lock / duplicate refusal | **No.** Two runs of one unit start without complaint. |
| per-run directory | **No.** Paths are `logs/codex/<unit>.{final.md,run.log}`, so a relaunch destroys the prior run. |
| byte counts | Only implicitly, as the size of `run.log` — which the orchestrator had to discover by `stat -c %s`. |
| progress markers or timestamps in the stream | **No.** |
| a completion signal | Only the appearance of `final.md`, and only if codex reaches the end. |

Everything it prints goes to stdout at two moments: three banner lines before the
run, and one line after.

**The exit line cannot report a failure.** The script opens with
`set -euo pipefail` and closes with:

```bash
    - > "${LOGS}/${UNIT}.run.log" 2>&1

echo "==> ${UNIT} exited $?"
```

Under `set -e`, a non-zero pipeline aborts the script before the `echo` runs, so
the only value `$?` can hold at that line is `0`. Both completed `.output` files
in the tasks directory say `exited 0`, and they would say so regardless. The exit
status is decorative.

**`logs/` is gitignored** (`.gitignore:17:logs/`), so neither the final report nor
the stream is tracked. The filesystem is the only record, and it does not survive
a rebuild.

**The project's own runbook does not fill the gap.** `docs/_internal/codex/RUN.md`
documents the invocation, the model, the dependency order, and a three-command
post-run verification (`git status --porcelain`, submodule reflog, guardrail
grep). It says nothing about observing a run in progress. Its only nod to
concurrency is a prose note that units 5 and 1 "are independent and may run
concurrently" — an assertion about two specific units, not a mechanism.

---

## 7. Delivery gap or design gap

**Delivery gap. Decisively.**

The project's toolkit checkout at `/workspaces/JR-MC/01_modules/SciAgent-toolkit`
is pinned at `e33b635` (Thu Aug 20 18:05:14 2026 −0500). Its `skills/delegate-cli/`
holds one file:

```
-rw-r--r-- 1 devuser devgroup 27620 Aug 20 20:37 SKILL.md
ls: cannot access '.../skills/delegate-cli/assets/': No such file or directory
```

416 lines, no `assets/`. The orchestrator loaded it at 02:08:15.703Z, immediately
after the owner's "Have codex do code implementations via delegate-cli with the
main model being gpt-5.6-sol". That is the guidance the session followed.

**What the old text told it about backgrounding.** Three passages, and all three
are Claude-Code-tool advice or verified-trap warnings, with no mechanism behind
either:

> line 168 — "2. **Run in background** if it's real implementation
> (`run_in_background: true`) — these take minutes."

> line 286 — "- **Background the slow ones** (`run_in_background`), each with its
> own scratch OUT path; read the OUT file on the completion notification, then
> relay — don't dump the raw stream."

> line 211 — "**Pattern for Anton to run:** `! codex exec --skip-git-repo-check
> -m gpt-5.5 -s danger-full-access - < /abs/prompt.md` — Anton runs in the `!`
> prefix; I prepare the prompt file and monitor output."

The advice at 168 and 286 is a Claude Code tool parameter. It is unavailable to
the owner typing behind `!`, which is precisely the mode line 211 and line 200
prescribe for this container — and line 211's "I … monitor output" names an
obligation with no instrument attached. The skill offers no way to observe a run
the orchestrator did not itself launch through the harness.

The old text *does* carry the two traps the session then walked into, at length:

> line 298 — "## Waiting on background workers — the self-match trap (verified
> 2026-07-27) … **The loop's own command line contains the pattern**, so
> `pgrep -f` matches the waiter's shell and the condition stays true after every
> worker has exited."

> line 342 — "- **Never `pkill -f <script>` to clean up** — same self-match
> problem, and here it killed the invoking shell (exit 144)."

Both were reproduced verbatim. Three false `RUNNING` verdicts came from
`pgrep -f 'codex exec -m gpt-5.6-sol'` self-matching (15:59:55 on u3b, 16:40:20 on
u3c, 16:46:27/16:46:39 after the u3d kill), and at 16:46:32 the orchestrator ran
`pkill -TERM -f 'codex exec -m gpt-5.6-sol'` and got back **`Exit code 144`** —
the same exit code the skill names. Documentation of a trap did not prevent the
trap; only at 16:46:52 did the orchestrator name it, and by 15:59 it had already
misread u3b's state for the previous four hours.

**The mechanism existed and did not arrive.** Commit `2a1c702`, "Let the assets
spell the flags the skill kept getting wrong" (Fri Aug 21 00:44:13 2026 −0500 =
05:44Z), added `skills/delegate-cli/assets/probe.sh` and
`skills/delegate-cli/assets/launch.sh` and cut the skill from 416 lines to 150.
`launch.sh` writes exactly what was missing: a per-run directory
`$run_root/$unit/$run_id`, a `pid` file, a `status.tsv` recording `codex_status`,
`status` and per-path `bytes` for the prompt, final, stream and every `--expect`
artifact, plus a `mkdir`-based lock at `$run_root/locks/$unit.lock` holding the
pid, which refuses a duplicate unit and refuses a *different* concurrent unit
unless `--parallel-ok` is passed. The rewritten SKILL.md is explicit:

> "Background completion is determined from the launcher's PID and status
> artifacts. Use those records instead of process-name matching."

> "A repeated active unit indicates overlapping work and must be resolved before
> relaunching."

That commit landed 05:44Z. The u3b failure was at 12:04Z, u3c at 16:02Z, u3d at
16:44Z, consult_plan at 18:05Z, u8 at 18:34Z, and the nag at 19:09Z. The fix
existed in the toolkit repository for **13 h 25 min** before the owner complained,
and never reached the project, because nothing re-pins a vendored
`01_modules/SciAgent-toolkit` mid-session. The design was right and undelivered.

One design residue is worth separating out. The old skill's "run in background"
advice was *ambiguous about who runs the command*, and the session's ambiguity
mattered: the owner's `!` runs got a task id by accident of timeout promotion and
therefore behaved well, while the orchestrator's `nohup` runs got nothing. The new
text's `--bg` plus status records is agnostic about the caller, which is the
correct shape. But the delivery channel is the failure here.

---

## 8. What a correct mechanism would have had to expose

Derived from the questions the orchestrator actually put to the filesystem and
could not get answered. Each row is a question it asked, how it asked, and what
came back.

| Question it asked | How it asked | What it got |
|---|---|---|
| Is unit X still running? | `pgrep -f 'codex exec -m gpt-5.6-sol'`, `ps -eo comm\|grep -c '^codex$'`, `ps -p <pid>`, `[ -d /proc/<pid> ]` — 51 calls | A boolean over a *pattern*, wrong three times by self-match, and never unit-scoped. No pid file exists to ask instead. |
| Which unit does this pid belong to? | at 05:12:27, `ps -o pid=,ppid=,etime=,cmd= \| grep 'codex e[x]ec' \| sed 's/\(-o [^ ]*\).*/\1/'` — reading the unit name back out of the `-o` argument | Worked, by scraping a command line. There is no registry mapping unit → pid. |
| When did this run last do something? | at 05:27:39, `grep -oE '^\[[0-9-]+T[0-9:]+\]' logs/codex/u5_references.run.log \| tail -2` | **Empty.** The stream carries no timestamps. The section header printed with nothing under it. |
| How far along is it? | `stat -c %s` on the run log (8 calls), `tail -c 400`/`tail -c 3000` on the stream (7 calls), and once the size of a `.building` temp file | A monotonically increasing byte count and a slice of raw reasoning text. No phase, no step, no percentage. |
| Did it exit, and how? | `ls -la logs/codex/<unit>.final.md` (19 calls) | Existence of a file as a proxy for success. No exit code on disk; the wrapper's `exited $?` line can only print `0`. |
| Is another run already holding this unit? | never asked — nothing to ask | The u1_build duplicate ran 2 min 48 s. |
| What is the in-flight run writing right now? | only after the fact, `git status --porcelain` | Nothing during the run. |
| Where did the *previous* run of this unit go? | implicit in reading `u08_annot_prep.run.log` after the relaunch | Overwritten. One path per unit, no run id. |

A mechanism that answers "is u05 still running, and how far along" in one cheap
command would therefore have to expose, from a single read and without touching
the process table:

1. **unit → pid binding on disk**, so liveness is `kill -0` on a recorded pid
   rather than a pattern match that can match the asker.
2. **a terminal state and an exit code**, written by the launcher after the child
   reaps, distinguishable from "never started" and from "still running".
3. **a monotonic progress quantity** — byte counts of the stream, the final
   report, and each declared expected artifact — so "how far along" has an answer
   that does not require reading megabytes of reasoning text.
4. **a last-activity timestamp**, so a hung run is distinguishable from a slow
   one. The session had no way to tell these apart and never tried.
5. **a per-run identity** (`unit/run_id`) so a relaunch does not destroy the
   record of the run it replaces.
6. **a lock keyed on the unit**, so a duplicate launch fails loudly instead of
   racing, and a launch beside a *different* live unit requires an explicit
   affirmation that the write scopes are disjoint.
7. **declared expected artifacts**, so "did it actually produce anything" is a
   question the launcher can answer rather than one the orchestrator must
   reconstruct.

Two further requirements come from the harness rather than the launcher, and the
evidence is unambiguous about both:

8. **The orchestrator cannot wake itself.** Every one of the nine idle gaps ends
   with a user message or a background-task notification; there is no third kind
   of resumption in the transcript. A detached run with no task id is therefore
   unobservable *in principle*, not merely in practice. Any launch the
   orchestrator wants to hear about must be either `run_in_background: true` or
   shadowed by a background watcher that the harness knows about. Both were
   available and both worked when used; the five times a watcher was armed, dead
   time stayed under two minutes.
9. **The observable must not be a pipe that buffers.** `| tail -5` is the failure
   in miniature: it converts a stream into a single deferred write. Whatever the
   orchestrator polls has to be a file the launcher updates as it goes.

---

## 9. Surprises, and what could not be determined

**Surprises.**

- *The owner's foreground `!` runs were the well-behaved ones.* Timeout promotion
  gave them a task id, so the harness tracked them and one produced the
  notification at 05:21:31. The orchestrator's own `nohup` launches — the
  deliberate, "correct-looking" backgrounding — were the invisible ones. The
  mechanism that looked like a limitation was the one that worked.
- *The orchestrator had the right pattern and dropped it.* A guardrail hook
  explicitly taught it the `until`-loop watcher at 05:08:39; it used the pattern
  five times, every use produced a dead time under two minutes, and then it
  stopped at 16:50:53 and launched six more runs bare. Nothing in the transcript
  explains the change.
- *`| tail -5` was inert and the launches were still the best in the session.*
  The 4-second dead times for u05 and u08 came entirely from the task id, not from
  the pipe. The orchestrator's own summary called them "monitored" and credited
  the wrong half of what it had done.
- *The `ls` on the same line already held the answer.* At 15:59:55 the probe
  printed `u3b RUNNING` and, in the next line of the same output,
  `4935 Aug 21 07:20 logs/codex/u3b_scrublet_threshold.final.md`. The truth and
  the falsehood arrived together and the falsehood won for eighteen seconds — and
  had won for the preceding four hours.
- *`grep -c 'git (add|commit|…)'` was run after every single unit.* The guardrail
  verification the house rules demand was executed consistently and without
  prompting. Observation discipline was not uniformly absent; it was absent
  precisely where the wrapper gave it nothing to observe.

**Could not be determined.**

- *Whether concurrent runs of different units wrote the same paths.* The run logs
  distinguish neither reads from writes nor when anything happened, so path
  mentions cannot be promoted to concurrent writes. Only the u1_build duplicate is
  provable, and that from the wrapper's own redirection rather than from any log.
- *Why the watcher-loop pattern was abandoned after 16:50:53.* No reasoning about
  it appears in the transcript.
- *Whether the owner was present during the 05:41→11:33 and 12:04→15:59 gaps.*
  Absence of a message is not evidence of absence, so the 349.8-minute u2_droplet
  figure should be read as a measure of how long a finished run went unread, and
  not as a measure of inattention.
- *Why some background tasks produced a visible notification and others did not.*
  Notifications appear for `bnzii60k9`, `bj4zvpzmh`, `bjb9zf6uj`, `bx3ohteug`,
  `buxipoqgp` and `bapxzb86w`, but `bs3gpurlh`, `b2lwx21kq`, `b27c8jk84`,
  `b47cnd64g`, `bngzbq4qf` and `bloee0e6l` completed without one appearing in the
  transcript. Several of those completed while the orchestrator was mid-turn, and
  a `/compact` ran at 05:19:33, either of which could account for it. The
  mechanism is outside this evidence.
- *The outcome of the two runs still in flight.* `u08_annot_prep`'s relaunch
  exited at 19:32:30Z during this investigation; `u05b_pilot_fix` was still
  running at 19:33:53Z and was left alone.
