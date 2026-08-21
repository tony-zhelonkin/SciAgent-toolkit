# Consultation — reorganise `delegate-cli` around observability

Read-only. Produce a written recommendation. **Modify nothing.**

## Read first, in this order

1. `docs/_internal/research/2026-08-21-delegation-monitoring/TRACE_jr-mc-session.md`
   — a forensic trace of one 17-hour session, from the Claude Code transcript and
   the run logs. 651 lines. **This is the evidence base; do not re-derive it.**
2. `skills/delegate-cli/SKILL.md` (150 lines, current) and
   `skills/delegate-cli/assets/{probe.sh,launch.sh}` — what ships today.
3. `docs/_internal/research/2026-08-21-delegation-invocation/CONSULT_skill_design.report.md`
   — the consult that produced those assets. Its contracts are what you are
   revising, not starting from.

## What the trace established

**The failure was delivery, not design.** The project ran the pre-asset 416-line
skill. `assets/` did not exist in its pinned checkout. The commit that added
`probe.sh`/`launch.sh` existed in the toolkit for **13 h 25 min** before the owner
complained, and nothing re-pins a vendored submodule mid-session.

**Dead time is entirely in launch-to-observation.** Worst case 349.8 minutes;
219.7 minutes with the owner demonstrably present. Observation-to-action was
negligible everywhere — 4 s, 4 s, 30 s. The split is mechanical: **every** dead
time above six minutes was a `nohup … &` launch with no watcher; **every** one
under two minutes was a run the harness knew about.

**Prose warnings did not work.** The old skill documented the `pgrep -f`
self-match trap at line 298 and the `pkill -f` trap at line 342, both marked
verified. The session reproduced both — three false `RUNNING` verdicts from a
pattern that matched the polling shell itself, and a `pkill -f` that returned
exit 144, the exact code the skill named. In one case the truth and the falsehood
arrived in the same command output and the falsehood won for four hours.

**`| tail -5` was inert.** Measured: 0 bytes observable after ten minutes while
the run log had grown to 425,956 bytes; then 199 bytes at exit, which were the
wrapper's own three banner lines plus a decorative `exited 0`.

**The wrapper's exit code is decorative.** `set -euo pipefail` plus
`echo "==> $UNIT exited $?"` as the last line means a failure aborts before the
echo, so that line can only ever print `0`.

**The orchestrator cannot wake itself.** All nine idle gaps end with a user
message or a harness background-task notification. There is no third resumption
mechanism in the transcript. A detached run with no harness task id is therefore
unobservable *in principle*.

**It had the right pattern and dropped it.** A watcher loop was used five times,
each time with dead time under two minutes, then abandoned at 16:50:53 with no
reasoning recorded, and six further runs went out bare.

**Observation discipline was not the problem.** The git-write guardrail grep ran
after every single unit, unprompted. Observation was absent precisely where the
wrapper gave it nothing to observe.

## The two design gaps this exposes in what we already ship

State plainly whether you agree, and say what changes.

**(a) `launch.sh --bg` reproduces the `nohup` failure.** It detaches, writes a
pid and a `status.tsv`, and returns. Nothing then wakes the orchestrator. By the
trace's requirement 8 that makes `--bg` insufficient on its own for a Claude
orchestrator: the launch must be either run through the harness's own background
mechanism (which holds a task id and notifies) or shadowed by a watcher the
harness knows about. Today's SKILL.md does not say this, and `--bg` reads like
the answer.

**(b) There is no cheap progress read.** `status.tsv` is written once, after the
child reaps. During a run, the only signal is the stream log's growing size,
which the trace shows the session polling with `stat -c %s`. There is no
last-activity timestamp, so a hung run and a slow run are indistinguishable — the
session had no way to tell them apart and never tried.

## Questions

1. **The one-command answer.** Specify the command that answers "is unit X still
   running, and how far along" from a single cheap read, without touching the
   process table. What does it print, in what format, and where does that state
   come from? Requirements 1–7 in the trace's §8 are the target. Is this a new
   `status.sh` asset, a `--status` mode on `launch.sh`, or something else — and
   why that shape?
2. **Who backgrounds, and how the skill says it.** Resolve gap (a). The owner
   sometimes runs launches himself behind `!`; the orchestrator sometimes launches
   them. The trace found the owner's foreground `!` runs were the *well-behaved*
   ones, because timeout promotion gave them a harness task id, while the
   orchestrator's deliberate `nohup` runs were invisible. What should the skill
   prescribe for each caller, and should `--bg` survive at all in its current
   form?
3. **Liveness and progress on disk.** What exactly should the launcher write
   while the child runs, so that a poll is one file read? Consider a heartbeat
   line, the stream's size and mtime, a terminal state distinguishable from
   "never started", and an exit code that is real rather than decorative. Bash
   only — no `jq`, no `yq`, no `python3`. Say what it must NOT write, and keep the
   file count small.
4. **Hung versus slow.** Can a bash predicate distinguish them without a policy
   knob the owner has to tune? The owner has rejected configuration knobs
   carrying intent, from experience: they produced churn and variables nobody
   could justify. If the honest answer is "report the numbers and let the reader
   judge", say so.
5. **What the skill should stop saying.** The old skill's traps were documented
   and reproduced anyway. Which of today's 150 lines are prose that an asset
   should own instead? Be specific by line. Deleting counts as progress.
6. **Delivery.** The mechanism existed and did not arrive. Beyond the fleet sweep,
   is there anything the *skill or its assets* can do so that a session running an
   older pinned copy discovers it is stale? Constraints: no new hook, bash only,
   nothing that phones home, and `lint` is the only enforcement surface. If the
   answer is "nothing, this is the sweep's job", say that plainly rather than
   inventing a version check.
7. **Concurrency, given the evidence.** The trace confirms one duplicate
   (`u1_build`, 2 min 48 s) and one dependency violation launched as parallelism.
   Today's lock refuses a duplicate unit and requires `--parallel-ok` beside a
   different live unit. Is that the right shape, and does anything need to
   express *dependency* rather than disjointness? Note the project expressed
   dependency order in a prose runbook, which is not a mechanism.

## Constraints — settled, do not reopen

- **No new hook.** Not Claude, not Codex, not git.
- **Bash only** in anything scio ships: no `jq`, `yq`, `python3`.
- **No hand-maintained metadata** and no intent-carrying config knobs.
- **Source shape equals mount shape** — anything added under
  `skills/delegate-cli/` is mounted into every project, so `assets/` is reachable
  with no new delivery mechanism and also unavoidable.
- The skill is 150 lines and should not grow. Assets own mechanics; prose owns
  judgement.
- Keep the agy parsing hazards and the delegation-specific verification
  checklist — both earned their place.

## Deliverable

1. The status-read contract (Q1, Q3): arguments, output format, exit codes, files.
2. The backgrounding prescription per caller, and the fate of `--bg` (Q2).
3. Your answer on hung-versus-slow (Q4).
4. Line-level deletions from `SKILL.md` (Q5).
5. The delivery answer (Q6).
6. The concurrency verdict (Q7).
7. Where you disagree with this brief or with the trace's conclusions.

Ground everything in the trace. Distinguish verified from inferred.
