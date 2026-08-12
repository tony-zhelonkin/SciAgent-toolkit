---
name: delegate-cli
description: >-
  Hand off or consult work to the local `codex` (OpenAI) and `agy` (Gemini)
  CLIs in headless mode. Trigger when Anton says "implement/do X via
  codex/agy", "consult codex/agy", "have codex/agy do/write X", "have agy
  web-research X", or similar. Covers headless flags, model/account quirks,
  sandbox fallbacks, installation, and multi-agent fan-out.
license: MIT
compatibility: "sciagent-scaffold: docs/_internal/reasoning/"
---

# Delegating to codex & agy (headless)

These peer coding CLIs may be available on the current host. Run them non-interactively, capture their
output, and relay the result. The account/model and devcontainer notes below are a verified **local
profile**, not assumptions to make on other operating systems or environments.
Pick by strength: **codex** = OpenAI GPT-5.x, sharp implementer/reviewer. **agy** = Gemini 3.x, large
context + strong web research.

**Two environments, different flags.** The macOS vault is the default described below. The Linux
devcontainers (`/workspaces/<project>`) differ on sandboxing, model selection, and git — see
`§Linux devcontainer` before delegating there.

## Preflight and installation

Check the requested CLI before preparing a long run (`command -v codex` / `command -v agy` on POSIX;
use the shell's equivalent executable lookup elsewhere). If it is missing, do not silently substitute a
different agent. Tell the user what is unavailable and offer the official installer:

```bash
# OpenAI Codex CLI
curl -fsSL https://chatgpt.com/codex/install.sh | sh

# Google Antigravity CLI (`agy`)
curl -fsSL https://antigravity.google/cli/install.sh | bash
```

Fetching and executing a remote installer changes the user's system and uses the network: ask before
running it. Afterward, verify with `codex --version` or `agy --version`; if the executable is not on
`PATH`, ask the user to reload/restart their shell. Follow the CLI's own sign-in flow and never ask the
user to paste credentials into chat. The commands above are POSIX-shell examples; give PowerShell or
other OS-appropriate instructions when that is the user's environment.

## codex (OpenAI) — `codex exec`

```
codex exec --skip-git-repo-check -m gpt-5.5 -s <sandbox> [-o OUT.txt] "PROMPT" < /dev/null
```

- **`-m gpt-5.5` is mandatory.** Verified on Anton's ChatGPT account: `gpt-5.5` and `gpt-5.4` work;
  **every `-codex` variant fails 400** ("not supported when using Codex with a ChatGPT account"), and the
  config default `gpt-5.2-codex` is one of them — so omitting `-m` errors out. Use `gpt-5.5` (newest).
- **`--skip-git-repo-check`** — this vault is not a git repo; codex refuses to run without it.
- **`< /dev/null`** — headless codex otherwise blocks on "Reading additional input from stdin…".
- **Sandbox `-s`** picks what it may touch: `read-only` (consult/review/plan), `workspace-write` (let it
  edit files here), `danger-full-access` (avoid). `exec` runs with approval `never` already.
- **Clean output:** `-o OUT.txt` writes only the final message to a file (stdout is noisy — header,
  reasoning, token count). Read that file. `--json` for event stream, `--search` to enable web search,
  `-C DIR`/`--add-dir DIR` to set/extend the writable root.

## agy (Gemini) — `agy -p`

Flag surface verified against agy 1.1.10.

```
agy -p "PROMPT" [--model "Gemini 3.1 Pro (High)"] [--dangerously-skip-permissions] \
    [--add-dir DIR] [--print-timeout 25m] > OUT.md 2>&1
```

- **`-p`, `--print`, and `--prompt` are one string flag** (`-p` short, `--prompt` an alias); the last
  occurrence wins, and it prints a clean answer to stdout with no header to strip. It consumes the
  next token as its value, so a bare `--print` in front of another flag eats that flag — in
  `--print --dangerously-skip-permissions … --prompt "…"` the permissions flag never takes effect.
  Give the prompt as the flag's value and put the remaining flags after it.
- A bare prompt with no `-p` starts interactive mode, which dies headless on `/dev/tty`. Flag parsing
  stops at the first bare argument, so flags after one are silently dropped.
- **`--model`** accepts either the display label (`"Gemini 3.1 Pro (High)"`) or the slug printed by
  `agy models` (`gemini-3.1-pro-high`). An unrecognized value exits 1 before the run and prints the
  label list. The session default is the `model` key in `~/.gemini/antigravity-cli/settings.json`,
  which an interactive session rewrites.
- **Omit `--effort`.** Every model in `agy models` bakes its effort into its name and rejects the
  flag: `invalid model selection (--model "Gemini 3.1 Pro (High)" --effort "high"): --effort is not
  supported for model "Gemini 3.1 Pro (High)"`, exit 1 before the run. Pick effort by picking the
  model.
- **`--print-timeout`** takes a Go duration with a unit — `25m`, `300s`; default 5m. A bare number
  exits 2 with `missing unit in duration "30"`.
- `--sandbox` restricts the terminal; `-c`/`--continue` resumes the last conversation.
- Print mode buffers: the redirect target stays 0 bytes until the run ends, then the whole answer
  lands at once.

## Linux devcontainer (`/workspaces/<project>`)

Verified 2026-07-28 in the Meta-Aging containers. Four differences, each found by execution.

**codex — the sandbox does not work.** The container disallows unprivileged user namespaces, so
codex's bwrap sandbox fails to start under any `-s` value. Use
`--dangerously-bypass-approvals-and-sandbox`, and get the safety back by writing the prohibitions
into the prompt and verifying afterwards:

```
codex exec -m gpt-5.6-sol -c model_reasoning_effort="high" \
  --dangerously-bypass-approvals-and-sandbox \
  -C /workspaces/<project> --skip-git-repo-check \
  -o /tmp/<unit>_last.md < /tmp/<unit>_prompt.md > /tmp/<unit>_run.log 2>&1
```

The projects **are** git repos here, so `--skip-git-repo-check` is belt-and-braces rather than
required. Because nothing sandboxes git, forbid write verbs in the prompt and confirm afterwards
with `grep -cE '\bgit (add|commit|push|checkout|restore|reset|stash|branch)' <log>` plus a
`rev-parse HEAD` / `reflog` check in **every** repo the agent could reach — a submodule is its own
repo.

**codex — usage limits kill runs instantly.** A limit-exhausted launch exits 1 within seconds having
written nothing. Agents also die mid-write on session limits, usually *after* verification. Always
check what landed on disk before relaunching, and remember a gitignored tree shows nothing in
`git status`.

**agy — headless auto-denies every tool.** Without `--dangerously-skip-permissions`, tool calls are
refused because headless mode cannot prompt, and the model answers from its own knowledge instead.
The failure is quiet: you get a fluent, plausible, entirely un-grounded reply. If an `agy` answer
never cites a file it was asked to read, suspect this first. With the flag, print mode logs
`--dangerously-skip-permissions set, auto-approving all tool permissions` — grep that line to
confirm the flag survived parsing.

**agy — read the model back after the run.** Each invocation writes a fresh
`~/.gemini/antigravity-cli/log/cli-<timestamp>.log` and repoints the `cli.log` symlink at it, so a
grep taken before a run reports the previous run's model. Read it afterwards:

```
grep -o 'label="[^"]*"' ~/.gemini/antigravity-cli/cli.log | tail -1
```

The full line is `Propagating selected model override to backend: label="…"`. It names the model the
session was configured with — `--model` when the flag parsed, the `settings.json` default otherwise,
which is also what a run that exits on an invalid `--model` logs. Say which model produced a
delegated result when relaying it.

```
agy -p "$(cat /tmp/<unit>_prompt.md)" --model "Gemini 3.1 Pro (High)" \
    --dangerously-skip-permissions --print-timeout 25m > /tmp/<unit>.log 2>&1
```

## Sizing the hand-off to the model

Match task size to the implementer. A unit sized for codex is too coarse for Flash.

- **Flash-class** (`gemini-3.6-flash-*`) — fast and capable, and it strays on style and on subtle
  assertion strength. Decompose further, state goals tightly, and fan out **more often** with
  verification between hops rather than handing over one large autonomous unit.
- **Pro / GPT-5.x-class** — can hold a multi-part unit with thirty discrete items.

Three rules that hold for every implementer, learned the expensive way:

1. **Put exact numeric targets in the prompt** (test counts, ledger sizes) and say *report a
   deviation rather than editing the artifact to match*. This is what surfaces an honest
   disagreement instead of a silent adjustment.
2. **Forbid weakening or deleting a test to make a change fit.** When an implementer relaxes an
   assertion, the relaxation is the finding. Ask explicitly for the disposition of every test it
   removed, by name.
3. **Review is not optional.** An independent reviewer on a stronger model has caught a Tier-1 in
   work whose own suite was green, more than once. Delegation replaces the typing, never the
   verification: re-run the tests yourself and check what the agent actually touched.

## How to drive them from here

1. **Compose a self-contained prompt** — they don't share our chat. State the task, the files/paths
   (absolute, quoted — spaces), and the exact deliverable shape.
2. **Run in background** if it's real implementation (`run_in_background: true`) — these take minutes.
3. **Relay, don't dump.** Summarize what came back; show their diff/answer. For codex edits in
   `workspace-write`, review the changed files yourself before telling Anton it's done.
4. **Match tool to ask:** implementation/code-review → codex; web research / very large context → agy.
   When Anton names one, use that one.

## codex sandbox constraints (workspace-write, verified firsthand)

**The clean pattern: codex writes files; the orchestrator (you, unsandboxed) builds and commits.**

- **No git writes.** `.git` is readable but not writable — `git add` fails with `Unable to create .../.git/index.lock: Operation not permitted`. Never put commit instructions in a codex prompt; they'll fail and waste a turn. After codex finishes editing files, `git add` / `git commit` in the main session yourself.
- **No build/probe steps that touch the kernel.** `quarto render` (and likely similar tools) does a macOS `sysctl` architecture probe that the seatbelt sandbox blocks → dies with "unrecognized architecture" before loading the project. Let codex write the files; run `quarto render` yourself in the main session. Generalize: any tool that probes the machine (arch detection, hardware queries) needs to run outside the sandbox.
- **Network is off.** No outbound connections under `workspace-write` — CDN, npm, font fetches all fail. Seed any external assets (fonts, data files) into the workspace before launching codex.
- **Writable roots are printed in the codex header**, e.g. `sandbox: workspace-write [workdir, /tmp, $TMPDIR]`. Reads elsewhere on disk are fine — you can point codex at absolute paths outside the workspace to read specs or docs.
- **Images via `-i`.** `-i IMAGE` attaches images to the prompt; gpt-5.5 is vision-capable. Repeat per file. Useful for design-to-code (e.g., feed Figma board screenshots).
- **Long prompts: use stdin.** `codex exec … - < /abs/prompt.md` (the `-` reads from stdin). Avoids shell-escaping hell for multi-paragraph prompts. Use `< /dev/null` only when the full prompt is in the inline arg and you just want to stop it blocking on stdin.

## devcontainer bwrap restriction (observed 2026-06-29)

This is one observed Linux/devcontainer failure mode, not a universal diagnosis. On another host,
preserve and report the exact sandbox error (for example namespace/bwrap, seatbelt, or container-policy
denials) and offer that platform's supported environment fix rather than prescribing Linux kernel or
Docker settings blindly.

**Problem:** In this devcontainer (Docker, seccomp filter active), `read-only` and `workspace-write` both fail completely. Every shell command and `apply_patch` exits with:
```
bwrap: No permissions to create a new namespace, likely because the kernel does not allow non-privileged user namespaces.
```
Codex uses bwrap (bundled at `~/.codex/packages/standalone/.../codex-resources/bwrap`) for ALL non-`danger-full-access` sandboxes. The Docker default seccomp profile blocks `clone(CLONE_NEWUSER)` — the syscall bwrap needs — even though `user.max_user_namespaces` is non-zero.

**What restores it:** the container must permit unprivileged user namespaces — `security_opt: [seccomp=unconfined]` on the compose service (or `--security-opt seccomp=unconfined` in `runArgs`), then a rebuild. The Meta-Aging containers run without it, so the restriction is live: `unshare -U true` exits `Operation not permitted` while `user.max_user_namespaces` reads non-zero.

**Workaround while the restriction stands:** Use `-s danger-full-access` — this skips bwrap entirely and runs without any sandbox. File writes work, shell commands work. **Anton's preference:** run `danger-full-access` commands himself via `! codex exec ...` in the Claude Code prompt, then tell me so I monitor artifact creation. The auto-mode classifier blocks me from invoking `danger-full-access` autonomously (it matches the `DANGEROUSLY_*` pattern and requires explicit user authorization).

A passing `read-only` probe does **not** prove `workspace-write` will run: read-only answers pure Q&A but dies the moment it inspects a file, so it clears the probe and fails the job.

Asking to “use codex” authorizes Codex, **not** disabling its sandbox. Never infer approval for
`danger-full-access`; require the user to run the exact command or explicitly approve that unsandboxed
retry. Do not change host security settings, container privileges, kernel settings, or sandbox mode
without approval. When unsandboxed execution is approved, prepend the allowed working directory/files
and re-forbid git, network, credentials/secrets, dependency installation, destructive commands, and
out-of-scope writes. Prompt constraints reduce risk but do not recreate a sandbox.

**Pattern for Anton to run:** `! codex exec --skip-git-repo-check -m gpt-5.5 -s danger-full-access - < /abs/prompt.md` — Anton runs in the `!` prefix; I prepare the prompt file and monitor output.

**Once the container permits user namespaces:** all three modes (`read-only`, `workspace-write`, `danger-full-access`) work normally.

## agy flag order (verified on agy 1.1.10)

**The prompt is `-p`'s value, and the first bare argument ends flag parsing.** A flag placed between
`-p` and the prompt therefore becomes the prompt, the directory that follows it terminates parsing,
and the real prompt never arrives — the run was observed to explore the added directory agentically
instead of executing the print-mode task.

```bash
# CORRECT — prompt is -p's value, flags follow
agy -p "$(cat prompt.txt)" --add-dir /path/to/dir --model "Gemini 3.1 Pro (High)" --print-timeout 10m

# WRONG — -p swallows --add-dir, /path/to/dir ends parsing, prompt.txt is dropped
agy -p --add-dir /path/to/dir "$(cat prompt.txt)"
```

**Also wrong:** Putting "Write your output to /path/to/file.md" in the prompt text also triggers agy's agentic file-write mode. Always capture stdout via shell redirect (`> output.md`) instead.

**Working pattern (from orchestrate.sh):**
```bash
agy -p "$(cat "$SCRATCHPAD/prompt.txt")" \
    --add-dir "$RESEARCH_DIR" \
    --print-timeout 10m \
    > "$RESEARCH_DIR/output.md" 2>&1 &
```

## Quick map

| Ask | Command skeleton |
|---|---|
| Consult codex (read-only) | `codex exec --skip-git-repo-check -m gpt-5.5 -s read-only -o OUT.txt "…" < /dev/null` |
| Have codex implement here | `codex exec --skip-git-repo-check -m gpt-5.5 -s workspace-write "…" < /dev/null` |
| codex web search | add `--search` |
| Long prompt via file | `codex exec … -m gpt-5.5 -s workspace-write - < /abs/prompt.md` |
| Attach image(s) | add `-i /abs/image.png` (repeat per file) |
| Consult agy | `agy -p "…"` |
| agy web research | `agy -p "research …; cite sources" --dangerously-skip-permissions` |
| Pick the agy model | `agy -p "…" --model "Gemini 3.1 Pro (High)"` (no `--effort`) |

## Fan-out orchestration defaults (multi-agent delegation)

When fanning out real work across these CLIs plus Claude subagents (verified 2026-07-14 on the JIA
heat-vs-hypoxia stage), these defaults hold up:

- **Keep fan-out bounded and independent.** Default to two or three workers with one concrete
  deliverable each. Parallelize read-only reviews/research freely within policy, but never let writers
  edit the same file or overlapping directories. Give each writer an exclusive path/worktree or ask
  for a patch-only deliverable, then integrate serially.
- **Separate permissions.** Network-enabled research, workspace writes, and unsandboxed execution are
  distinct authorizations; approval for one does not imply either of the others. Give every run unique
  prompt/output paths and do not reuse `agy --continue` across independent workers.

- **Split by model strength, and keep reviewer ≠ author.** codex (GPT) implements → a *different*
  model family (a Claude/Opus subagent) reviews the code → agy (Gemini) does the domain-knowledge /
  large-context interpretation (e.g. classifying gene-set biology) → the orchestrator integrates,
  RE-RUNS to verify, and commits. Cross-model review only adds signal when the reviewer did not write
  the code — never have one Claude both author and review.
- **Probe before the long job.** A tiny read-only call (`codex exec … -s read-only -o probe.txt "say
  CODEX_OK" < /dev/null`, or `agy -p "reply OK"`) confirms the CLI runs here before you launch a
  multi-minute unattended job. Caveat: in the devcontainer a passing `read-only` probe does NOT prove
  `workspace-write` works (see the bwrap section) — read-only dies as soon as it shells out.
- **Self-contained prompts, composed from files, piped via stdin.** These CLIs don't see your chat:
  state the task, absolute quoted paths, the infra to reuse, the exact output contract, the house
  rules, and where to save reasoning. Keep a reusable spec file and prepend a small per-run note:
  `cat prepend.md spec.md | codex exec … -`.
- **Clean capture + orchestrator owns the artifacts.** codex: `-o OUT.txt` (final message only) and
  redirect the noisy stream separately (`> stream.log 2>&1`). agy: `-p` prints clean to stdout →
  redirect to a file. (Superseded 2026-07-28: agy *can* be asked to write and edit files — see the
  worker-edits/orchestrator-commits section. What hangs a headless run is a **shell command** outside
  `permissions.allow`, not file writing.) For
  stateful research, have the orchestrator save each agent's stdout under `docs/_internal/reasoning/`
  (or tell codex, which can write, to save there) so the multi-model chain is traceable and replayable.
- **Background the slow ones** (`run_in_background`), each with its own scratch OUT path; read the OUT
  file on the completion notification, then relay — don't dump the raw stream.
- **codex must never run git.** Under `workspace-write`, Git metadata writes are denied and can waste
  the run while leaving worktree edits behind; under `danger-full-access` it genuinely could commit or
  push. Put "do NOT run git" in the prompt and commit yourself after reviewing
  the diff. Under `danger-full-access`, re-assert ALL guardrails in the prompt (no git, no network,
  scope to one dir) — the sandbox is no longer enforcing them.
- **Verify, don't trust.** Headless codex prints `reasoning effort: none` by default; for hard tasks
  that is thin. Raise it in-flight with `-c model_reasoning_effort=high` (verified: the header then
  prints `reasoning effort: high`), or use a higher-reasoning model (`-m gpt-5.6-sol`). Always re-run
  the produced artifact yourself before believing its numbers or committing.

## Waiting on background workers — the self-match trap (verified 2026-07-27)

The obvious way to block until a fan-out finishes is a `pgrep` poll loop:

```bash
# BROKEN — never exits
while pgrep -f "codex exec --skip-git-repo-check" > /dev/null; do sleep 30; done
echo "ALL WORKERS EXITED"
```

**The loop's own command line contains the pattern**, so `pgrep -f` matches the waiter's shell and
the condition stays true after every worker has exited. It runs until it is killed or times out, and
the timeout surfaces as a *task failure* that looks like the workers died. They did not. Check
`out.txt`/`stream.log` before concluding anything about the workers.

Fixes, cheapest first:

- **Wait on PIDs, not on a name.** Capture `$!` at launch and `wait "${pids[@]}"`. Only works if the
  waiter is the launching shell — across shells, poll the PIDs with `kill -0 "$pid"`.
- **Wait on the artifact.** The workers' output files are the real completion signal:
  `while [ ! -s w1.out.txt ] || [ ! -s w2.out.txt ]; do sleep 30; done`. Robust across shells, and it
  waits for what you actually care about.
- **Exclude self** if you must match by name: `pgrep -f "codex exec" | grep -v "^$$\$"`, or break the
  literal so the waiter's own cmdline does not contain it (`"codex e""xec"`).

Generalizes to any `pgrep`/`ps | grep` poll over a pattern that appears in the polling command.

### Killing a parent orphans its forked workers (verified 2026-07-28)

The mirror image of the same trap, on the way out. Stopping a fan-out by killing the PID you
captured at launch kills only the parent. Anything it forked — R's `parallel::mclapply`, Python's
`multiprocessing`, a `&`-backgrounded loop — is reparented to init and **keeps running at full
tilt**. Nothing reports this: the log stops growing, the parent is gone, and the job looks stopped.

It surfaced here as a later run taking 20 minutes for work that had taken 13, on a 72-core box
sitting at load 87. Two earlier killed runs had left 105 orphaned workers between them, each still
holding a core.

- **Kill the group, not the process:** launch under `setsid`, then `kill -- -"$PGID"`.
- **Audit before assuming a job is stopped:** `ps -o pid=,ppid=,cmd= -C R | awk '$2==1'` lists
  exactly the orphans. Reparenting to `ppid 1` is the signature.
- **Never `pkill -f <script>` to clean up** — same self-match problem, and here it killed the
  invoking shell (exit 144). Filter `ps` output by `ppid == 1` and pipe PIDs to `kill`, or break the
  literal with a bracket (`'probe[4]\.R'`) so the pattern cannot match its own command line.

## Two failure modes in how the *prompt* is written (verified 2026-07-27)

Both were the orchestrator's fault, not the CLI's, and both are cheap to avoid.

- **A wrong acceptance test gets satisfied, not reported.** An acceptance command that was subtly
  wrong (it invoked an import from the wrong working directory) was met by codex *creating a shim
  file* so the command would pass, rather than saying the test was wrong. It disclosed the new file
  under "what changed", so this is not deception — but the artifact was pure test-satisfying scaffold
  and had to be deleted. **Run every acceptance command yourself before putting it in the prompt**,
  and add a standing line: *"if an acceptance command here is wrong, report that instead of changing
  the repository to satisfy it."*
- **The acceptance test sets the blast radius.** A criterion phrased as a repo-wide
  `grep -rn <banned-string> <root>/` licensed edits anywhere in that root, and the worker duly
  rewrote an unrelated notebook to clear the grep. The work was defensible, but it was not the task,
  and it silently carried a stale number into prose nobody had reviewed. **Scope the verification
  command to the paths you scoped the task to.**

Corollary for both: give expected numbers as a **check with an explicit instruction to report a
mismatch rather than reconcile to it**. Workers otherwise treat a stated number as the target and
quietly make the artifact agree.

## The optimal split: the worker edits, the orchestrator commits (verified 2026-07-28)

**agy edits and regenerates; it has no git write permission, so the orchestrator reviews against disk
and commits.** That is the right split rather than a limitation — in a two-task pass, *both* tasks
needed a correction caught in review that the model had reported as clean.

Why it is a feature and not a workaround:

- **Withholding git write is not a downgrade, it is where the review gate lives.** A worker that can
  commit has already published before anyone looked. One that cannot leaves its work in the worktree,
  where `git status --porcelain <scope>` is a complete and cheap statement of what it touched. Give
  the worker read-only git (`git diff`, `git log`, `git show`, `git ls-files`) so it can verify its
  own work, and nothing that writes.
- **The worker's self-check has a blind spot its report cannot expose.** In the verified case a
  caption was rewritten from two lines to three; the numbers were right, the linter passed `0 hard`,
  and the banned string was gone — every check the worker could run said clean. The third line grew
  upward into the title and collided with it, because the text element was anchored `va="bottom"`.
  **Nothing in a text-level report can reveal a layout regression.** Whoever reviews must open the
  regenerated artifact, not read about it.
- **Generalize that:** after any worker edit that changes the *length* of rendered text, re-render and
  look at the image. Same for row counts in a table that has a fixed-height panel, and label text that
  feeds an axis.

Practical mechanics, both verified:

- **Small decomposed tasks beat one broad implementation pass.** At one-file / one-edit granularity
  the orchestrator can open every artifact the worker touched. Across a six-item pass it cannot, and
  the review degrades to reading the worker's own summary — which is exactly the signal that fails.
- **Do the locating yourself and hand over an exact edit.** Reading the target artifact *before*
  writing the prompt changed what the correct fix was in the verified case. A prompt that names
  `file:line`, the exact replacement text, and an explicit do-not-touch list gets a 2-line diff; a
  prompt that says "find and fix" gets a search, a judgment call, and a wider blast radius.
- **Plant a trap in the do-not-touch list.** A nearby string that *looks* like the defect but is
  correct (same wrong-looking number, different subject) tests whether the worker scoped its edit or
  pattern-matched. Cheap to include, and it either passes silently or catches a real over-reach.
- **agy writes files fine — but `cd` is auto-denied in headless mode and aborts the whole run.**
  Its permissions live in `~/.gemini/antigravity-cli/settings.json` under `permissions.allow` as
  `command(<prefix>)` rules. Headless (`-p`) cannot prompt, so any command outside the list is denied
  and the run ends with *"no output produced — a tool required the command permission"* and nothing
  else. Add the prefixes the task needs (`cd`, `mkdir`, `quarto`, …); deliberately omit `rm`, `mv`,
  `bash`, `sh`, `pip`, `curl`, `wget`, and every `git` write verb. Back the file up before editing it.
  Most generator scripts self-`chdir`, so check before granting `cd` at all.
- **Probe writes separately from reads.** A passing read-only probe says nothing about whether the
  worker can edit: write a throwaway file under a `_scratch/` path and read it back.

---

## See also

- `reasoning-trace`
- `architecture-first-dev`
