---
name: delegate-cli
description: >-
  Hand off or consult work to the local `codex` (OpenAI) and `agy` (Gemini) CLIs in headless mode.
  Trigger when Anton says "implement/do X via codex/agy", "consult codex/agy", "have codex/agy do/write
  X", "have agy web-research X", or similar. Tells you the exact flags, the working models (codex
  `-codex` variants are REJECTED on this account — must pass `-m gpt-5.5`), and that this vault is NOT
  a git repo so codex needs `--skip-git-repo-check`. Also covers the Linux devcontainer, where
  codex's sandbox is unavailable and `agy --model` is ignored.
license: MIT
metadata:
  scope: concept
  requires: []
  skill-author: SciAgent-toolkit
  last-reviewed: 2026-07-28
  category: workflow
  tier: standard
  tags:
  - tooling
  - planning
  complementary-skills:
  - reasoning-trace
  - architecture-first-dev
---

# Delegating to codex & agy (headless)

Two peer coding CLIs. Run them non-interactively, capture their output, relay the result.
Pick by strength: **codex** = OpenAI GPT-5.x, sharp implementer/reviewer. **agy** = Gemini 3.x, large
context + strong web research.

**Two environments, different flags.** The macOS vault is the default described below. The Linux
devcontainers (`/workspaces/<project>`) differ on sandboxing, model selection, and git — see
`§Linux devcontainer` before delegating there.

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

```
agy -p [--model "Gemini 3.1 Pro (High)"] [--sandbox] [--add-dir DIR] "PROMPT"
```

- **`-p` (`--print`)** runs one prompt and prints a clean answer to stdout — no header to strip.
- Default model **Gemini 3.1 Pro**; good default for web research and big-context reads. Other names from
  `agy models` (e.g. `Gemini 3.5 Flash (High)`, `Claude Opus 4.6 (Thinking)`); pass via `--model "<exact
  display name>"`.
- `--sandbox` restricts the terminal; `--dangerously-skip-permissions` auto-approves; `-c`/`--continue`
  resumes the last conversation; `--print-timeout` default 5m.

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
never cites a file it was asked to read, suspect this first.

**agy — `--model` is ignored.** `agy models` advertises `gemini-3.1-pro-high` and the server
confirms the entitlement, yet every run logs
`Propagating selected model override to backend: label="Gemini 3.6 Flash (High)"` regardless of the
flag. No model key exists in `~/.gemini/antigravity-cli/settings.json` or
`~/.gemini/config/config.json`. Check which model you actually got:

```
grep 'Propagating selected model' ~/.gemini/antigravity-cli/cli.log | tail -2
```

Changing it likely needs one interactive `agy` session, after which headless runs inherit the
selection. **Verify the model before trusting a delegated result**, and say which model produced it
when relaying.

```
agy --print --model=<id> --effort=high --dangerously-skip-permissions \
    --print-timeout 25m --prompt "$(cat /tmp/<unit>_prompt.md)" > /tmp/<unit>.log 2>&1
```

`--print` buffers: the log stays 0 bytes until the run ends. Watch `git status` for progress.

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

**Problem:** In this devcontainer (Docker, seccomp filter active), `read-only` and `workspace-write` both fail completely. Every shell command and `apply_patch` exits with:
```
bwrap: No permissions to create a new namespace, likely because the kernel does not allow non-privileged user namespaces.
```
Codex uses bwrap (bundled at `~/.codex/packages/standalone/.../codex-resources/bwrap`) for ALL non-`danger-full-access` sandboxes. The Docker default seccomp profile blocks `clone(CLONE_NEWUSER)` — the syscall bwrap needs — even though `user.max_user_namespaces` is non-zero.

**Root cause:** `.devcontainer/devcontainer.json` was missing `--security-opt seccomp=unconfined`. **Fixed** (added to `runArgs`). Requires devcontainer rebuild to take effect.

**Workaround until rebuild:** Use `-s danger-full-access` — this skips bwrap entirely and runs without any sandbox. File writes work, shell commands work. **Anton's preference:** run `danger-full-access` commands himself via `! codex exec ...` in the Claude Code prompt, then tell me so I monitor artifact creation. The auto-mode classifier blocks me from invoking `danger-full-access` autonomously (it matches the `DANGEROUSLY_*` pattern and requires explicit user authorization).

**Pattern for Anton to run:** `! codex exec --skip-git-repo-check -m gpt-5.5 -s danger-full-access - < /abs/prompt.md` — Anton runs in the `!` prefix; I prepare the prompt file and monitor output.

**After rebuild:** All three modes (`read-only`, `workspace-write`, `danger-full-access`) should work normally.

## agy --add-dir quirk (observed 2026-06-29)

**Problem:** Passing `--add-dir DIR` BEFORE the prompt triggers agy's agentic mode — it treats the directory as a Claude Code project workspace and starts exploring it autonomously instead of executing the print-mode task.

**Correct order:** PROMPT must come immediately after `-p`, then flags after:
```bash
# CORRECT — prompt first, add-dir after
agy -p "$(cat prompt.txt)" --add-dir /path/to/dir --model "..." --print-timeout 10m

# WRONG — add-dir before prompt triggers agentic exploration
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
| agy web research | `agy -p "research …; cite sources"` |
