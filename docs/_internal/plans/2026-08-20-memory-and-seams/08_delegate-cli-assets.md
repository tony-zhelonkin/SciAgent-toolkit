# Phase 08 — delegation becomes a mechanism, not a habit

**Repo:** scio · **Two workers, disjoint files** · Evidence:
`docs/_internal/research/2026-08-21-delegation-invocation/`

## Why

`delegate-cli` is 416 lines of prose asserting CLI flags. Five verified failures:
it prescribes `--search`, which exists in **neither** codex 0.147.0 nor 0.149.0
(`SKILL.md:58`, `:246`) and therefore *caused* a field failure; its
backgrounding advice is a Claude Code tool parameter (`:168`) that does not
exist when the operator runs with `!`, so a run blocked a session for 1m42s;
it says nothing about concurrency, and two overlapping `u1_build` runs wrote the
same paths for ~2.5 minutes, making the worker's self-report unreliable; a long
inline invocation wrapped and split `-m gpt-5.6-sol`, so codex reported a missing
model value and bash tried to execute the model name; and it prescribes `/tmp`
paths as the convention.

Meta-Aging works only because the operator *habitually* uses file-backed prompts
through stdin. A habit that works is one keystroke from not working. JR-MC built
the right wrapper — untracked, inside a gitignored tree, unreconstructable from a
clone.

**Flag spelling belongs in an asset that probes, not in prose that asserts.**

## Worker 1 — the assets

Create `skills/delegate-cli/assets/probe.sh` and
`skills/delegate-cli/assets/launch.sh`. Bash, no `jq`/`yq`/`python3`, `chmod +x`,
`bash -n` clean, `set -eu` at minimum.

### `probe.sh`

Args: `--workdir DIR`, `--model MODEL`, `--sandbox MODE`, `--web`,
`--file-read-check`.

Emits shell-readable `KEY=VALUE` lines: `CODEX_VERSION`,
`HAS_OUTPUT_LAST_MESSAGE`, `HAS_STDIN_PROMPT`, `HAS_CONFIG`, `HAS_ENABLE`,
`HAS_BYPASS`, `HAS_SEARCH_FLAG`, `WEB_MODE=config_enable|unsupported`.

Derives every fact from `codex exec --help`. **Never `strings`** — an operator
burned rounds on that and it is not a capability API.

`--file-read-check`: write a temp file holding a nonce, ask codex to read it,
assert the nonce appears in the `-o` output. That is what proves tool use rather
than recall, and it is what caught the container's file-inspection failure.

Exit: `0` ok · `10` codex absent · `11` help/version failed · `12` a required
flag missing · `13` web requested but unsupported · `14` file-read check failed.

Refuses to install codex, edit config, change sandbox posture, or assume network.

### `launch.sh`

Args: `--unit NAME --workdir DIR --prompt FILE`, plus `--rules FILE`,
`--model MODEL`, `--effort LEVEL`, `--sandbox MODE`, `--bypass`, `--web`,
`--expect PATH` (repeatable), `--parallel-ok`, `--bg`.

- Run dir under `${TMPDIR:-/tmp}/scio-delegate/<unit>/<run-id>/`, holding
  `prompt.md`, `final.md`, `stream.log`, `status.tsv`, `pid`.
- Concatenate `--rules` then `--prompt` into `prompt.md`; feed codex via stdin
  with trailing `-`; final message to `-o`; noisy stream redirected separately.
  **The typed command must stay short** — that is the wrap fix.
- `--web` expands to `-c tools.web_search=true --enable web_search_request`, and
  only after the probe reports `WEB_MODE=config_enable`.
- **Lock:** atomic `mkdir` lock dir. Same unit already active → refuse. A
  *different* unit active → refuse unless `--parallel-ok`, which makes
  disjointness an explicit operator claim rather than an accident. Parallel
  workers are wanted; unnoticed duplicates are not.
- `--bg`: return after writing pid and paths; the child records final status and
  byte counts. Without it, foreground exits with codex's status — **unless**
  `final.md` is empty or an `--expect` path is missing, which is nonzero.
- Print byte counts for the rendered prompt, final message, stream log and every
  `--expect` path.
- Refuse: empty or missing prompt, unsupported web capability, a stale lock whose
  pid is live.
- Must **not**: run git write commands, install anything, mutate project config,
  create project-local scaffolds, use `pgrep`/`pkill` name matching, or infer
  approval for unsandboxed execution.

### Tests — `tests/test_delegate_cli_assets.sh`

Both assets exist, are executable, `bash -n` clean. `probe.sh` on a stub `codex`
on `PATH` reports `HAS_SEARCH_FLAG=0` and the right `WEB_MODE`; exits 10 when
codex is absent; exits 13 for `--web` against a stub lacking the config flags.
`launch.sh` refuses a missing prompt, refuses an empty prompt, refuses a second
run of the same unit while a lock is held, **accepts** a different unit with
`--parallel-ok`, exits nonzero when an `--expect` path is absent, and honours
`--bg` by returning promptly with a pid file. Use a stub `codex` script; never
invoke the real CLI.

## Worker 2 — the skill surgery

`skills/delegate-cli/SKILL.md` only. Delete or replace, by line:

| Lines | What | Action |
|---|---|---|
| 47, 50–52 | `gpt-5.5` as the mandatory recipe | drop the hardcoded model |
| **58, 246** | `--search` | **delete — exists in no live version** |
| 93 | flags as a current guarantee | keep only as dated evidence |
| 100–105, 139–140 | the `/tmp/<unit>` convention | replace with the assets |
| 168, 286 | `run_in_background` | delete — Claude-tool-specific |
| 240–251 | quick map of direct CLI snippets | cut hard |

New shape: *when to delegate → how to prepare a prompt → run `probe.sh` → run
`launch.sh` → review the result.* **Assets own flag spelling; prose owns
judgement.** Keep the agy parsing hazards (real, and no asset covers agy), the
prompt-construction rules, and the three implementer rules at 152–160 — those
earned their place.

Keep the delegation-specific verification checklist — re-run the tests, inspect
touched files, treat a submodule as its own repo, verify against disk. The
broader stance (*"delegation replaces the typing, never the verification"*) is
user-global and belongs to dev-env; do not restate it here.

Cite `scbio-docker/docs/ai-integration.md` for the bwrap seccomp *cause* rather
than restating it. Frontmatter description stays ≤350 chars.

## Gates — both workers

```bash
bash tests/run-all.sh                     # 62 → 63, 0 failing
bin/scio lint --check toolkit --strict    # clean
bash -n skills/delegate-cli/assets/*.sh
grep -n '\-\-search\|run_in_background' skills/delegate-cli/SKILL.md   # empty
```

## Do not

- Do not add a hook. Do not touch `lib/scio/` (phases 01/02/#37 own it).
- Do not create anything under a project's `docs/_internal/`.
- Do not invoke the real `codex` from a test.
- Do not vendor JR-MC's wrapper verbatim — generalize it. It hardcodes one
  project root, one model, and has neither lock nor background mode.
