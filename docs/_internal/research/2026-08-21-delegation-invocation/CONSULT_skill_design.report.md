Verified: host `codex-cli-exec 0.147.0`; `codex exec --help` has `-c`, `--enable`, `-m`, `-s`, `-C`, `--add-dir`, `--skip-git-repo-check`, `-o`, `--json`, stdin via `-`, and the bypass flag. It has no `--search`.

**1. Capability Drift**
Use runtime probing, not prose flag claims, not a hand-maintained version table. Pinning Codex is a scbio-docker/devcontainer policy decision; scio should still probe because the user can upgrade, downgrade, or run outside that container.

Delete or replace these claims in [SKILL.md](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/skills/delegate-cli/SKILL.md:47):

- Lines 47, 50-52: hard-coded `gpt-5.5` as mandatory/newest recipe.
- Line 58 and line 246: `--search`.
- Line 93 as an active-current flag guarantee; keep it only as dated evidence.
- Lines 100-105 and 139-140: fixed `/tmp/<unit>` convention.
- Lines 168 and 286: `run_in_background`, which is Claude-tool-specific.
- Lines 240-251: quick map built from stale direct CLI snippets.

The skill should say: run `assets/probe.sh`, then use `assets/launch.sh`; if `--web` is requested, launch only when the probed capability supports the configured web-search mechanism. Verified current mechanism is JR-MC’s `-c tools.web_search=true --enable web_search_request`.

**2. Asset Contracts**
`probe.sh` contract:

- Args: `--workdir DIR`, `--model MODEL`, `--sandbox MODE`, optional `--web`, optional `--file-read-check`.
- Reports shell-readable facts: `CODEX_VERSION`, `HAS_OUTPUT_LAST_MESSAGE`, `HAS_STDIN_PROMPT`, `HAS_CONFIG`, `HAS_ENABLE`, `HAS_BYPASS`, `HAS_SEARCH_FLAG=0`, `WEB_MODE=config_enable|unsupported`.
- Runs `codex exec --help`; never uses `strings`.
- With `--file-read-check`, creates a temp nonce file, asks Codex to read it, and verifies the nonce appears in `-o` output. This catches the devcontainer bwrap/file-inspection failure.
- Exit codes: `0` ok; `10` codex missing; `11` help/version failed; `12` required flag missing; `13` web requested but unsupported; `14` file-read check failed.
- Refuses to install Codex, change config, change sandbox posture, or enable network by assumption.

`launch.sh` contract:

- Args: `--unit NAME --workdir DIR --prompt FILE`; optional `--rules FILE`, `--model MODEL`, `--effort high`, `--sandbox MODE|--bypass`, `--web`, `--expect PATH` repeatable, `--parallel-ok`, `--bg`.
- Builds a temp run dir under `${TMPDIR:-/tmp}/scio-delegate/.../<unit>/<run-id>/`.
- Concatenates `rules` plus prompt into a rendered prompt file, feeds Codex through stdin with `-`, writes final message via `-o`, and redirects noisy stream separately.
- Uses an atomic lock directory. Same unit active: refuse. Any other active unit: refuse unless `--parallel-ok` is supplied, making disjointness an explicit operator claim.
- Unique outputs: `prompt.md`, `final.md`, `stream.log`, `status.tsv`, `pid`.
- Foreground mode exits with Codex’s status unless final/expected outputs are absent, in which case it exits nonzero. Background mode returns after writing PID and paths; the child writes final status and byte counts.
- Prints byte counts for rendered prompt, final message, stream log, and every expected artifact.
- Refuses empty prompt, missing files, unsupported web capability, missing final output, absent expected artifacts, and inherited stale lock with a live PID.
- Must not run git write commands, install dependencies, mutate project config, create project-local prompt scaffolds, use `pgrep`/`pkill` name matching, or infer approval for unsandboxed execution.

Sketch:

```bash
lock="$run_root/locks/$unit.lock"
mkdir "$lock" || exit 20
trap 'rm -rf "$lock"' EXIT

if other_locks_live && [ "$parallel_ok" != 1 ]; then exit 21; fi
rendered="$run_dir/prompt.md"
cat ${rules:+"$rules"} "$prompt" > "$rendered"

codex_args=(exec -m "$model" -C "$workdir" --skip-git-repo-check -o "$final")
[ "$effort" ] && codex_args+=(-c "model_reasoning_effort=$effort")
[ "$web" = 1 ] && codex_args+=(-c tools.web_search=true --enable web_search_request)
[ "$bypass" = 1 ] && codex_args+=(--dangerously-bypass-approvals-and-sandbox)

codex "${codex_args[@]}" - < "$rendered" > "$stream" 2>&1
status=$?
[ -s "$final" ] || status=30
for p in "${expects[@]}"; do [ -e "$p" ] || status=31; done
wc -c "$rendered" "$final" "$stream" "${expects[@]}" > "$bytes"
exit "$status"
```

**3. Project-Local State**
Durable delegation state lives in the existing memory grammar: `docs/_internal/<stage-stem>/reasoning/<topic>.md` or `session.md`, created only when there is real content. The durable record holds worker partition, acceptance baseline, decisions, deviations, and outcome.

Runtime prompts, rendered prompts, logs, final messages, byte counts, and wrappers stay temporary by default. Keep a raw prompt/transcript only when it is itself the authoritative spec, evidence, or audit record.

Generic house rules belong in the skill/assets. Project-specific constraints belong in AGENTS or the stage memory record, not in a new `docs/_internal/codex/00_house_rules.md` convention.

**4. Convergence Plan**
The surviving convention is JR-MC’s wrapper shape, generalized as shipped assets: short command, file-backed prompt, stdin, capability-aware web flag expansion, unique outputs, and locks.

Fleet sweep delivers the assets. Existing projects then converge as follows:

- JR-MC: replace untracked `.devcontainer/scripts/codex_unit.sh` with calls to mounted `delegate-cli/assets/launch.sh`; keep project-specific unit notes as inputs, not as a scio convention.
- Meta-Aging: stop using long inline `codex exec`; call `launch.sh --prompt existing_prompt.md`; record durable decisions in stage memory.
- 14782-DM: stop recreating `/tmp` wrappers; use shipped `probe.sh` and `launch.sh`.
- Toolkit examples: stop prescribing `/tmp/<unit>` and direct CLI snippets as the primary path.

No mass migration of ignored historical artifacts.

**5. Enforcement Map**
Asset-enforced: capability probe, missing CLI, missing/empty prompt, stdin invocation, web capability gate, unique output paths, lock discipline, background PID/status files, byte counts, absent final/expected output.

`scio lint`-enforceable: `assets/probe.sh` and `assets/launch.sh` exist, executable, `bash -n` clean; skill text contains no exact stale `--search` or `run_in_background`; frontmatter description stays within cap.

Operator judgement: whether units are truly disjoint, whether unsandboxed execution is acceptable, whether web access is authorized, model choice by task risk, acceptance-test validity, reviewing diffs/artifacts, visual inspection, and final trust.

Not enforceable mechanically: model semantic quality, “no git write” under full bypass in all possible forms, whether a raw prompt is audit-worthy, and whether a successful final message is true.

**6. Ownership Of The Verification Stance**
The human stance follows the human, so it belongs in dev-env/global habit: delegation replaces typing, not verification.

The skill should carry the delegation-specific checklist because it follows the repo with the assets: re-run tests, inspect touched files, check submodules as separate repos, and verify outputs against disk.

The asset owns evidence production only: paths, status, bytes, locks, and optional syntactic/log checks. It cannot own belief.

**7. What To Delete**
Cut the direct command cookbook down hard. Remove stale flag snippets, `/tmp/<unit>` examples, Claude-specific backgrounding, and long incident narratives once assets encode the mechanics. Keep concise rules for prompt construction, permissions, review, and agy parsing hazards.

Net shape: skill prose should become “when to delegate, how to prepare a prompt, run probe, run launch, review result”; assets handle flag spelling.

**8. Disagreement**
No major disagreement with the brief. The one boundary I would keep explicit: `launch.sh` should not try to prove semantic safety or replace review. It should fail on mechanical absence and record evidence. The final verification remains an operator action, with dev-env owning the habit and scio owning the delegation-specific checklist.