## Delegation invocation investigation

Legend: **Verified** = observed in the two projects or toolkit source. **Brief-established** = supplied as settled fact. **Inferred** = reconstruction from those facts.

### 1. Artifact persistence

| Artifact | Path | Tracked / ignored / tmp | Survives container recreate? |
|---|---|---|---|
| 14782 Codex bundle: `probe.sh`, `launch.sh`, prompts, outputs | `/tmp/claude-<pid>/…/scratchpad/codex/` | tmp; absent from project | **No** |
| 14782 worker outcome records | `14782-DM/docs/_internal/plans/…` | ignored; 0 tracked files under `docs/_internal/` | **Yes**, because the project is bind-mounted; absent from a fresh clone |
| 14761 prompt bundle: 17 prompts | `14761-DM/03_results/_scratch/codex_handoff/*_prompt.md` | ignored by `03_results/_scratch/*`; 0 tracked | **Yes**, because it is inside the bind mount |
| 14761 output bundle: 21 `*_last.md`, including four `*_run_last.md` | Same directory | ignored; 0 tracked | **Yes**, because it is inside the bind mount |
| 14761 wrapper | Project-wide search excluding the toolkit submodule | **Absent** | N/A |
| Pinned `delegate-cli` | Both projects: `01_modules/SciAgent-toolkit/skills/delegate-cli/` | parent tracks the submodule gitlink; skill tracked inside the submodule | **Yes** |
| Claude/Agents skill mounts | Both projects: `.claude/skills/`, `.agents/skills/` | ignored; 0 tracked files | **Yes** in the existing bind mount; absent from a fresh clone |
| Codex installer | Both projects: `.devcontainer/scripts/setup_ai_env.sh` | tracked | **Yes**; reruns after container start |
| Project-local Codex configuration | `.codex/` in either project | **Absent** | N/A |
| Container-home Codex configuration | `~/.codex/config.toml` | container home | **No** after rebuild |

Ignored files under `/workspaces` survive a container recreation even though Git cannot recreate them. `/tmp` and container home do not.

### 2. Invocation shape and length

#### 14761-DM: long command held

**Verified:**

- No wrapper exists.
- The handoff directory contains 17 separate prompt files, ranging from **3,449 to 11,170 bytes** and 49–149 lines.
- Successful final captures exist, including `10B_run_last.md` through `10E_run_last.md`.
- Absolute `/workspaces/Meta-Aging/14761-DM` links occur throughout prompts and outputs.
- No handoff artifact contains the actual `codex exec` command: `rg 'codex exec'` found **0 files**.

**Brief-established:** invocation used `codex exec … - < prompt.md`. Thus the prompt was passed through stdin, with `-` as the prompt operand. Its 3–11 KB body was outside the shell command and required no shell quoting.

A reconstruction using the toolkit’s Linux flags and the observed 10B–10E paths is **327 characters**, or **329** including the leading `! `. This is a reconstruction, because no run log preserved the literal command:

```text
codex exec -m gpt-5.6-sol -c model_reasoning_effort="high" --dangerously-bypass-approvals-and-sandbox -C /workspaces/Meta-Aging/14761-DM --skip-git-repo-check -o /workspaces/Meta-Aging/14761-DM/03_results/_scratch/codex_handoff/10B_run_last.md - < /workspaces/Meta-Aging/14761-DM/03_results/_scratch/codex_handoff/10B_prompt.md
```

It held because the prompt body was on stdin and the 327-character command remained one logical shell line. It would visually occupy several ordinary terminal rows, so it retained the same hard-wrap exposure that later broke JR-MC; that break simply did not occur here.

The current project devcontainer declares `/workspaces/14761-DM`, which reduces the same reconstruction to **294 characters**. The successful artifacts use the older umbrella path, indicating those runs occurred in an environment where `/workspaces/Meta-Aging/14761-DM` existed.

#### 14782-DM: wrappers held

**Brief-established:** Claude created `probe.sh` and `launch.sh` beside file-backed prompts in a `/tmp/claude-…/scratchpad/codex/` directory.

The wrapper absorbed the long CLI expression: the operator invoked a script, while model, effort, sandbox, working-directory, output, and stdin syntax lived inside it. That removes terminal hard wrapping from the CLI command itself.

The scripts and logs are gone, as expected for `/tmp`, so I could not verify:

- their exact command lengths;
- whether `launch.sh` used explicit `- < prompt.md` or implicit stdin;
- their exact flags;
- the exact absolute scratch path.

No project-resident fallback wrapper exists.

For scale, the toolkit’s literal Linux command skeleton is **222 characters**; its shorter operator pattern is **86 characters**. These are measured toolkit examples, not recovered 14782 invocations.

### 3. Version and capability drift

**Installed version: undetermined.** Both projects dynamically install Codex from the official installer, falling back to unversioned `npm install -g @openai/codex`; nothing pins or records a release ([14761 setup script](/data1/users/antonz/projects/Meta-Aging/14761-DM/.devcontainer/scripts/setup_ai_env.sh:38), [14782 setup script](/data1/users/antonz/projects/Meta-Aging/14782-DM/.devcontainer/scripts/setup_ai_env.sh:38)). A recreated container may therefore receive a different version.

I could not query a running project container: Docker API access was denied. No project log records `codex --version`, and no project-local `.codex/` exists. Consequently, the active `~/.codex/config.toml`, including web-search configuration, is also undetermined.

Outside the pinned toolkit, tracked project text contains:

- **0** `codex exec` invocations;
- **0** Codex-specific flag declarations.

The pinned and current `delegate-cli` bodies both assume:

- `gpt-5.6-sol`;
- `model_reasoning_effort="high"`;
- `--dangerously-bypass-approvals-and-sandbox`;
- stdin prompt files;
- `--search` for web search.

The brief establishes that `--search` is absent from Codex 0.147.0 and 0.149.0 and that web search is configuration-gated. Therefore the skill’s `--search` claim is demonstrably stale ([current skill](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/skills/delegate-cli/SKILL.md:55)). The Linux skeleton also omits the explicit `-` used by the long-prompt guidance elsewhere in the same skill.

A separate path drift exists: both current devcontainers declare `/workspaces/<dataset>` ([14761](/data1/users/antonz/projects/Meta-Aging/14761-DM/.devcontainer/devcontainer.json:5), [14782](/data1/users/antonz/projects/Meta-Aging/14782-DM/.devcontainer/devcontainer.json:5)), while recorded handoffs use `/workspaces/Meta-Aging/<dataset>`. The old paths demonstrably worked in the run environment, but they are not supplied by either project’s current compose definition.

### 4. Skill reachability

| Project | Pinned skill directories | Reachable from `.claude/skills` | Reachable from `.agents/skills` | `delegate-cli` reachable? |
|---|---:|---:|---:|---|
| 14782-DM | 86 | 42 | 42 | **No** |
| 14761-DM | 86 | 60 | 60 | **No** |

Both pin the same toolkit commit, `5e5347e46e60b925453ed3acb546befa835fad90`, and both pinned toolkits contain `delegate-cli`.

In both projects:

- `.claude/skills` and `.agents/skills` are **real directories**, not category symlinks.
- They contain legacy per-item relative links such as `../../01_modules/SciAgent-toolkit/skills/<name>`.
- Those relative links resolve on both the host and inside the project container because the pinned toolkit retains the same project-relative position.
- Therefore no observed category link is “host-only” or “container-only.”

This differs from current [link.sh](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/lib/scio/link.sh:140), which creates six category-level links using `$SCIO_TOOLKIT`. Neither project has that current layout. The brief’s absolute-link side rule applies to newly materialized category links, not these legacy per-item relative links.

### 5. What a correct, reachable `delegate-cli` would have prevented

It would specifically have prevented:

1. **The invalid `--search` launch.** Accurate capability guidance would identify web search as config/capability-gated instead of supplying a nonexistent flag. The current skill would actually reproduce this failure.
2. **Putting prompt text on the command line.** Its stdin pattern, `- < /abs/prompt.md`, keeps multi-paragraph prompts out of shell quoting and argv.
3. **The split model operand.** A wrap-safe scripted or otherwise structurally multiline invocation would keep `-m` and `gpt-5.6-sol` in the same shell command. 14782’s `launch.sh` supplied that protection; 14761’s 327-character logical line happened to remain intact.
4. **Repeated capability discovery.** A version/help preflight would establish the supported model, sandbox, reasoning, and web-search surface before the long job.
5. **Mistaking `/tmp` artifacts for durable records.** Accurate environment guidance would state that the 14782 prompt, wrapper, logs, and final capture disappear on container recreation.

Reachability alone was insufficient here: both pinned toolkits contained the skill, but neither harness exposed it. Accuracy alone was also insufficient because the skill’s `--search` statement was stale.

### 6. Surprises and limits

- The two projects pin the same 86-skill toolkit but expose materially different subsets: **42 versus 60**.
- `delegate-cli` is present in both pinned copies and absent from all four reachable skill surfaces.
- The actual links are portable relative per-item links, unlike the current toolkit’s claimed six absolute category links.
- 14761 accumulated durable-across-recreate handoffs under an ignored project path; 14782 kept its whole launch mechanism in disposable `/tmp`.
- No actual invocation string survived in either project. The **327-character** 14761 figure is explicitly reconstructed; exact 14782 command lengths cannot be determined.
- The installed Codex version and home configuration cannot be determined from project evidence. The installer is unpinned, project logs omit the version, and Docker access was unavailable.