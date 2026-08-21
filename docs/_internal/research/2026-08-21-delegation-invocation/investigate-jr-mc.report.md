## JR-MC delegation investigation

Project: `/data2/users/JCRLab/JR-MC-Tonsill/JR-MC`

### 1. Artifact durability

“Survives recreate” means recreating the container while remounting this same host checkout. A fresh clone is materially different because Git only restores tracked content.

| Artifact | Project-relative path | Git/tmp status | Survives recreate? |
|---|---|---|---|
| Codex runner | `.devcontainer/scripts/codex_unit.sh` | Untracked, not ignored | Yes in the current bind-mounted checkout; absent from a fresh clone |
| Runtime CLI installer | `.devcontainer/scripts/setup_ai_env.sh` | Tracked | Yes; reinstalls an unpinned Codex release into ephemeral home |
| Runbook | `docs/_internal/codex/RUN.md` | Untracked, ignored | Yes in this checkout; absent from a fresh clone |
| Shared house rules | `docs/_internal/codex/00_house_rules.md` | Untracked, ignored | Yes in this checkout; absent from a fresh clone |
| Unit prompts | `docs/_internal/codex/u{1,2,5}_*.md` | Untracked, ignored | Yes in this checkout; absent from a fresh clone |
| Noisy streams | `logs/codex/*.run.log` | Untracked, ignored | Yes in this checkout; absent from a fresh clone |
| Final-message outputs | `logs/codex/*.final.md` | Absent; path would be ignored | No artifact existed during inspection |
| Claude/Agents category links | `.claude/{skills,agents,commands}`, `.agents/{skills,agents,commands}` | Untracked, ignored symlinks | Yes in this checkout; absent until relinked in a fresh clone |
| Pinned toolkit | `01_modules/SciAgent-toolkit` | Tracked gitlink at `e33b635…` | Yes when the submodule is initialized |
| Project Codex config | `.codex/` | Absent | No |
| Toolkit-prescribed `/tmp` artifacts | `/tmp/<unit>_*` | JR-MC does not use this convention | No; `/tmp` is ephemeral |

Verified: the tracked project contains no Codex invocation text. The runner is untracked; the runbook, prompts, harness links, and logs are ignored. Thus the convention survives a same-checkout container recreation through `/workspaces`, while Git cannot reconstruct it in a new checkout.

### 2. Invocation failure and measured lengths

The project’s recovery wrapper uses:

```text
cat HOUSE_RULES UNIT_PROMPT | codex exec ... -
```

The terminal-facing commands and reconstructed expanded pipelines measure:

| Unit | Operator command | Typed characters | Expanded pipeline | Prompt through stdin |
|---|---|---:|---:|---:|
| `u1_build` | `.devcontainer/scripts/codex_unit.sh u1_build` | 44 | 373 | 10,514 |
| `u5_references` | `.devcontainer/scripts/codex_unit.sh u5_references --web` | 55 | 441 | 9,755 |

The operator command is relative. The wrapper derives `/workspaces/JR-MC` and uses absolute paths internally. Its final `-` makes Codex read the concatenated prompt from stdin.

The brief establishes the original failure sequence:

1. A long inline command wrapped between `-m` and `gpt-5.6-sol`. Codex received `-m` without a value, and Bash later interpreted `gpt-5.6-sol` as a command.
2. The retry supplied nonexistent `--search`.

The wrapper eliminates the first mechanism by absorbing the 373–441-character pipeline and carrying the roughly 10,000-character prompts over stdin.

I could not measure the failed inline command itself: no copy exists in the allowed project files. The run logs record Codex’s resulting session, not the launching shell command. The exact failure lengths therefore remain undetermined.

The wrapper launches held successfully: both logs contain valid Codex headers. Unit 5 later encountered an HTTP 404 inside its task; that is an application-level failure rather than an invocation failure.

### 3. Version and capability drift

Verified versions:

- Host: `codex-cli 0.147.0`.
- Container: `OpenAI Codex v0.149.0`, recorded in both run logs.
- Container runs used `gpt-5.6-sol`, reasoning effort `high`, `/workspaces/JR-MC`, and `danger-full-access`.

The project’s tracked setup script installs whatever Codex version is current when an ephemeral container home is rebuilt. It has no version pin, making version drift across rebuilds expected.

Capability findings:

- Host 0.147.0 help supports `-m`, `-c`, `--enable`, `-C`, `--skip-git-repo-check`, `-o`, the bypass flag, and stdin through `-`.
- `--search` is absent from host help. The brief establishes that it is also absent from container 0.149.0.
- The wrapper translates `--web` into `-c tools.web_search=true --enable web_search_request`.
- Container 0.149.0 accepted that invocation but warned that `[features].web_search_request` is deprecated, web search is enabled by default, and current configuration uses top-level `web_search`.
- The runbook contradicts itself: its invocation block says `[--search]`—53 characters—while the following paragraph correctly directs the operator to `--web`.
- No project `.codex/` exists.

The current toolkit and pinned-project `delegate-cli` files have the same SHA-256. Both advertise `--search` in the Codex guidance and quick map. That claim disagrees with both inspected CLI versions. The skill’s separate long-prompt guidance—stdin plus `-` and absolute paths—is accurate.

The explicit project-only scope excluded `~/.codex/config.toml`, so I could not inspect user-level configuration. Docker API access was also denied, preventing a fresh `codex exec --help` inside the running container.

### 4. Skill reachability

The pinned toolkit contains 85 skill directories, including `delegate-cli`.

| Perspective | `.claude/skills` status | Reachable skills |
|---|---|---:|
| Host at `/data2/users/.../JR-MC` | Symlink exists but its absolute `/workspaces/JR-MC/...` target dangles | 0 through the link |
| Container at `/workspaces/JR-MC` | Target maps to the pinned toolkit’s `skills/` directory | 85 inferred, including `delegate-cli` |

All six category links are absolute:

```text
/workspaces/JR-MC/01_modules/SciAgent-toolkit/{skills,agents,commands}
```

Consequently, all six dangle on the host. Based on the verified compose mount and the logs’ `/workspaces/JR-MC` working directory, all six resolve inside the container. Container resolution is inferred from that exact mapping because live Docker inspection was unavailable.

Therefore `delegate-cli` was reachable to Claude running inside the container. The allowed files do not show whether Claude actually loaded or followed it.

### 5. What a correct, reachable `delegate-cli` would have prevented

A correct skill would have prevented both observed interface failures:

- Its stdin/file-prompt pattern would have prevented the terminal wrap from separating `-m` from `gpt-5.6-sol`.
- An accurate 0.149.0 capability description would have excluded `--search` and supplied the supported web-search configuration, preventing the unknown-flag retry.
- An accurate capability map would have removed the need to probe the binary with repeated `strings` executions.

The existing skill could have prevented the wrapping failure. It would have encouraged the second failure because it explicitly prescribes `--search`.

### 6. Surprises and limits

- Both run logs were actively growing during inspection. Unit 1 grew from 311,251 to 509,264 bytes; Unit 5 grew from 111,509 to 174,239 bytes. At `2026-08-21T00:12:54-05:00`, neither `.final.md` existed, so completion and final exit status were still undetermined.
- The recovery machinery is durable across recreation of this specific bind-mounted checkout, yet almost entirely absent from Git. Its operational persistence and source persistence differ.
- The pinned and current toolkit copies carry identical `delegate-cli` content despite different repository commits, so the stale `--search` guidance is present on both sides.
- I could not determine the exact failed inline command length, user-level Codex configuration, or whether Claude invoked the reachable skill because those records were unavailable within the permitted project path.