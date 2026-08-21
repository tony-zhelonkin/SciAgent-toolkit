## 1. Artifact durability

| Artifact | Path | Git state / location | Survives container recreate? |
|---|---|---|---|
| IG1–IG4 prompts | `integration/03_results/_scratch/codex_handoff/*_prompt.md` | Gitignored by `integration/03_results/_scratch/` | Yes. They are under the `/workspaces/Meta-Aging` bind mount. |
| Codex clean final messages | `IG1_consult_last.md`, `IG2_run_last.md`, `IG3_run_last.md`, `IG4_run_last.md` beside the prompts | Gitignored | Yes, through the bind mount. |
| Codex-authored reports | `IG1_scoring_design.md`, `IG2_last.md`, `IG3_last.md`, `IG4_last.md` beside the prompts | Gitignored | Yes, through the bind mount. |
| Codex noisy stream/run log | — | Absent from `codex_handoff/` | No artifact to preserve. |
| Meta-Aging wrapper | — | Absent, as established by the brief and consistent with the sampled scripts/text | No. |
| Root delegation documentation | [AGENTS.md](/data1/users/antonz/projects/Meta-Aging/AGENTS.md:203) | Tracked | Yes. |
| Root CLI installer | [.devcontainer/scripts/setup_ai_env.sh](/data1/users/antonz/projects/Meta-Aging/.devcontainer/scripts/setup_ai_env.sh:61) | Tracked; installs the current release without a version pin | Script: yes. Installed binary: no; container home is recreated. |
| Root pinned toolkit | `Meta-Aging/01_modules/SciAgent-toolkit`, commit `5e5347e…` | Tracked gitlink | Yes when submodules are populated. |
| Root `.claude/skills`, `.claude/agents`, `.claude/commands` | — | All absent | No mounts to preserve. |
| 14616 pinned toolkit | `14616-DM/01_modules/SciAgent-toolkit`, commit `106f59f…` | Tracked gitlink | Yes when submodules are populated. |
| 14616 six category links | `.claude/{skills,agents,commands}`, `.agents/{skills,agents,commands}` | Gitignored absolute symlinks | Link bytes: yes. Container usability: no; their targets are host paths. |
| Container Codex installation/config | normally container `$HOME/.local/bin/codex` and `$HOME/.codex/config.toml` | Ephemeral container home | No; setup reinstalls Codex on the next start. |
| Project `.codex/` configuration | Meta-Aging and 14616-DM | Absent in both | N/A. |

A significant distinction is that “gitignored” does not mean ephemeral here: the handoff tree remains on the bind-mounted workspace after a container rebuild, although Git cannot reproduce it in a fresh clone.

## 2. Invocation mechanism and measured lengths

Verified from the surviving artifacts:

- Prompts were external files ranging from 4,924 to 7,713 bytes.
- The brief records the invocation form as `codex exec … - < prompt.md`.
- Codex 0.147.0 help confirms that positional `-` reads the prompt from stdin.
- The prompt contents therefore were not embedded in the shell command.
- Outputs beside the prompts prove all four handoffs completed.
- The tasks explicitly prohibited network access, so these runs had no reason to pass the invalid `--search` flag.
- No wrapper absorbed the option/path length.

The exact executed command was not preserved in a script, tracked document, or run log. Consequently, its actual character count and whether the redirection operands were absolute or relative cannot be verified honestly.

Using the project-pinned Linux command skeleton and the surviving filenames gives these measured reconstructions:

| Unit | Prompt bytes | Absolute-path command | Relative-path command |
|---|---:|---:|---:|
| IG1 | 7,713 | 336 characters | 290 characters |
| IG2 | 4,924 | 324 characters | 278 characters |
| IG3 | 5,237 | 324 characters | 278 characters |
| IG4 | 4,979 | 324 characters | 278 characters |

These counts include `-m gpt-5.6-sol`, high reasoning, sandbox bypass, `-C /workspaces/Meta-Aging`, `--skip-git-repo-check`, `-o`, positional `-`, and input redirection. They are reconstructed bounds, not claimed transcripts.

Why Meta-Aging held:

1. The 5–8 KB prompts remained outside the command line.
2. `-` and stdin redirection formed a valid prompt-input pair.
3. The model option remained a single `-m gpt-5.6-sol` token pair.
4. The handoffs avoided `--search`.
5. The documented unsandboxed form bypassed the known container bwrap restriction.

Externalizing the prompt substantially shortens the command, but the remaining 278–336-character invocation still wraps visually in a typical terminal. Visual wrapping alone does not alter shell parsing; an inserted physical newline does. The current convention offers no wrapper-level protection if copied text acquires a real newline between `-m` and its value.

Inside the prompt text, `/workspaces/Meta-Aging` is absolute while most repository artifacts are named relative to that working directory. The shell redirection path style remains undetermined.

14616-DM contained no sampled Codex prompt/output convention, invocation log, or wrapper of its own.

## 3. Version and capability drift

Verified:

- The accessible host executable is `codex-cli 0.147.0`.
- Its `codex exec --help` supports `-m`, `-c`, `-s`, `--dangerously-bypass-approvals-and-sandbox`, `-C`, `--add-dir`, `--skip-git-repo-check`, `-o`, and positional `-`.
- It does not expose `--search`.
- Host `~/.codex/config.toml` sets:

  - `model = "gpt-5.6-sol"`
  - `model_reasoning_effort = "high"`

  It contains no `web_search`, `web_search_request`, or corresponding feature setting.
- Neither project has a `.codex/` directory.
- Both devcontainers install Codex dynamically using the official installer with an npm fallback. Neither pins a version.
- Meta-Aging uses image `scdock-r-dev:v0.5.10`; standalone 14616-DM uses `v0.5.13`. Codex is installed at container start rather than baked into either image.

The exact Codex version and config inside either project container could not be verified because Docker API access was denied and the project retains no version/run header. A recreate may install a different Codex version while the ignored prompt files persist unchanged.

The project-pinned skill drift is concrete:

| Copy | Commit / hash relation | Capability claim |
|---|---|---|
| Meta-Aging root toolkit | `5e5347e…`; skill differs from current source | Claims `--search` enables web search |
| 14616-DM toolkit | `106f59f…`; skill SHA-256 matches current toolkit source | Also claims `--search` |
| Current toolkit source | Same skill hash as 14616-DM | Also claims `--search` |

Thus both pinned skills disagree with verified 0.147.0 help and with the investigation brief’s established 0.149.0 result. The correct surface is config/feature-gated web search, including `--enable web_search_request`; `--search` is not a valid exec flag in the verified versions.

## 4. Skill reachability

| Project context | Pinned skills | Host-reachable through `.claude/skills` | Container-reachable | `delegate-cli` |
|---|---:|---:|---:|---|
| Meta-Aging root | 84 | 0 | 0 | Unreachable: `.claude/skills` is absent |
| 14616-DM | 84 | 84 | 0 inferred | Reachable on host; unreachable in container |

For Meta-Aging root, all six category links are absent, so they resolve on neither side. This is especially notable because root [AGENTS.md](/data1/users/antonz/projects/Meta-Aging/AGENTS.md:203) directs Claude to read `delegate-cli`.

For 14616-DM, all six links point to absolute paths beginning:

`/data1/users/antonz/projects/Meta-Aging/14616-DM/01_modules/SciAgent-toolkit/...`

They resolve on the host. The standalone container mounts the project at `/workspaces/14616-DM`; the umbrella container mounts it under `/workspaces/Meta-Aging/14616-DM`. Neither compose file mounts the host project at its `/data1/...` pathname. Therefore the six links dangle from either container perspective. This container conclusion is inferred from the compose mounts because a live Docker check was unavailable.

This matches [link.sh](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/lib/scio/link.sh): it loops over two harnesses and three categories and creates each link directly to `$SCIO_TOOLKIT/$category`, preserving an absolute toolkit path when `SCIO_TOOLKIT` is absolute.

## 5. What a correct, reachable `delegate-cli` would have prevented

A correct and reachable skill would have prevented the `--search` failure specifically by stating that the verified CLI versions have no such flag and that web search is enabled through configuration or `--enable web_search_request`. It would also have avoided the subsequent binary-`strings` discovery exercise.

It would have told the operator before execution that:

- long prompts belong in files and should be supplied with `- < /absolute/prompt.md`;
- the Meta-Aging container’s bwrap sandbox cannot inspect project files, requiring an explicitly authorized unsandboxed invocation;
- `gpt-5.6-sol` and high reasoning should be selected explicitly or verified from configuration;
- clean final output belongs under `-o`, with a separate stream log if command provenance is needed.

The current skill does explain stdin prompting and the bwrap restriction, but it fails the “correct” condition for web search because it prescribes `--search` twice. It also leaves a 278–336-character typed invocation and provides no wrapper, so it cannot be credited with preventing the observed physical split between `-m` and `gpt-5.6-sol`; it only removes the much larger prompt body from that line.

## 6. Surprises and limits

- The successful Meta-Aging handoffs coexist with **0 of 84** skills reachable through the root Claude harness.
- 14616-DM has all six expected links, but they are host-bound and unusable in either documented container layout.
- The ignored handoff artifacts are more durable across a container recreate than Codex itself.
- The current toolkit skill and both pinned project copies retain the invalid `--search` guidance.
- No exact invocation, run stream, container Codex version, or container Codex config survives in the authorized project paths. Exact executed command length, redirection path style, and live container link resolution therefore remain unverified.