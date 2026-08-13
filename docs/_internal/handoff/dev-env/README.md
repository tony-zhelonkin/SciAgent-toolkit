# dev-env handoff: Claude preferences and user provisioning

Adopt these files in `/data1/users/antonz/pipeline/dev-env`. They are staged here so
the dev-env owner can choose how its rollout layer installs user preferences and
container-global context.

| Handoff file | Destination in dev-env | Previous job |
|---|---|---|
| `claude-settings-defaults.sh` | `scripts/lib/claude-settings-defaults.sh` | Project/user defaults and statusline ensure/teardown function bodies extracted from `lib/sciagent/claude_settings.sh`. |
| `statusline.sh.template` | `claude/statusline.sh` | Claude Code model/context statusline used at project and user scope. |
| `project-settings.json.template` | `claude/project/settings.json.template` | Full former project seed, including the gold defaults and the two SciAgent hook registrations. Dev-env should retain the preference keys; SciAgent still owns its project hook registrations. |
| `user-settings.json.template` | `claude/user/settings.json.template` | User-scoped gold defaults and the user-scoped statusline command. |
| `provision.sh` | `scripts/provision-ai-user.sh` | Complete former `sciagent provision --user` verb: harness selection, dry-run, settings seed, and global context stamp. |
| `harness.sh` | `scripts/lib/harness.sh` | Detection and global-context path map used by `provision.sh`. |
| `global-AGENTS.md.template` | `claude/global/AGENTS.md.template` | Shared global context body stamped by `provision.sh`. |

The extracted shell preserves the implementation as shipped. Its destination may
source or rename the helpers to match dev-env's rollout conventions.
