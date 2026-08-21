# Bug root causes (evidence)

## Bug 1 — host "session persistence is disabled" (binary-proven)

**Cause:** `~/.bashrc:163` → `export CLAUDE_CODE_SKIP_PROMPT_HISTORY=1` (under a
"Reduce local prompt/history persistence" comment block, lines 152–163).

**Proof:** in the compiled CLI `~/.local/share/claude/versions/2.1.198`, the guard predicate:

```js
function _U(){ return cfc()==="test"&&!TEST_ENABLE_SESSION_PERSISTENCE
              || P3() || st(process.env.CLAUDE_CODE_SKIP_PROMPT_HISTORY) || LGe() }
```

`_U()===true` gates the exact user-facing strings:
- "Cannot open agents — session persistence is disabled, so this conversation cannot be backgrounded."
  (fires `tengu_left_arrow_blocked`, `reason:"persistence"`)
- "Cannot background — session persistence is disabled, so the forked job would have nothing to resume."

So in this build `CLAUDE_CODE_SKIP_PROMPT_HISTORY` forces the whole disabled state, not just input
history. `P3()` = the `--session-persistence`/no-persistence CLI flag path (not in play). `LGe()` = a
further OR term (not needed; the env var alone suffices and is what's set).

**Red herring:** the recently-added statusLine has zero connection to persistence. Diff of
`settings.json.bak`→`settings.json` shows only `model/statusLine/editorMode/teammateMode/
showThinkingSummaries/autoMemoryEnabled/effortLevel` changed — none touch persistence.

**Host-only** because `.bashrc` is host-account rc; containers run as `devuser`, never source it, and no
compose mounts the host home.

**Fix:** comment `~/.bashrc:163`. Neighbors (`DISABLE_TELEMETRY`, `CLAUDE_CODE_DISABLE_NONESSENTIAL_TRAFFIC`,
`DISABLE_ERROR_REPORTING`, `CLAUDE_CODE_DISABLE_AUTO_MEMORY`, `DISABLE_AUTOUPDATER`) can stay.

---

## Bug 2 — container Claude has only default dark theme (two gaps)

**Gap (a): nothing seeds container user-level `~/.claude/settings.json`.**
`lib/sciagent/claude_settings.sh` writes only the **project** file
(`_CLAUDE_PROJECT_SETTINGS=".claude/settings.json"`, line 170). `ensure_project_defaults` (198–233) and
`ensure_statusline` (177–190) target the project cwd, never `${HOME}/.claude`. The gold defaults
(vim, `effortLevel xhigh`, `alwaysThinkingEnabled`, `autoMemoryEnabled false`, `teammateMode auto`,
`statusLine`, hooks) live in `templates/project/_common/.claude/settings.json.template`.

**Gap (b): project defaults only land on `sciagent activate`** (sole callers `activate.sh:97–98`), and
**no devcontainer lifecycle runs activate**:
- Umbrella `postCreateCommand`→`postcreate.sh`: only `/workspaces/<dataset>` symlinks + `si` alias.
- 4 dataset devcontainers: **no `postCreateCommand`** at all; only `postStartCommand`→
  `poststart_sanity.sh`+`configure_pi_models.sh` (neither touches `~/.claude`).
- `14616-DM/.devcontainer/scripts/setup_claude_mcp.sh` mentions claude but only edits `~/.bashrc` PATH;
  writes no settings; and isn't wired into any lifecycle command.

**No compose mounts `~/.claude`** (all 5 files, both services each inspected — only workspace bind,
`SSH_AUTH_SOCK:ro`, raw-data `:ro`, atlases/refcache `:ro`, neuroimmune externals `:ro`; tmpfs-home
commented out). Container home is ephemeral → each rebuild empty → Claude first-run writes only the
theme default.

**Precedence** (user < project < local < managed): an *activated* project's `.claude/settings.json`
overrides the user theme, but only inside that cwd. Fresh install / non-activated dir → bare user
default → the symptom.

**Fix:** add `claude_settings_ensure_user_defaults()` (reverse `jq` merge, existing-wins) writing
`${CLAUDE_CONFIG_DIR:-$HOME/.claude}/settings.json`; invoke from a devcontainer `postCreateCommand`
(add one to the 4 datasets; extend umbrella). Do **not** bind-mount host `~/.claude` (recouples to
host; breaks satellite independence).
