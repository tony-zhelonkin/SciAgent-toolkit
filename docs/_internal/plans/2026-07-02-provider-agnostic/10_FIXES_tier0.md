# Tier 0 — Minimal fixes (the two concrete bugs)

**Effort:** ~1 hour. **Risk:** low. **No architectural commitment** — ships independently of Tiers 1–3.
Both root causes are proven (`research/03_bug_root_causes.md`).

---

## FIX 1 — Host Claude Code: "session persistence is disabled"

### Root cause (binary-proven)

Host `~/.bashrc:163` exports:

```bash
export CLAUDE_CODE_SKIP_PROMPT_HISTORY=1
```

Added under a "Reduce local prompt/history persistence and memory accumulation" comment block
(`.bashrc:152–163`). In the compiled CLI (`~/.local/share/claude/versions/2.1.198`) the guard is:

```js
function _U(){ return cfc()==="test"&&!TEST_ENABLE_SESSION_PERSISTENCE
              || P3() || st(process.env.CLAUDE_CODE_SKIP_PROMPT_HISTORY) || LGe() }
```

`_U()` returning true is *exactly* the "session persistence is disabled" state that gates:

- *"Cannot open agents — session persistence is disabled, so this conversation cannot be backgrounded."*
- *"Cannot background — session persistence is disabled, so the forked job would have nothing to resume."*

So `CLAUDE_CODE_SKIP_PROMPT_HISTORY=1` is **not** merely "don't save up-arrow input history" in this
build — it forces the whole session-persistence-disabled state, which disables backgrounding and the
agents/teams panel. The recently-added statusline is a **red herring** (it does not touch persistence).

### Fix

Remove or comment `~/.bashrc:163`:

```bash
# export CLAUDE_CODE_SKIP_PROMPT_HISTORY=1   # <-- this disables session persistence + agents/teams
```

The neighboring privacy flags (`DISABLE_TELEMETRY`, `CLAUDE_CODE_DISABLE_NONESSENTIAL_TRAFFIC`,
`DISABLE_ERROR_REPORTING`, `CLAUDE_CODE_DISABLE_AUTO_MEMORY`, `DISABLE_AUTOUPDATER`) are fine — none
gate `_U()`. Only `CLAUDE_CODE_SKIP_PROMPT_HISTORY` does. Open a new shell (or `source ~/.bashrc`) and
`env | grep CLAUDE_CODE_SKIP` returns nothing → backgrounding works.

### Why host-only (and why containers are unaffected)

`.bashrc` is the host account's rc. Containers run as `devuser` with their own home and never source
it; no compose mounts the host home. So the error is host-only. **The tension to record:** the user
*wanted* reduced local persistence for privacy, but this specific flag is too blunt — it also kills a
feature they use. If reduced retention is still desired, use `cleanupPeriodDays` in
`~/.claude/settings.json` (bounded retention) instead of the all-or-nothing env var. See ADR-P1.

### This should not silently recur

Add a `sciagent doctor`-style check (folds into `validate`, tier 1) that warns if
`CLAUDE_CODE_SKIP_PROMPT_HISTORY` is set in the environment, since it silently disables a core feature.

---

## FIX 2 — Container Claude Code: only the default dark theme, no power-user settings

### Root cause (two independent gaps)

**Gap (a) — nothing seeds container *user-level* `~/.claude/settings.json`.** The toolkit's gold
defaults template lives at `templates/project/_common/.claude/settings.json.template` (vim,
`effortLevel xhigh`, `alwaysThinkingEnabled`, `autoMemoryEnabled false`, `teammateMode auto`,
`statusLine`, hooks). But `lib/sciagent/claude_settings.sh` only ever writes the **project** file
(`_CLAUDE_PROJECT_SETTINGS=".claude/settings.json"`, line 170) — never `${HOME}/.claude/settings.json`.

**Gap (b) — even the project defaults only land on `sciagent activate`**, whose sole callers are
`activate.sh:97–98`. No devcontainer lifecycle step runs `activate`:

- Umbrella `postCreateCommand` → `postcreate.sh` only makes `/workspaces/<dataset>` symlinks + appends
  the `si` alias. It never touches `~/.claude`.
- The four dataset devcontainers (`14616-DM`, `14782-DM`, `14761-DM`, `neuroimmune-receptor-atlas`)
  have **no `postCreateCommand` at all** — only `postStartCommand` → `poststart_sanity.sh` +
  `configure_pi_models.sh`, neither of which touches `~/.claude`.
- **No compose mounts `~/.claude`.** Container home is ephemeral (tmpfs-home block commented out), so
  each rebuild starts empty and Claude's first-run writes only `{"theme":"dark"}`-class defaults.

Precedence context: Claude merges **user < project < local < managed**. So an *activated* project's
`.claude/settings.json` *would* override the user theme — but only inside that project's cwd. A fresh
install, or any non-activated dir, sees only the bare user default → the reported symptom.

### Fix (minimal, tier-0 form — hardened properly in tier 1)

Add a user-level seeding function mirroring the existing project logic, and call it from postcreate.

1. In `lib/sciagent/claude_settings.sh`, add:

```bash
# claude_settings_ensure_user_defaults
# Seed / backfill ${CLAUDE_CONFIG_DIR:-$HOME/.claude}/settings.json from the gold template.
# Idempotent, non-clobbering: reverse jq merge, existing user values always win.
claude_settings_ensure_user_defaults() {
    local dst="${CLAUDE_CONFIG_DIR:-$HOME/.claude}/settings.json"
    local src="$SCIAGENT_TOOLKIT/templates/project/_common/.claude/settings.json.template"
    [[ -f "$src" ]] || return 0
    mkdir -p "$(dirname "$dst")"
    if [[ ! -f "$dst" ]]; then cp "$src" "$dst"; echo "wrote: $dst"; return 0; fi
    command -v jq >/dev/null 2>&1 || { echo "sciagent: jq missing; cannot seed $dst" >&2; return 0; }
    local tmp; tmp=$(mktemp)
    if jq -s '.[0] * .[1]' "$src" "$dst" > "$tmp" && [[ -s "$tmp" ]]; then
        cmp -s "$tmp" "$dst" && rm -f "$tmp" || { mv "$tmp" "$dst"; echo "updated: $dst"; }
    else rm -f "$tmp"; fi
}
```

   **Note** the user-level `statusLine.command` must point at a user-level path
   (`$HOME/.claude/statusline.sh`), not `$CLAUDE_PROJECT_DIR/...` — so seed a user-level
   `statusline.sh` too, or use a template variant whose statusline path is `~/.claude/statusline.sh`.
   The hooks in the template are `$CLAUDE_PROJECT_DIR`-scoped and are harmless at user level only if
   the referenced scripts exist; **tier-0 recommendation: seed a hooks-free user template** and let
   project-level `activate` add hooks. (Tracked as ADR-P2.)

2. Wire it in. Add a `postCreateCommand` to each of the four dataset `devcontainer.json` (they lack
   one), and extend the umbrella `postcreate.sh`, to run:

```bash
"$WORKSPACE/01_modules/SciAgent-toolkit/bin/sciagent" provision --harness claude --user   # tier 1 verb
# tier-0 stopgap until `provision` exists: source claude_settings.sh and call the function directly.
```

### Recommended tier-0 stopgap vs tier-1 proper

- **Tier-0 stopgap (today):** call `claude_settings_ensure_user_defaults` from postcreate directly.
  Fixes the symptom in every container on next rebuild.
- **Tier-1 proper:** fold user-level seeding into the new `sciagent provision` verb across all
  harnesses, and consider `CLAUDE_CONFIG_DIR` → a mounted persistent volume so login + settings +
  sessions survive rebuilds (see `20_provisioning_tier1.md` and ADR-P3). Do **not** bind-mount the
  host `~/.claude` — it recouples containers to host state and breaks the satellites' "independently
  releasable" design.

---

## Tier-0 acceptance checks

- [ ] Host: new shell, `env | grep CLAUDE_CODE_SKIP` empty; `claude` can background a conversation /
      open agents.
- [ ] Container rebuild: `~/.claude/settings.json` contains `editorMode: vim`, `effortLevel`,
      `alwaysThinkingEnabled`, etc. before any `si activate`.
- [ ] Re-running postcreate twice does not change an already-seeded `~/.claude/settings.json` that the
      user has since edited (non-clobber verified).
- [ ] `sciagent validate` (or doctor) warns when `CLAUDE_CODE_SKIP_PROMPT_HISTORY` is set.
