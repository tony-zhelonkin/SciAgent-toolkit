# Tier 1 — Minimal implementation (AGENTS.md-first substrate + `provision` verb)

**Effort:** ~2–3 days. **Risk:** low–medium. **Depends on:** Tier 0 (or ships alongside it).
**Delivers:** the portable substrate + user-level, cross-provider provisioning, with no per-harness
adapter refactor yet (that's Tier 2). This is the smallest change that makes SciAgent context reach
*any* AGENTS.md-reading harness and makes fresh containers come up correctly for *every* harness.

---

## T1.1 — Make `AGENTS.md` the canonical injection target (flip the default)

**What.** Today `activate` renders CRAFT/ROLES blocks into `AGENTS.md` **and** materializes a `.claude/`
tree. Formalize `AGENTS.md` as the *primary* surface and Claude's `.claude/` + `CLAUDE.md` shim as the
*first adapter*, not the default.

**Concretely:**
- Guarantee `activate`/`new` always ensure `CLAUDE.md` = `@AGENTS.md` (import shim) — the umbrella
  already does this; make it universal and idempotent.
- Materialize `GEMINI.md` = `@AGENTS.md` **only if** legacy Gemini CLI support is requested (agy reads
  `AGENTS.md` natively, so it's usually redundant — the orphan `templates/GEMINI.md.template` should
  stay orphan unless a harness needs it). Record as ADR-P4.
- No behavior change for Claude users; this is a framing + guarantee change that unblocks Tier 2.

**Why now.** Costs ~10 lines, and every subsequent tier assumes `AGENTS.md` is the source the harnesses
read. See `01_ARCHITECTURE.md:B.2`.

---

## T1.2 — New verb: `sciagent provision`

**Contract.** `provision` seeds **user-level / global** context + settings for one or more harnesses,
idempotently and non-clobberingly. It is the verb a devcontainer runs once at create time. It does
**not** activate a role (that's per-project); it establishes the baseline every project inherits.

```
sciagent provision [--harness claude,pi,codex,agy,opencode|all]
                    [--user]           # seed ~/.<harness> global config  (default)
                    [--context]        # drop global AGENTS.md-class context file
                    [--settings]       # seed power-user defaults (claude settings.json, etc.)
                    [--dry-run]
```

**Tier-1 scope:** implement the `claude` harness fully (folds in `claude_settings_ensure_user_defaults`
from Tier 0) and a **generic AGENTS.md context drop** for the others (see T1.3). Full per-harness
settings adapters are Tier 2. So Tier-1 `provision`:

- `--harness claude --settings` → seed `${CLAUDE_CONFIG_DIR:-~/.claude}/settings.json` + user
  `statusline.sh` (hooks-free user template, ADR-P2).
- `--harness all --context` → write the global AGENTS.md-class context file to each installed
  harness's global root (T1.3).

**Idempotency.** Same mechanisms as project settings: reverse `jq` merge (existing-wins), copy-if-absent
for files, `.sciagent/provision.state` ownership tag. Re-runs are no-ops unless the template changed.

**Detection.** Extend `status.sh`'s harness probe (already detects `.pi/`) into a reusable
`harness_detect()` that reports which of the five are installed (binary on PATH and/or config dir
present: `~/.claude`, `~/.pi`, `~/.codex`, `~/.gemini/antigravity-cli`, `~/.config/opencode`).
`provision --harness all` only touches detected harnesses.

---

## T1.3 — Global AGENTS.md-class context drop (the portable win)

The single highest-leverage, lowest-cost move. Author one global context file (the user's standing
reproducible-science preferences + a pointer to per-project CRAFT) and fan it to each harness's global
root, all of which read AGENTS.md-class files natively:

| Harness | Global context path | Native? |
|---------|--------------------|---------|
| pi | `~/.pi/agent/AGENTS.md` | ✅ |
| codex | `~/.codex/AGENTS.md`-class (walks up; also `project_doc_fallback_filenames`) | ✅ |
| agy | `~/.gemini/antigravity-cli/` Global Customizations Root `AGENTS.md` | ✅ |
| opencode | `~/.config/opencode/AGENTS.md` (+ `instructions[]`) | ✅ |
| claude | `~/.claude/CLAUDE.md` (import `@~/.../AGENTS.md` or inline) | shim |

Source this global file from a new toolkit template (e.g. `templates/global/AGENTS.md.template`) so
it is single-sourced and re-propagated by `provision`. **Do not LLM-author it** (the ETH finding —
LLM-written instruction files reduce success; `00_ARCHITECTURE.md:214`). Keep it terse.

---

## T1.4 — Devcontainer wiring

- Add a `postCreateCommand` to the four dataset `devcontainer.json` that currently lack one, and extend
  the umbrella `postcreate.sh`, to run `sciagent provision --harness all` (detected-only).
- Keep it out of `postStartCommand` (create-time, not every-start) so it's a one-time seed.
- Guard with the toolkit-locality rule already in `bin/sciagent` (a project shipping its own
  `01_modules/SciAgent-toolkit` provisions from *its* copy).

**Optional (ADR-P3): `CLAUDE_CONFIG_DIR` → persistent volume.** Setting `CLAUDE_CONFIG_DIR` to a mounted
named volume makes Claude login + settings + session transcripts survive container rebuilds (and lets
`provision` target a stable path). Note the caveat: `~/.claude`/config dir must be a **real dir, not a
symlink** — Claude writes to it. This is a devcontainer/compose change, orthogonal to the toolkit; do
it if rebuild-churn is painful.

---

## T1.5 — What Tier 1 deliberately does NOT do

- No provider dispatch table refactor of `activate` (Tier 2). `activate` stays Claude-only for
  *project* materialization; `provision` is the only cross-provider surface in Tier 1.
- No codex/agy/opencode *settings* schema writing (Tier 2 adapters). Tier 1 only drops the shared
  **context file** into their global roots — which is 80% of the value for 20% of the work, because
  they all read AGENTS.md natively.
- No pi extension (Tier 3).

---

## Tier-1 acceptance checks

- [ ] `sciagent provision --harness all --context` drops the global AGENTS.md-class file into each
      *installed* harness's global root; skips uninstalled ones.
- [ ] Fresh container: `postCreateCommand` runs `provision`; Claude comes up with power-user defaults;
      pi/codex/agy/opencode (if installed) see the global context file.
- [ ] Re-running `provision` after the user edits any seeded file is a no-op on the user's keys.
- [ ] `harness_detect()` correctly reports installed vs absent for all five.
- [ ] Existing Claude-only `activate` behavior is byte-for-byte unchanged.
