# Decision 1 — make docs/_internal/ its own git repository

**Owner said go, 2026-08-24.** Not yet executed. Do this before decision 2, because
the rename touches the same tree and a clean split first keeps the two apart.

## Current state, measured

- `docs/_internal/` is **not** a repository: no `.git` inside it.
- The parent (scio) **tracks 104 files** under it, force-added past the bare
  `_internal/` line in `.gitignore`. Tracked paths are never ignored, which is why
  `lint --check docs-layout` correctly reports "NOT gitignored".
- 141 files total, **107 MB**, of which **16 are `.log`**: `codex-step5b.log` 34 MB,
  `codex-step3.log` 20 MB, `codex-step2.log` 12 MB, `codex-step3b.log` 8.2 MB,
  `codex-step8.log` 6.2 MB, rest smaller.
- 7 files under `handoff/dev-env/` are shell and template assets staged for the
  dev-env repo, not memory: `provision.sh`, `harness.sh`,
  `claude-settings-defaults.sh`, `statusline.sh.template`,
  `global-AGENTS.md.template`, `project-settings.json.template`,
  `user-settings.json.template`.
- `handoff/user-level-background.md` is staged for dev-env too — see
  `interaction-stance-home.md`.

## Steps

1. **Decide the logs first.** They are raw codex transcripts from the step0–step9
   refactor. Their findings are already in `*-report.md` files beside them. Ask the
   owner: delete, or move to `/data1/users/antonz/pipeline/_transcripts/`? Do not
   commit 107 MB into a fresh repository by default.
2. `git -C <scio> rm -r --cached docs/_internal` — untracks all 104 without
   deleting any file. Commit that in scio.
3. Confirm `.gitignore` already carries `_internal/`. It does, line 11.
4. `git init` in `docs/_internal/`, then an initial commit of the Markdown.
5. Move `handoff/dev-env/` and `handoff/user-level-background.md` out to the
   dev-env repo, or to a staging path outside the memory tree.
6. `bin/scio lint --check internal-memory --project-dir .` — expect the "not its
   own repository", "parent tracks", and every "is not memory" finding to clear.
7. `bin/scio lint --project-dir .` — the `docs-layout` warning should clear too,
   because the path stops being tracked and the ignore takes effect.

## Traps

- The plan documents under `docs/_internal/plans/2026-08-22-bulkirna-api-crosscheck/`
  and every `_project/` note written this session are among the tracked 104. They
  are safe: `rm --cached` leaves files on disk. Verify before committing anyway.
- `docs/_internal/handoff/session-state.md` is the pre-grammar continuity file. The
  grammar now wants `_project/session.md` plus `_project/session-history/`. Migrate
  it rather than leaving two continuity records.
- After the split, committing memory needs `git -C docs/_internal commit`. The
  `/handoff` command already says so.
