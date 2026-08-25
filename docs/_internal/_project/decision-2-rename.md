# Decision 2 — rename everything to scio

**Owner said go, 2026-08-24: folders and GitHub, `gh` is authed.** Not yet executed.

**Date:** 2026-08-24 · **Scope:** ADR-D9

## Premise it started from

ADR-D9 decided the rename in August and sequenced it into the fleet sweep. The
sweep is deferred until scio is parked, so the rename looked blocked. It was not:
what actually blocked it was that the rename was a **flag day**.

## What was wrong

Six places in the mechanism hardcoded `SciAgent-toolkit`, and two shipped
executable templates did too. The day a project moved to `01_modules/scio/`:

- `scio link` would not find the project's own toolkit
- `scio lint --check freshness` would stop checking the pin
- the pre-commit hook would stop linting
- the caption sweep's fallback would stop running

All four silently. That is why the sweep had to be atomic, and why it kept
looking bigger than it was.

## What is done

- `lib/scio/common.sh::_SCIO_TOOLKIT_DIRS=(scio SciAgent-toolkit)` — one list,
  sourced first by every verb. `bin/scio` decides module order; modules do not
  source each other.
- Discovery is staged: declared in `.gitmodules`, then `01_modules/`, then one
  level down. **Evidence outranks the name**; the name only breaks ties inside a
  stage. A candidate must carry `bin/scio` or `craft.yaml` to count.
- Freshness requires a candidate with its own `.git`, because `git -C` walks
  upward and would otherwise compare the project's repository and report nothing.
- The two templates find the toolkit by globbing `01_modules/*/bin/scio`, so they
  are name-independent and cannot desynchronize from the list.
- The caption hook called `validate --check captions`. That verb was deleted with
  validate.sh; its fallback path could never have worked. It calls `lint`.
- Shipped prose, `README.md`, and a new project's `analysis_config.yaml` name
  `01_modules/scio/`.

**A project may sit at either name. The sweep no longer has to be atomic.**

## What remains, measured 2026-08-24

**`gh auth status`: logged in as `tony-zhelonkin`, ssh.** So every step below is
executable without the owner present.

There are FOUR names to move, not one. This is the part the earlier note missed:

| # | Thing | Current | Notes |
|---|---|---|---|
| 1 | GitHub repository | `tony-zhelonkin/SciAgent-toolkit` | `gh repo rename scio -R tony-zhelonkin/SciAgent-toolkit`. GitHub redirects the old URL, so nothing breaks the moment it flips. |
| 2 | Local bare hub | `/data1/users/antonz/git/SciAgent-toolkit.git` | The `hub` remote. Rename the directory, then fix the remote URL. |
| 3 | This working copy | `scbio-docker/toolkits/SciAgent-toolkit` | A submodule of scbio-docker. |
| 4 | Consumer vendor dirs | `01_modules/SciAgent-toolkit` etc. | See `decision-3-fleet-repin.md`. |

Order that avoids a broken fetch: rename GitHub first (redirect covers the gap),
then update `origin` URLs, then the hub directory and its URL, then the
directories.

### This repository's own remotes

```
hub     /data1/users/antonz/git/SciAgent-toolkit.git
origin  git@github.com:tony-zhelonkin/SciAgent-toolkit.git
```

Both carry the old name. `git remote set-url` for each after the renames.

### The parent, scbio-docker

`.gitmodules` has `[submodule "toolkits/SciAgent-toolkit"]` with `path`, `url`,
and `branch = dev`. The rename edits the section name, the path and the url. Use
`git mv toolkits/SciAgent-toolkit toolkits/scio`, which moves the gitlink and the
`.gitmodules` path together, then set the url.

**The parent is on branch `feat/bulkirna-v0.5.0`, and its recorded pin is
`3ab37688`, which is behind.** The standing instruction has been not to bump that
pin. The rename changes the submodule PATH, not the pin — keep it that way unless
the owner says otherwise, and say plainly in the commit that the pin is untouched.

`.git` gitdir pointers inside a moved submodule follow automatically; `git mv`
rewrites them. Verify with `git -C toolkits/scio status` afterwards.

## Deliberately not renamed

- `docs/extractor-import-fidelity/**` and `docs/implementation-kickoff.md` carry
  `repo_root:` frontmatter recording where that work happened. A record of the
  past does not move.
- Tests that build the legacy name are fixtures, not claims.
- `docs/_internal/**` is memory.

## Trap worth carrying

Editing a hook template obliges `tools/gen-template-provenance.sh`. Without the
rehash, every consumer holding the new bytes is read as user-authored and never
refreshed again. The manifest is append-only, so older hashes stay recognized.
