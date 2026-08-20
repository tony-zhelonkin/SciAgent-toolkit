# Phase 05 — seven docs point at a deleted verb

**Repo:** scbio-docker (NOT scio) · **Parallel-safe**

## Goal

`sciagent new project --type analysis <dir>` was deleted in the demolition.
Seven files still document it as the second half of the two-step setup, so every
reader and every agent following them runs a command that does not exist.

## Files

`AGENTS.md`, `README.md`, `QUICKSTART.md`, `docs/devcontainer.md`,
`docs/roadmap.md`, `docs/repo-structure.md`, `docs/ai-integration.md`.

Confirm the set before editing:

```bash
cd /data1/users/antonz/pipeline/scbio-docker
grep -rln "sciagent" --include="*.md" . | grep -v "^toolkits/"
```

## What replaces it

The real workflow is:

```bash
./init-project.sh <dir>                       # render the container (scbio-docker)
cd <dir> && ./01_modules/SciAgent-toolkit/bin/scio link   # bind the catalog
                                              ./bin/scio craft  # render CRAFT
```

Two facts to get right, because both are easy to state wrongly:

1. **There is no scaffolding verb.** The three verbs are `link`, `craft`,
   `lint`. Project scaffolding was deleted with `new`. If a doc needs to
   describe how a project tree comes to exist, describe vendoring the submodule
   and running `link` — do not invent a verb.
2. **The CLI is `scio`, not `sciagent`**, and lives at
   `<project>/01_modules/SciAgent-toolkit/bin/scio`.

Keep each file's existing voice and length. This is a correction, not a rewrite.
`docs/changelog.md` is the one place that spells out versions — do not add
version numbers elsewhere.

## Verify

```bash
grep -rn "sciagent" --include="*.md" . | grep -v "^toolkits/" | grep -v "^docs/changelog.md"
```

Must print nothing except deliberate historical references in
`docs/changelog.md` (which records history and may legitimately name the old
verb in a past entry).

## Do not

- **Do not bump the submodule pin.** This repo has live pin drift (records
  `3ab3768`, checked out at the current toolkit HEAD) and resolving it is a
  separate owner step in the fleet sequence.
- Do not touch `toolkits/SciAgent-toolkit/` — that is the scio submodule and
  other phases own it.
- Do not rename anything to `scio` beyond the CLI invocation. The repo and
  vendor-path rename is ADR-D9, an undecided owner call, and doing it here would
  make this wording wrong twice.
