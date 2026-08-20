# Phase 04 — restore the seam in the artifact that defines it

**Repo:** scio (file staged for dev-env) · **Parallel-safe**

## Goal

`docs/_internal/handoff/dev-env/global-AGENTS.md.template` is the **user-global**
instruction layer — it travels with the human into every repo. It currently
hardcodes one project's path grammar. Remove the paths, keep the habit.

## The violation

Lines 10–12 today:

```
- Planning: plans under `docs/_internal/plans/<date-slug>/`; …
- Reproducibility: … log non-trivial decisions under `docs/_internal/reasoning/` …
- Memory: durable memory lives only in tracked files (`AGENTS.md`, `docs/_internal/`, …)
```

Those are scio's project grammar written into a file that will be loaded in
repos with no scio at all, and that can never be corrected by re-pinning. The
routes named are also the ones phase 06 retires, so adopting this as-is would
propagate a convention we have already decided against, globally and
unfixably.

## The rule that settles it

> If it must follow the repo across humans and machines → scio.
> If it must follow this human across unrelated repos → dev-env.

Corollary, from the wider audit: **dev-env may name a path for a tool. It may
not name a path for an agent.** Watching, mounting, credentials — legitimate
(`vscode-sync-settings.sh` excluding `03_results/checkpoints` from inotify is
correct and stays). Telling an agent where to write — scio's, always.

## Changes

Rewrite lines 10–12 to state habits with no paths. The habit survives; the map
is deferred to whatever the repo declares. For example, in your own words:

- record durable reasoning **before** proceeding past a non-trivial decision,
  in whatever location this repo declares
- leave a handoff that names where the work stopped before ending substantive
  work
- durable memory means tracked files; harness auto-memory is not project state

Then audit the remaining 7 files in that directory for the same violation.
`harness.sh`, `claude-settings-defaults.sh`, `provision.sh`,
`user-settings.json.template`, `project-settings.json.template`,
`statusline.sh.template`, `README.md`. Fix any that name a project path in an
agent-facing instruction; leave paths used for tooling.

## Verify

```bash
grep -nE 'docs/_internal|02_analysis|03_results|01_modules' \
  docs/_internal/handoff/dev-env/global-AGENTS.md.template
```

Must print nothing. Then the same grep across the directory, inspecting each hit
and justifying it as tool-scoped or removing it.

The bundle is **adopted nowhere** (`#33`), so this costs nothing now and cannot
be done cheaply later.

## Do not

- Do not adopt the bundle into `~/.claude` or dev-env as part of this phase.
  That is task #33 and a separate owner decision.
- Do not edit the `dev-env` repository itself. It is clean; the audit found its
  only project-path references are legitimate tool scoping.
