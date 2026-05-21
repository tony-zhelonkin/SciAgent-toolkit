# Meta-map

Produce the portfolio-level cartography artifact `docs/_meta/map.md` from the per-feature maps + reviews + designs that already exist.

This is the entry point of the meta-layer. Use it when ≥2 features are in flight and you want one place that names their shared touch-points, convergent concerns, and inter-feature dependencies.

## Phase 0: Parse arguments

- `$ARGUMENTS` — optional comma-separated list of feature slugs (e.g. `normalisation,umap,theme-bundles`)
  - Empty / absent → auto-detect: every subdirectory of `docs/` that contains a `map.md`, except `docs/_meta/`
  - Explicit list → use exactly those slugs

Announce the resolution before proceeding:

```
Resolved meta-map scope:
  features = [{slug}, {slug}, ...]
  output    = docs/_meta/map.md
```

If a slug was passed that has no `docs/{slug}/map.md`, warn and ask whether to drop it or run `/map {slug}` first. WAIT.

## Phase 1: Verify prerequisites

The meta-layer is only useful when ≥2 features have been mapped. Count features in scope that have `docs/{slug}/map.md`.

- **<2 features have a map.md** → reject:

  ```
  Meta-layer requires ≥2 mapped features. Currently: {N}.

  Run /map <slug> on at least one more feature, then re-run /meta-map.
  ```

  Stop.

- **≥2** → continue.

## Phase 2: Check existing state

Always check whether `docs/_meta/map.md` already exists — this phase always runs, just with different outcomes:

- **Exists** → ask the user to choose:

  ```
  docs/_meta/map.md already exists ({date}, covers features [...]).

  Choose:
    1. Use existing — skip dispatch, surface what's there
    2. Regenerate — overwrite from current per-feature artifacts
    3. Append updates — re-dispatch with instruction to update only changed sections
  ```

  WAIT for choice. On (1), skip Phase 3 and go straight to Phase 4.

- **Does not exist** → state `No existing _meta/map.md — proceeding to meta-architect dispatch.` and continue. (Do NOT describe this phase as "skipped" — the check ran; there was simply nothing to prompt about.)

## Phase 3: Dispatch meta-architect

Use the Agent tool with `subagent_type: "meta-architect"`.

Prompt template:

```
Produce docs/_meta/map.md.

Features in scope: {comma-separated slugs}

For each feature read in this order:
  1. docs/{slug}/map.md (REQUIRED)
  2. docs/{slug}/synthesis.md (preferred if present)
  3. docs/{slug}/review/*.md (read all if no synthesis)
  4. docs/{slug}/design/*.md (if present, especially 03-decisions.md)
  5. docs/{slug}/audit.md or scope.md (if present)

Also read docs/_meta/design.md and docs/_meta/plan.md if they already exist (do not contradict accepted MADRs without flagging).

Output: docs/_meta/map.md following the _meta/map.md template in your system prompt. Cover:
  - Features in scope (status table)
  - Shared touch-points (files/modules touched by >1 feature, with file:line nature)
  - Convergent concerns (the cross-feature theme; ≥2 features must surface each)
  - Inter-feature dependencies (explicit blocking relationships)
  - Open cross-feature questions (OQ-M-N)
  - Reader map

Facts and citations only. Every claim cites a path inside docs/. Do not read source code.
```

The meta-architect agent has tools `Read, Grep, Glob, Write` — Grep and Glob are restricted to `docs/` (do not let it loose on the codebase). It will produce exactly one file.

## Phase 4: Surface the meta-map

After meta-architect completes (or after the user chose "use existing" in Phase 2):

1. Read `docs/_meta/map.md`
2. Print to user:

```
Meta-map ready: docs/_meta/map.md ({N} lines)

Features in scope:    {count}
Shared touch-points:  {count}
Convergent concerns:  {count}
Inter-feature deps:   {count}
Open questions (OQ-M): {count}
```

3. Quote the Convergent concerns headlines (titles only, not full entries).
4. Quote any `OQ-M-N` titles.
5. Recommend next step:

```
Recommended next steps:
  /meta-design                  — draft inter-feature ADRs (MADRs) from the convergent concerns
  /review {slug} --as <list>    — if a specific feature still needs review before meta-design
  Stop                          — if this was a snapshot only
```

## Rules

1. **`_meta/map.md` is produced only by `meta-architect`.** The main agent never writes it directly.
2. **Meta-architect reads artifacts, not code.** If meta-architect's output cites a path outside `docs/`, ask it to revise.
3. **The directory `docs/_meta/` is created if missing.** No manual setup required.
4. **Do not auto-progress to `/meta-design`.** Even if the meta-map looks complete, stop and let the user choose.
5. **`OQ-M-N` are distinct from per-feature `OQ-N`.** If meta-architect reuses a per-feature numbering scheme, ask it to renumber.
