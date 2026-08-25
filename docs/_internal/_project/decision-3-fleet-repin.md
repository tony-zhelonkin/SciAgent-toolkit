# Decision 3 — push the work out to consumers, surgically

**Owner said go, 2026-08-24**, with named targets and an explicit method: update
the pin surgically, re-run `link`/`craft` only where needed, and edit `AGENTS.md`
**preserving local context**. Subagent fan-out is authorised for monitoring and
review. Not yet executed.

This is the first authorised write outside this repository. Everything before now
was read-only under `/workspaces`, `/data2` and `/scratch`.

## Named targets, measured 2026-08-24

Toolkit HEAD at the time of measurement: **`715feb83`**. Every copy below is an
ancestor of it, so each re-pin is a clean fast-forward.

| Project | Vendor path | Pin | Behind |
|---|---|---|---|
| Meta-Aging/14616-DM | `01_modules/SciAgent-toolkit` | `106f59f6` | yes |
| Meta-Aging/14761-DM | `01_modules/SciAgent-toolkit` | `5e5347e4` | yes |
| Meta-Aging/14782-DM | `01_modules/SciAgent-toolkit` | `5e5347e4` | yes |
| Meta-Aging/neuroimmune-receptor-atlas | `01_modules/SciAgent-toolkit` | `5e5347e4` | yes |
| DC-nexus | `01_modules/SciAgent-toolkit` | `5e5347e4` | yes |
| 14839-DM-cGAS | `01_modules/SciAgent-toolkit` | `5e5347e4` | yes |
| 13036-DM_DMlab_summer_2025 | **`01_Scripts/SciAgent-toolkit`** | `5e5347e4` | yes |
| STING-JR | `01_modules/SciAgent-toolkit` | `5e5347e4` | yes |

Roots: `/data1/users/antonz/projects/Meta-Aging/…` and
`/scratch/current/antonz/projects/…`.

`Meta-Aging/integration` has no `01_modules/`. `Meta-Aging` itself vendors
`RNAseq-toolkit` but no toolkit copy — confirm before touching.

**`13036-DM_DMlab_summer_2025` uses `01_Scripts/` with a capital S.** Discovery
handles it through the stage-3 glob rather than the `01_modules/` check, so `link`
works, but any command written as `01_modules/...` will miss it. It is also the
project that makes the four-spelling problem real rather than theoretical.

## Method

Per project, in this order:

1. `git -C <project> submodule update --init <vendor-path>` then fetch and check
   out the target commit in the submodule. Commit the gitlink in the project.
2. `<vendor-path>/bin/scio link --project-dir <project>` — convergent, so it is
   safe to re-run. It sweeps legacy child links and refreshes the gitignore block.
3. `<vendor-path>/bin/scio craft --project-dir <project>` — rewrites only the
   `SCIO:CRAFT` managed block. **Text outside the markers is preserved**, which is
   what "preserving local context" relies on; a hand-edited block body needs
   `--force`, and that refusal is the signal to read it first rather than override.
4. `<vendor-path>/bin/scio lint --project-dir <project>` and record the findings.
   Expect new ones: `harness-links` now reports vendor-path skill citations, and
   `internal-memory` is opt-in so it stays quiet unless asked.

## What the consumers gain from this session

Worth stating in each project's commit body, because it is why the re-pin matters:
the pre-commit hook and the caption-sweep hook both stopped depending on the
toolkit's directory name, and the caption hook stopped calling the deleted
`validate` verb — under the old copies its fallback path could not run at all.

## Traps

- **`craft` on a hand-edited block refuses rather than clobbers.** Do not reach
  for `--force` to clear the refusal; read the local edits and fold them in.
- **A hook body the user edited is user-owned** and will not refresh. Hash-and-cede
  warns once. Expect that on projects where someone tuned a hook.
- Re-pinning changes only the gitlink. `link`/`craft` are separate hops, and
  skipping them leaves the new bytes unmounted — that is the two-hop model.
- Do the pilot on ONE project, review it fully, then fan out. `DC-nexus` is a good
  pilot: standard `01_modules/` layout and not a Meta-Aging sibling, so a mistake
  is contained.
