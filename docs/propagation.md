# How a change reaches a project

This document answers one question: **you changed something in the toolkit — what
has to happen for a project to see it?**

The answer depends on how that particular artifact lives in the project. There
are three ways, and each has its own propagation rule.

---

## 1. The two hops

Every change travels the same two hops to reach a project.

```
   toolkit commit
        │
        │   HOP 1 — VERSION
        │   git moves bytes into the project's own copy of the toolkit
        ▼
   <project>/01_modules/SciAgent-toolkit/     (submodule, at a pinned commit)
        │
        │   HOP 2 — BINDING
        │   sciagent activate wires that copy into what the harness reads
        ▼
   <project>/.claude/  .agents/  AGENTS.md
```

**Hop 1 moves bytes.** A consumer holds the toolkit as a git submodule pinned to
one commit. The submodule working tree contains that commit's files. New commits
in the toolkit reach the project when the project re-pins.

Measured 2026-08-11 across the three fleet roots (`/data1/users/antonz`,
`/data2/users/JCRLab`, `/scratch/current/antonz`): 24 toolkit checkouts, every one
of them held by a gitlink (mode `160000`).

| Checkouts | Pinned at |
|---:|---|
| 22 | `5e5347e` — the baseline every consumer still runs |
| 1 | `cf19c6d` — `pansci-rna` |
| 1 | `fb6012a` — `Gama_Vivian_DRP1_bulkRNAseq`, excluded from propagation |

**Hop 2 wires bytes.** `sciagent activate` reads the pinned copy and produces the
project's binding: symlinks, managed blocks in `AGENTS.md`, hook bodies, settings.

**A change reaches a project when it has completed every hop it requires.** Which
hops a change requires follows from the artifact's class, below.

> **Development shortcut.** When `$SCIAGENT_TOOLKIT` points at a checkout you are
> editing directly, hop 1 is already complete on every save. This is why toolkit
> work feels instant and fleet delivery takes a re-pin.

---

## 2. The three artifact classes

### Class A — Linked content

The project holds a **symlink** into the toolkit. Reads follow the link and land
on whatever bytes the pinned copy currently has.

| Project path | Kind | Toolkit target |
|---|---|---|
| `.claude/skills/<name>` · `.agents/skills/<name>` | directory symlink | `skills/<name>/` |
| `.claude/agents/<name>.md` · `.agents/agents/<name>.md` | file symlink | `agents/**/<name>.md` |
| `.claude/commands/<name>.md` · `.agents/commands/<name>.md` | file symlink | `commands/**/<name>.md` |
| `02_analysis/helpers/figure-style` · `interactive-style` | directory symlink | `lib/<name>/` (analysis repos) |

Skills mount as **directory** symlinks, so every file inside a skill — `SKILL.md`,
references, scripts, assets — is linked by that one link.

The two helper-lib mounts are hyphenated, which no Python `import` can name, so
each is paired with an **underscored shim module** written next to it
(`figure_style.py`, `figure_style.R`, `interactive_style.py`). The mount is Class
A; the shim is Class B, immediately below.

**Rule: Class A propagates on hop 1 alone.** Once the pinned copy holds the new
bytes, the next read through the link sees them.

### Class B — Materialized content

The project holds **its own copy** of the bytes, written by `activate`. The copy
keeps its content until a command rewrites it.

| Project path | How it is produced |
|---|---|
| `AGENTS.md` → `SCIAGENT:ROLES` block | rendered from the role YAML and the effective catalog |
| `AGENTS.md` → `SCIAGENT:CRAFT` block | rendered from the toolkit's `craft.yaml` |
| `.claude/hooks/*.sh` | plain `cp` of `templates/project/_common/.claude/hooks/*.sh.template` |
| `.claude/statusline.sh` | plain `cp` of the statusline template |
| `02_analysis/helpers/figure_style.{py,R}` · `interactive_style.py` | plain `cp` of the analysis shim templates (analysis repos) |
| `.claude/settings.json` | `cp` when absent; missing top-level keys backfilled when present |
| `.gitignore` → `SCIAGENT:GITIGNORE` block | rendered from a fixed list |
| `.sciagent/manifest.json`, `.sciagent/*.state`, `.sciagent/hook_state/*`, `.sciagent/helper_shim_state/*` | ownership records |

**Rule: Class B propagates on hop 1 and hop 2.** The pinned copy needs the new
template, and `activate` needs to run to rewrite the project's copy.

### Class C — Toolkit code

`bin/sciagent` and `lib/sciagent/*.sh` execute from the pinned copy. A project
runs the CLI version it has pinned.

**Rule: Class C propagates on hop 1.** Behaviour changes take effect on the next
invocation after the re-pin. Run `activate` afterwards so the project's binding is
rebuilt by the new code.

Retired `.claude/output-styles/` mounts have a teardown-only compatibility path.
`deactivate` removes toolkit-targeting symlinks from that directory even when their
former source directory is gone, while preserving regular files and outside-pointing
symlinks.

---

## 3. What each change requires

| Change in the toolkit | Hop 1 · re-pin | Hop 2 · activate |
|---|:---:|:---:|
| Edit a skill body (`SKILL.md`, references, scripts) | required | — |
| Edit an agent or command body | required | — |
| Edit a helper lib under `lib/figure-style` | required | — |
| Add a skill, agent, or command to the catalog | required | required |
| Remove a skill, agent, or command | required | required |
| Rename a skill, agent, or command | required | required |
| Change a role YAML under `roles/` | required | required |
| Change `craft.yaml` | required | required |
| Change a hook or statusline template | required | required |
| Change a `02_analysis/helpers` shim template | required | required |
| Change `settings.json.template` | required | required |
| Change `lib/sciagent/*.sh` or `bin/sciagent` | required | recommended |

The pattern is simple: **content that is linked needs one hop, and content that is
copied needs two.** Membership in the catalog is itself copied information — it
lives in the set of symlinks and in the `AGENTS.md` block — so every add, remove,
and rename needs the second hop.

---

## 4. The commands that perform each hop

### `git submodule update` — hop 1 only

Moves the project's copy to a new commit. Linked content becomes current
immediately. Materialized content keeps its existing bytes.

### `sciagent update [--to <ref>]` — hop 1 then hop 2

The propagation verb. It re-pins the submodule, then re-execs the freshly pinned
`bin/sciagent` so the new library code performs the re-activation. This is the
normal way a project takes a toolkit release.

`--no-pin` performs hop 2 alone, which suits a toolkit reached by PATH or by a
development symlink.

### `sciagent activate <role> [overlay]` — hop 2 only

Rebuilds the binding from the currently pinned copy: tears down the previous
stack, recreates every symlink, rewrites the `ROLES` and `CRAFT` blocks, and
refreshes the materialized bodies.

Activation is idempotent. Running it against an unchanged toolkit produces an
unchanged project.

### `sciagent status` — read the current state

Reports the active stack, the effective catalog, block presence and hash, and
symlink health. Use it to confirm a hop landed.

### `install.sh` — planned, for the audience without git

A collaborator or a paper reader receives a release tarball and runs `install.sh`.
That single step performs both hops at once, because the tarball is the pinned
copy. The 22 fleet projects continue to use gitlinks.

---

## 5. Ownership discipline

`activate` refreshes materialized content under a rule that protects local edits.

For each managed body the toolkit answers one question: **did the toolkit write
these exact bytes?**

| Evidence | Action |
|---|---|
| The file is absent | Write it, and record the hash |
| The content matches the current template | Adopt it silently, and record the hash |
| The content matches any version the toolkit has ever shipped | Refresh it to the current template, and say so |
| The content matches an ownership record in `.sciagent/` | Refresh it to the current template |
| The content matches something else | Keep it, warn once, and cede the file permanently |

The third row carries the weight. `templates/PROVENANCE.sha1` lists every hash each
managed template has ever had. Bodies are materialized by a plain `cp` with zero
substitution, so a project's copy is byte-identical to the template version that
wrote it. Content equal to a shipped version therefore proves toolkit authorship,
which is what licenses the overwrite.

This matters because ownership records postdate most provisionings. Measured
2026-08-11: zero consumers have a `.sciagent/hook_state` directory, so a
record-based rule alone would refresh nothing in exactly the population that
needs it. Content-provenance repairs that population.

`tools/gen-template-provenance.sh` regenerates the manifest, and
`tests/test_template_provenance.sh` fails the build when a managed template's
current hash is missing from it — or when one grows a `{{PLACEHOLDER}}`, which
would silently break the byte-identity the whole argument rests on.

The rule lives in `lib/sciagent/ownership.sh` and is used by two callers:
`claude_settings.sh` for `.claude/hooks/*.sh` and `.claude/statusline.sh`, and
`symlinks.sh` for the `02_analysis/helpers` shims. It is deliberately not a
Claude-specific function — an R helper module is not a Claude artifact — and both
callers are listed in the `VERB_MODULES` closures of every verb that loads them.

---

## 6. Symlink targets and portability

`activate` writes a **relative** symlink target whenever the toolkit resolves
inside the project. A relative target survives a change of mount point, so the
same checkout binds correctly on the host and inside a container.

Older activations wrote **absolute** container paths. Measured 2026-08-11 across
all three fleet roots — 579 skill mounts in 22 activated locations:

| Mounts | Target style | Resolves on this host |
|---:|---|---|
| 347 | relative | yes |
| 232 | absolute, every one under `/workspaces/…` | no — `/workspaces` is absent on the host |

The split is clean per location: 9 locations carry relative targets throughout, and
13 carry absolute targets throughout. Those 232 mounts resolve inside the dev
container alone.

The relative group is the recently activated work — the three Meta-Aging
sub-projects and the six STING-JR projects. The absolute group holds everything
older, including two nested `pathway-explorer` checkouts and two nested toolkit
checkouts that were themselves activated.

**Re-activating a project converts its targets to relative form.** This is one of
the concrete repairs a re-pin delivers.

---

## 7. Worked examples

### You fixed a typo in `skills/anndata/SKILL.md`

Class A. Commit in the toolkit, then re-pin the project. The next read through
`.claude/skills/anndata` sees the fix. No activation step applies.

### You added `skills/velocity-analysis/`

Class A content behind a Class B membership change. Commit, re-pin, and run
`sciagent activate`. Activation creates `.claude/skills/velocity-analysis` and
lists the skill in the `AGENTS.md` block.

### You fixed a guardrail in `no_ephemeral.sh.template`

Class B. Commit, run `tools/gen-template-provenance.sh`, commit the manifest, then
re-pin and activate. Activation compares the project's hook body against the
manifest, recognises it as an unmodified older copy, and refreshes it.

### You fixed a bug in `figure_style.py.template`

Class B, and the case that motivated managing these shims at all. Until
2026-08-12 only `sciagent new project` ever wrote them, so an already-provisioned
repo held the hyphenated mounts with no importable module beside them and no way
to receive one. Now: commit, run `tools/gen-template-provenance.sh`, commit the
manifest, then re-pin and activate. Activation recognises the project's copy as an
unmodified older one and refreshes it. A project without `02_analysis/` is not an
analysis repo and is left untouched, directory included.

### You changed the locality guard in `lib/sciagent/symlinks.sh`

Class C. Re-pin the project. The next `sciagent` invocation runs the new guard.
Run `activate` so the binding is rebuilt by the new code.

### A user customized `.claude/hooks/no_ephemeral.sh`

Their content matches no shipped version, so `activate` keeps their file, warns
once, and writes a ceded marker. Every later activation leaves the file alone in
silence. `sciagent deactivate` removes the marker along with the rest of the
project state.

---

## 8. Summary

1. Bytes travel by git. `activate` wires them in.
2. Linked content needs the re-pin. Copied content needs the re-pin and the activation.
3. Catalog membership is copied information, so adds, removes, and renames need both hops.
4. `sciagent update` performs both hops in the right order.
5. Refreshing a copy is safe because content-provenance proves who wrote it.

See also: `docs/architecture.md` §5 (CLI surface), §7 (symlink topology),
§9 (state files), §10 (idempotency).
