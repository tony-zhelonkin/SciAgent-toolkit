# ADRs — offline-first distribution

Uses the toolkit's decision taxonomy (`kickoff.md §2`): **[Arch]** hard to reverse, think slowly;
**[Eng]** reversible at refactor cost, pick the cheapest option that meets the bar; **[Taste]** pick
fast, document why.

**D-series** = decisions this plan settles. Each records **Decision** (owner-approved 2026-08-11) or
**Open** (still needs input). Nothing here is a proposal awaiting review unless marked Open.

---

## ADR-D1 — Git is the version authority [Arch]

**Context.** Nothing in the repo carries a version today: no `package.json`, no `plugin.json`, no
`VERSION`, no `.github/workflows/`. Every candidate version model was therefore available. The
observed alternatives duplicate a version string across two or three JSON manifests and need a CI
job doing `jq` surgery to keep them consistent (`marimo-pair/.github/workflows/release.yml`).

**Decision.** Git commits and tags are the version model. Tags name releases; the **full** commit SHA
is the identity that appears in release metadata and in the content-addressed install directory. No
version string is duplicated into a manifest that can disagree with it.

**Consequences.** No `npm version`, no registry, no `prepublishOnly`. The corollary is that any
manifest SciAgent ever ships for compatibility (ADR-D2) must **derive** its version from Git rather
than store an independent copy — and that derivation is a test, not a CI job, because a 75-test suite
is already the place this project enforces invariants.

**Blocks:** `build-release.sh`. **Reversible:** yes, cheaply — adding a manifest later is additive.

---

## ADR-D2 — Offline-first Bash transport; marketplaces and `npx` are optional compatibility paths, not architecture [Arch]

**Context.** Three distribution channels were surveyed against the vendored references in
`sciagent-rna/_ref/` (google-deepmind/science-skills, marimo-team/marimo-pair,
vercel-labs/skills). Findings that decided it, all verified and recorded in `00_INDEX.md` §5:
the Agent Skills spec defines **no** package manifest or marketplace protocol; a root `plugin.json`
is metadata and is not what makes a repo discoverable; the reference client hardcodes telemetry and
audit endpoints; and `npx` executes a downloaded installer on every invocation.

**Decision.** Distribution is offline-first and implemented in Bash: local bare Git hub, `git
bundle`, or deterministic `tar.gz` + SHA256. Claude Code marketplaces and `npx skills add` are
**optional compatibility paths** — supported if a stranger prefers them, never the package authority,
never a dependency of the fleet.

**Rationale beyond privacy.** An opt-out env var (`DISABLE_TELEMETRY`, `DO_NOT_TRACK`) is a weaker
guarantee than an installer with no network code path. For a lab that may put an install line into a
protocol document, "cannot reach the network" is checkable; "asks politely not to" is not.

**Blocks:** everything else in this plan. **Reversible:** yes — compatibility manifests can be added
without unwinding anything.

---

## ADR-D3 — Packaging must not select or configure harnesses [Arch]

**Context.** The July 2026 provider-agnostic plan allowed provisioning to detect and seed harnesses.
In practice that merges two layers: transport (which bytes exist on this machine) and harness
configuration (what a given tool is told about them). Once merged, an installer starts making
decisions that belong to the user.

**Decision.** `install.sh` performs no harness detection and writes no harness configuration. The
harness set is always something the user typed — surfaced as `sciagent bind --harness
portable[,claude,…]` (ADR-D5). This is an amendment to the July plan.

**Consequences.** Makes the rule checkable rather than aspirational: grep `install.sh` for any
harness name and the invariant either holds or it doesn't.

**Blocks:** `install.sh`, the `bind` verb. **Reversible:** yes, but re-merging the layers would
reintroduce exactly the confusion this separates.

---

## ADR-D4 — User-global provisioning becomes an opt-in personal bootstrap [Arch]

**Context.** `docs/architecture.md:11` lists among **non-goals**: "installing harnesses, managing
MCPs, managing API keys, multi-provider posture, anything global to the user's home directory."
`lib/sciagent/provision.sh:1-8` does user-level cross-harness provisioning for five harnesses,
including seeding global context files and settings defaults. **These two statements cannot both be
true.** It is a live doc/code contradiction, not a stylistic drift.

**Decision.** User-global provisioning is a **separate opt-in personal bootstrap**, outside the
normal project-activation path. It stays in the toolkit; it stops being something a project verb can
trigger. `architecture.md`'s non-goals list is amended to say so precisely, rather than denying a
capability that ships.

**Consequences.** Resolves the contradiction in the direction of honesty about what exists. Keeps
`activate` project-scoped, which is what makes the per-project reproducibility claim meaningful.

**Blocks:** the `architecture.md` reconciliation in the merge gates. **Reversible:** yes.

---

## ADR-D5 — Extract project binding behind `common` + harness adapters, narrowed to project scope [Arch]

**Context.** The Tier-2 adapter layer from the July plan was never built. There is no
`lib/sciagent/harness/{common,claude,codex,pi,…}.sh`; project activation is Claude-shaped
throughout. The toolkit is only partially harness-agnostic today: neutral source directories,
`AGENTS.md` + `.agents/skills`, detection and global-context paths for five harnesses — but full
project materialization for Claude plus `.agents` mirrors only.

The enabling fact is that per-skill symlink binding is already portable: **Codex reads project
`.agents/skills` directly and follows symlinked skill directories; Claude Code supports project-local
`.claude/skills` and follows per-skill symlinks.** Both harnesses that matter already consume what
`activate` produces, so no marketplace is needed for either.

**Decision.** Revive the Tier-2 adapter split, **narrowed to project scope**: a `common` binder
(`AGENTS.md` + `.agents/skills` + ownership receipt) plus thin per-harness adapters. Codex typically
needs no extra skill binding. Other harnesses get only their supported native surfaces.

**Sequencing.** After the merge and the pilot re-pin, not before. This is a refactor of the code path
that Phases 1–5 just hardened and that 75 tests now cover; doing it while the branch is unmerged
would put the hardening and its refactor in the same reviewable unit.

**Blocks:** the `bind` verb. **Reversible:** at refactor cost — it is a re-shaping of existing
behaviour, and the test suite is the safety net.

---

## ADR-D6 — Coupled skills declare compatibility; no curated public subset [Eng]

**Context.** The concern was that skills installed standalone lose the guardrails they assume — the
`no_ephemeral` PreToolUse hook, `02_analysis/stages/`, `02_analysis/helpers`. A static scan over the
83 active skills flags **15** with direct references to `02_analysis/`, helpers, or Claude paths
(list and two-flavour breakdown in `00_INDEX.md` §5). The remaining 68 show no obvious coupling.
**The count in this ADR's first draft was 8 — too low**; a wider re-scan, reproduced independently,
found 7 more that the original pattern set missed. Note also that this ADR originally said "or the
toolkit hooks": that clause has **zero** hits — no skill in the corpus references the enforcement
hooks at all.

**Decision.** Do not curate a public subset. The 8 coupled skills use the Agent Skills
specification's compatibility field to state that they require a SciAgent project scaffold; the
installer may additionally warn. A semantic audit of the 15 is still owed and is cheap. Most of them
need a **scaffold**-compatibility declaration (they want `02_analysis/config|stages|notebooks`, a
project layout) rather than a **toolkit**-compatibility one (only ~3–4 genuinely need toolkit code at
runtime) — so the declaration vocabulary needs both, not one.

**Rationale.** A curated subset creates a second corpus to keep in sync — the same drift tax ADR-D1
refuses for version strings. Declaring a requirement is self-maintaining; maintaining a hand-picked
list is not.

**Blocks:** nothing. **Reversible:** trivially.

---

## ADR-D7 — Release naming and starting version [Taste] — **DECIDED 2026-08-11**

**Context.** Two values were needed before `build-release.sh` could emit a filename, and both are
user-visible. A third consideration arrived late and settled the first: **`sciagent` is already taken**
by another project, so the name was not merely a matter of taste.

**Decision.** **`scio`**, starting at **`0.1.0`**. Artifacts: `scio-0.1.0.tar.gz` + `.sha256`.

*Why `scio`.* Latin *sciō* — "I know (how to)", the practical-knowledge sense, fourth conjugation.
Four characters, unambiguous to type, and it keeps the `sci-` stem the fleet already reads as this
tool. The owner's fuller reference is the Socratic **"scio me nihil scire"** — *I know that I know
nothing* — which fits this codebase's actual engineering posture more than a slogan usually does:
ownership records so teardown reverses only what it can prove it wrote; `rc=3` for a hand-edited
block rather than a guess; "cannot verify → leave it alone"; warn-only lints over false confidence.
Put the phrase in the README as an epigraph, not in the binary name. Rejected: `nescio` / `nescire`
("I do not know") — a disclaimer is the wrong signal for a tool people must trust with their
`.claude/` directory; and `niscire`, which is not a Latin form.

*Why `0.1.0`.* The honest signal about an intentionally-unfinished corpus (warn-only lints,
`--strict` opt-in, corpus-as-backlog) is the more useful one, and a major bump is free later.
References: DeepMind `1.1.0`, marimo `0.0.18`.

**Scope — read this before renaming anything.** This ADR decides the **release artifact name only**.
Renaming the **CLI** from `sciagent` to `scio` is a separate, much larger change and is NOT approved
here: it touches ~22 consumer checkouts, `$SCIAGENT_TOOLKIT`, the `si` alias, hook paths, and the
hardcoded `./01_modules/SciAgent-toolkit` in `_guard_toolkit_locality` (`bin/sciagent`). It needs its
own ADR and its own gated pass, keeping `sciagent` as a compatibility symlink for at least one
release. Pleasingly, `si` already being the short alias makes `si` → `scio` a widening rather than a
break.

**Blocks:** nothing further — `build-release.sh` can be written now. **Reversible:** the version
trivially; the name awkwardly once published, which is why the CLI rename is deliberately deferred
rather than bundled in.

---

## ADR-D8 — Fleet operations belong to `module-vendor`, not to this toolkit [Arch]

**Context.** While explaining the distribution design, a `sciagent fleet status` / `sciagent fleet pin`
surface was proposed to solve two stated pains: "updates require 21 manual visits" and "you cannot
see which project runs which version." **The proposal was reinvention.** `module-vendor` — already
built, already in production on this workstation — implements all of it:

| Proposed | Already shipped in `module-vendor` |
|---|---|
| `fleet status` | `status` / `detail` / `json` / `watch` / `debt` |
| `fleet pin <sha>` | `pin bump` |
| bulk propagate | `sync up` (= `sync push` → `sync align` → `pin bump`) |
| exclude frozen repos | the documented **LIVING vs FROZEN** decision policy |
| `--dry-run` | dry-run **by default**; `--go` to apply; `refs/module-vendor-undo/<branch>-<ts>` |
| discovery | worktree sweep across `MODULE_VENDOR_ROOTS` (`config.sh:5`) |

More significantly, `module-vendor`'s own mental model already states the separation this plan spent
a long conversation re-deriving: *"Two things move INDEPENDENTLY: the submodule WORKING TREE (the
files) and the superproject GITLINK (the pinned SHA)."* That is precisely the bytes-vs-pin-record
split. It was not a new insight; it was an existing, implemented one.

**Decision.** Fleet-scope operations are **out of scope for this toolkit, permanently**. The
workstation's tool boundary is three-way and already drawn — `provc --help` names it explicitly as
its layer 4:

| Tool | Scope | Owns |
|---|---|---|
| `provc` | the container substrate | base image pin, dev-env layer, AI-harness provisioning |
| `module-vendor` | the fleet of toolkit worktrees | working-tree sync, gitlink pins, hubs, mirrors, drift |
| **this toolkit** | **one project** | mounts, managed blocks, ownership records, teardown |

**Consequences.** Removes two of the five motivating pains from the distribution problem entirely —
they are solved, not unsolved. What remains is narrower and should be judged on its own: **duplicated
bytes** (5 checkouts of one identical commit measured at 584 MB under the Meta-Aging umbrella, against
5.6 MB of tracked content), and offline installability for someone outside the fleet.

**Process note, recorded deliberately.** The failure here was not the design — it was proposing a
design without first reading the `--help` of a tool the owner had already written for the same
problem. Check the existing workstation tooling before designing fleet-scope anything.

**Blocks:** nothing. **Reversible:** yes, but reversing means re-litigating a boundary the owner has
already drawn twice.
