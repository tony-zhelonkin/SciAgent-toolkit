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

**Demolition note, 2026-08-13.** The release gate now runs `scio lint
--check toolkit` plus the full test suite. The activation stack's
`.sciagent/manifest.json` was deleted; catalog-link ownership is carried by
each symlink target. Git remains the version authority described by this ADR.

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

**Demolition note, 2026-08-13.** The shipped project verb is `link`. It binds
the six fixed catalog surfaces and materializes project guardrails; the local
installer still performs no project or harness mutation. The explicit harness
set proposed here was removed with the adapter premise.

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

**Demolition note, 2026-08-13.** User-global provisioning moved completely
outside this toolkit, and `provision` and `activate` were retired. The toolkit
now exposes the project-scoped verbs `link`, `craft`, and `lint`.

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

**Demolition note, 2026-08-13.** The adapter split and `bind` surface were
superseded by one convergent `link` operation. It creates one whole-tree link
per category under both `.agents/` and `.claude/`, with hook materialization
at the Claude boundary. Per-entry links and the project manifest were deleted.

---

## ADR-D6 — Coupled requirements live in skill prose; no curated public subset [Eng]

**Owner ruling, 2026-08-12:** the machine-readable declaration proposed below was rescinded.
Frontmatter stays on the catalog's ordinary authoring surface, while genuine module, version,
companion-skill, and layout requirements live in `## Prerequisites`, `## When not to use`, or
`## See also`. The coupling audit remains useful evidence; its drift guard was retired with the
declarations.

**Context.** The concern was that skills installed standalone lose the guardrails they assume — the
`no_ephemeral` PreToolUse hook, `02_analysis/stages/`, `02_analysis/helpers`. A static scan over the
83 active skills flags **15** with direct references to `02_analysis/`, helpers, or Claude paths
(list and two-flavour breakdown in `00_INDEX.md` §5). The remaining 68 show no obvious coupling.
**The count in this ADR's first draft was 8 — too low**; a wider re-scan, reproduced independently,
found 7 more that the original pattern set missed. Note also that this ADR originally said "or the
toolkit hooks": that clause has **zero** hits — no skill in the corpus references the enforcement
hooks at all.

**The owed semantic audit — done 2026-08-11, and it moved the set.** The static scan measured
*mentions*, not *requirements*. Reading all 15 skills end-to-end (SKILL.md plus `scripts/`,
`references/`, `assets/`, `checks/`) gives a different picture, and the count staying near 15 hides
that it is a **different set**:

| | count | |
|---|---|---|
| flagged by the static scan | 15 | |
| …of those, genuinely coupled | **8** | 2 need toolkit code, 6 need only the project layout |
| …of those, false positives | **6** | `anndatar-seurat-scanpy-conversion`, `mllmcelltype-consensus-annotation`, `peak-atlas-framework`, `peak-atlas-unpaired`, `coresh-signature-search`, `scrna-pipeline-conventions` — they *mention* the layout, they do not *require* it |
| …of those, mis-flagged but coupled elsewhere | **1** | `peak-atlas-multiome` is not scaffold- or toolkit-coupled; it sources `peak-atlas-framework/scripts/` — a **sibling skill** |
| coupled but never flagged | **6** | `cellranger-multi-to-anndata`, `consensus-nmf-multirun`, `scrna-cxg-host`, `bulk-rnaseq-pathway-explorer`, `iterative-peak-merging`, `delegate-cli` |
| reader-facing requirements to preserve in body prose | **15** | 8 + 1 + 6 — the same size as the flagged set, overlapping it in only 9 places |

Two earlier numbers in this ADR are therefore superseded, and both are kept visible so the error is
not re-proposed:

* **"The 8 coupled skills"** — the Decision's first-draft wording, already corrected to 15 by this
  ADR's own Context. Neither number is the audited one: the reader-facing set is 15, but **6 of
  the originally flagged 15 were false positives** and **6 skills the scan never saw are coupled**.
* **"only ~3–4 genuinely need toolkit code at runtime"** — a guess. The audited number is **2**:
  `figure-style` (`lib/figure-style`) and `interactive-breakpoint-explorer` (`lib/interactive-style`).
  Those are the only two directories `symlink_create_helper_lib` mounts
  (`symlinks.sh:596-604`, the literal `for libdir in figure-style interactive-style`), so "needs
  toolkit code" cannot exceed two by construction.

Four skills depend on `RNAseq-toolkit`, `TE-RNAseq-toolkit`, or `pathway-explorer`; one depends on
the sibling `peak-atlas-framework` skill; two use SciAgent contract libraries; the remainder depend
on analysis-repo paths. Those distinctions now appear directly in the relevant skill bodies.

**Decision.** Do not curate a public subset. State genuine coupling in the skill body. Accept the
loss of a machine drift guard rather than maintaining a second frontmatter contract.

**Rationale.** A curated subset creates a second corpus to keep in sync — the same drift tax ADR-D1
refuses for version strings. The semantic audit showed that a static scan was 40% wrong in both
directions, so prose reviewed alongside each skill is the chosen source for its requirements.

**Blocks:** nothing. **Reversible:** trivially.

**Demolition note, 2026-08-13.** `symlinks.sh` was deleted, while `link.sh`
continues to materialize the two audited helper-library contracts. Citations
above to `symlinks.sh` remain the evidence for the semantic audit as performed;
their historical paths are intentionally preserved.

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

**Demolition note, 2026-08-13.** The artifact/CLI naming decision still holds.
The locality check cited above now lives in `lib/scio/link.sh`, and the
retired activation and teardown modules no longer contribute paths to a future
CLI rename.

**Rebrand note, 2026-08-13.** Owner ruling 7 supersedes the CLI half of this
decision. The artifact and installed command are both `scio`; the executable
is `bin/scio`, internal modules live under `lib/scio/`, and runtime variables
use the `SCIO_` prefix. The old command has no compatibility shim.

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

**Demolition note, 2026-08-13.** The fleet boundary remains in force. Within
one project, Scio owns six catalog links, the `SCIO:CRAFT` block, guardrail hook
bodies and registrations, and their narrow ownership records. Project removal
is an explicit manual operation; the retired teardown verb contributes no
toolkit state.

---

## ADR-D9 — Repository and consumer vendor paths become `scio` [Arch] — **DECIDED 2026-08-13**

**Context.** The in-repo rebrand establishes `scio` as the shipped identity, while the GitHub
repository and consumer vendor paths still carry the earlier name. Owner ruling 9 supersedes doc 08
§5.3's recommendation to keep that name.

**Decision.** Rename the GitHub repository to `scio` and each consumer vendor path to
`01_modules/scio/`. The repository rename is pending and is sequenced as fleet work with the pilot
re-pin. GitHub redirects the old repository URL, so existing submodule fetches continue to work
during the migration.

**Consequences.** Each consumer migration requires a `.gitmodules` edit plus
`git mv 01_modules/SciAgent-toolkit 01_modules/scio`. `README.md` and the docs deliberately retain
`01_modules/SciAgent-toolkit/` because that is where the vendored copies live until the fleet moves.

**Blocks:** the pilot re-pin and coordinated fleet migration. **Reversible:** yes, at fleet
coordination cost.

**Owner ruling, 2026-08-21 — sequenced into the sweep.** The rename runs in the same pass as the
fleet sweep rather than before or after it. The sweep already rewrites the paths this rename touches:
four vendor-dir spellings are live (`01_modules`, `01_scripts`, `01_Modules`, `01_Scripts`) and each
resolves individually, so a separate rename pass would walk 25 copies twice. Docs written between
now and the sweep keep naming `01_modules/SciAgent-toolkit/`, which stays correct until the copies
move.

---

## ADR-D10 — Project memory mirrors analysis and is enforced from the tree [Arch]

**Context.** The toolkit's always-on text says durable memory lives in tracked files under
`docs/_internal/`, while `link` writes that path into `SCIO:GITIGNORE` and `docs-layout` warns when
the ignored directory is absent. By the instruction's own definition, the mandated location
produces no durable record in the parent repository. The measured evidence is recorded in
`docs/_internal/research/2026-08-20-internal-skeleton/`: among 24 projects, 10 have no
`docs/_internal/`; two independently made it a nested Git repository, with 126 and 210 commits and
no remote; empty `handoffs/` scaffolds coexist with real handoffs in `sessions/`; the fleet uses ten
handoff spellings; the named root `_scratch/` convention has no populated instance; and one internal
tree reached 568 MB after absorbing a virtual environment. The instruction, scaffold, and mechanism
describe different systems.

**Decision.**

1. **Memory mirrors the analysis.** Memory for `02_analysis/stages/NN_<stem>.*` lives at
   `docs/_internal/NN_<stem>/`. Memory that spans or precedes stages lives under
   `docs/_internal/_project/`. The stage number is the key agents and humans already share, so this
   mapping requires no intent metadata. There is no mirror for `02_analysis/helpers/`: a helper can
   serve several stages, and its rationale belongs with a consuming stage or under `_project/`.

   **Scope is the only structure.** A stage directory holds `session.md` and flat topic notes
   beside it — `30_grn/network-selection.md`, not `30_grn/reasoning/network-selection.md`. A
   per-document category directory classifies rather than scoping: it promises a collection, usually
   delivers one file, and is the same promise-shaped thing that made `handoffs/` ship empty in five
   projects. The fleet's competing `reasoning/` spellings are the evidence that the category became
   a vocabulary problem.

   A repository with no stages has no observable work key for lint to validate, so it uses
   `_project/` alone. Inventing a second grammar would name a vocabulary nothing checks; a software
   repository's issues and branches already partition its work.

2. **A plan is a third durable form.** The two ordinary forms are a flat `<topic>.md` and
   `session.md`; `session.md` is updated in place, which keeps one current record instead of a dated
   pile. That trade is only safe where the tree has history: an in-place update in an ignored,
   non-nested tree destroys what it replaces. Decision 6 states the topology that supplies it,
   and `lint --check internal-memory` reports a populated tree that has none — the condition is
   observed rather than asserted away. Plans are also
   durable. The strongest working example in the fleet is 14782-DM's
   `docs/_internal/plans/2026-08-14_consensus-migration/`, whose 93-line `00_STATE.md` coordinates six
   numbered phase files. This toolkit ships `templates/plan/` and uses the same pattern for this
   work.

   **Owner ruling, 2026-08-21.** A plan lives at `_project/plans/<date-slug>/`, because it spans
   stages. `session.md` is the one spelling: `00_STATE.md` was the same form under another name, so a
   plan directory holds `session.md` beside its `00_INDEX.md` and `NN_<slug>.md` phase files, and
   `session.md` is updated in place. A plan is created with its first substantive file. This keeps
   one grammar — `session.md` means "where this work stands" at every level of the tree — and it is
   what the `internal-memory` lint already checks for.

   A plan earns a *directory* by being phased: the lint requires a non-empty `00_INDEX.md`, because
   the phase map is what a directory buys over a single file. A plan that fits in one file is a
   topic note beside the stage it serves.

3. **Scratch is disposable work.** Scio names no sanctioned scratch location. Agents and tools may
   use an appropriate disposable workspace without turning it into durable memory.

4. **Directories arrive with content.** Scaffolding creates a directory with its first authentic
   file. A populated example carries the structure agents can copy; `.gitkeep` carries no usable
   contract.

5. **Tree lint is the enforcement surface.** Enforcement is the opt-in, harness-blind
   `scio lint --check internal-memory`, which reads the repository tree. Codex has a full hook
   surface — `PreToolUse`, `Stop`, and `SessionStart` in `~/.codex/hooks.json`, trusted by definition
   hash — and using it is rejected. Project enforcement would require a project `.codex/`, a
   directory present in zero of the 24 surveyed projects, and would duplicate a Claude-shaped
   mechanism. No Claude, Codex, or Git hook is added. A `pre-commit` invocation of the lint remains
   available as a later escalation if evidence shows agents skip it.

6. **Memory is its own repository, always.** `docs/_internal/` is a Git repository in its own
   right, in every project. The parent already ignores the path, so the nesting needs no additional
   wiring. The parent tracks `docs/internal-memory.md`, which names the memory location and makes
   the nested repository visible.

   Two properties follow, and both are the point. An in-place update of `session.md` is safe,
   because the repository holds what the update replaced — decision 2 depends on this. And the
   memory's **visibility is a separate decision from the code's**: the owner chooses whether the
   reasoning is published, on its own schedule, rather than inheriting the parent's answer.

   **Publication is a submodule, at the moment of publication.** When the reasoning should ship —
   at paper submission, or when a repository goes public — the parent embeds the memory as a
   submodule and the whole history travels with a durable pin. A submodule is *rejected as the
   day-to-day durability mechanism*, because durability would then depend on pin bumps and 20 of 24
   surveyed copies share one stale pin. It is the right form for a deliberate release, where the pin
   is the citation. Force-adding the tree into the parent as plain files is rejected outright: it
   gives one tree two histories, floods the human log with agent churn, and destroys the choice
   above.

   **No verb creates this.** `link` writes no repository and no pointer: a binding verb that
   initialises a Git repository as a side effect would surprise, and the two projects that adopted
   this topology did it by hand without help. What the toolkit does instead is *notice*:
   `lint --check internal-memory` reports a populated tree that is not its own repository, and
   reports plain-tracked files sitting beside one.

   **Still open: reachability.** Both nested repositories found in the survey have no remote, so
   they carry history and still die with the disk. That is a backup question rather than a topology
   question, and it is the subject of
   `docs/_internal/research/2026-08-21-memory-publication/`.

   *Owner ruling, 2026-08-23.* The nested repository is not a recommendation to be weighed against
   alternatives — it is the topology, in every project. Publication follows the parent's public /
   private knob, at the owner's discretion, via the submodule route above. An earlier draft of this
   decision ranked three options and left publication open, which is why the design had to be
   re-derived from the gitignore rule more than once.

7. **Ownership follows the seam.** Scio owns path grammar, examples, always-on router text, and lint
   predicates. dev-env owns user-global habits. The test is: *follows the repository → Scio; follows
   the human → dev-env.* dev-env may name a path for a tool, such as a watcher or credential mount;
   agent-facing paths belong to Scio.

**Consequences.** The unsupported tracked-memory promise, root `_scratch/` directive, and
`docs-layout` remediation claim become deletable, as do the pre-created empty memory categories and
their `.gitkeep` files. Migration does not bulk-add or relocate the 8,317 existing files. Existing
memory is handled opportunistically when its project next becomes active. No `docs/_internal/` is
created solely to conform: absence is the majority state and is legitimate. This ADR records the
nested-repository topology and parent pointer; it does not implement either or determine where
memory is pushed.

The review question generalises beyond memory: *for every claim the toolkit makes, what enforces it,
and can that thing actually deliver?* Four defects found in one day had the same shape: an
instruction named something its mechanism did not guarantee.

**Blocks:** nothing. **Reversible:** yes; the deletions are recoverable from history and no data is
destroyed.
