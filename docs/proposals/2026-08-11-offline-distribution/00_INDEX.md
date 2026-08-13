# Offline-first distribution — plan index

**Date:** 2026-08-11
**Status:** approved shape, not yet implemented
**Supersedes:** "Phase 6 — package" as specified in `sciagent-rna/docs/10_sciagent-plan-of-record.md` §3.2
**Narrows:** `docs/_internal/plans/2026-07-02-provider-agnostic/` **(untracked — see §7)** — its core insight survives; two amendments in §4 below
**Owner decision, verbatim:** *"Git is the version authority, the submodule is the scientific lock, a deterministic tarball is the portable package, and Bash performs installation and harness binding without owning network transport."*

---

## 1. The decision in one paragraph

SciAgent distributes itself **offline-first, in Bash, with Git commits and tags as the version
model**. There is no npm package, no registry, no marketplace as package authority, and no installer
that reaches the network. Claude Code marketplaces and `npx skills add` remain *optional
compatibility paths* — things a stranger may use if they already live in those ecosystems — never
SciAgent's architecture. For the analysis fleet the project-local submodule stays canonical and
overrides any global installation, because that is what preserves exact-commit reproducibility.

## 2. The three layers, kept separate

The July plan had the right core insight — **neutral source, harness-specific materialization** —
and went wrong later by mixing three concerns that must stay apart: transport, project binding, and
harness configuration.

```
SciAgent source
  skills/ + craft.yaml + hooks + helpers + agents/ + commands/
                    │
                    ▼
Offline transport                     ← owns bytes, never projects
  local bare Git hub | git bundle | tar.gz + SHA256 | explicit git clone
                    │
                    ▼
Project binder                        ← owns projects, never transport
  common:  AGENTS.md + .agents/skills + ownership receipt
  claude:  CLAUDE.md + .claude/* + hooks/settings
  codex:   usually no extra skill binding (reads .agents/skills directly)
  pi / opencode / agy: only their supported native surfaces
```

Two rules follow, and they are the load-bearing ones:

1. **Packaging must not select or configure harnesses.** An installer that decides you are a Claude
   user has merged layers 2 and 3.
2. **User-global provisioning is a separate opt-in personal bootstrap**, not part of normal project
   activation. See ADR-D4 — this is also the fix for a live doc/code contradiction.

## 3. Why per-skill symlink binding remains the portable basis

Codex reads project `.agents/skills` directly and follows symlinked skill directories; Claude Code
supports project-local `.claude/skills` and follows per-skill symlinks. Both of the harnesses that
matter therefore consume exactly what `activate` already produces. **Marketplaces are unnecessary
for either.** This is the strongest argument that Phases 1–5 got the mount shape right, and it is
independently corroborated: `vercel-labs/skills` converged on the same design from scratch —
canonical copy in `.agents/skills`, per-harness symlinks into it, a lock file for ownership.

Refs:
- Codex skill locations — https://learn.chatgpt.com/docs/build-skills
- Claude skill locations — https://code.claude.com/docs/en/skills

## 4. Two amendments to the July 2026 provider-agnostic plan

`2026-07-02-provider-agnostic/01_ARCHITECTURE.md:72` still states the durable part correctly. What
changes:

| July plan | Amendment |
|---|---|
| Tier 1 user-global provisioning as a normal tier | Becomes opt-in personal bootstrap, out of the activate path (ADR-D4) |
| Packaging/provisioning may detect and seed harnesses | Packaging must not select or configure harnesses (ADR-D3) |
| Tier 2 adapters, whole-toolkit scope | Kept, **narrowed to project scope** (ADR-D5) |

**What does not exist in code today.** There is no `lib/sciagent/harness/{common,claude,codex,…}.sh`.
Project activation remains Claude-shaped. The toolkit is only *partially* harness-agnostic: neutral
source dirs, `AGENTS.md` + `.agents/skills`, detection and global-context paths for five harnesses,
but full project materialization for Claude plus `.agents` mirrors only.

## 5. Verified findings that motivated this (with evidence)

Recorded because three of them correct claims made earlier in the design conversation, and a plan
that silently drops a correction invites re-proposing the error.

**A root `plugin.json` is not an Agent Skills distribution standard.** The Agent Skills spec
(https://agentskills.io/specification) defines a `SKILL.md` directory layout and frontmatter and
*deliberately* no repository package manifest or marketplace protocol — it is transport-neutral by
design. Verified in the vendored client at `sciagent-rna/_ref/skills/src/plugin-manifest.ts`: it
reads **only** `.claude-plugin/marketplace.json` (lines 79, 125) and `.claude-plugin/plugin.json`
(lines 104, 167). It never reads a root `plugin.json`. Discovery of google-deepmind/science-skills
works because the client scans the conventional `skills/` directory
(`plugin-manifest.ts:74`, `skills.ts:266`), not because of its seven-line root manifest — that file
is metadata. **Earlier claim that a root `plugin.json` is what makes a repo discoverable: wrong.**

**The telemetry concern is concrete, and opt-out is weaker than not using the client.**
`_ref/skills/src/telemetry.ts:1-2` hardcodes `https://add-skill.vercel.sh/t` and `/audit`; it
transmits CLI version, detected agent, and install parameters (`:151-167`). Opt-out exists —
`:87` gates on `!process.env.DISABLE_TELEMETRY && !process.env.DO_NOT_TRACK` — but under a policy of
"do not depend on a corporate installer that attempts outbound traffic," an opt-out flag is a weaker
guarantee than an installer that has no network code path at all.

**`npx` is not network-neutral.** No SciAgent npm publication is needed, but `npx skills add` still
downloads and executes the installer package through npm on every invocation.

**Version pinning is possible in the distribution channel, so that is not the reason to reject it.**
`_ref/skills/src/source-parser.ts` handles `owner/repo#ref` and `tree/<ref>/<path>`; the lock entry
carries `ref` explicitly for ref-aware updates. The real reason the submodule stays canonical is
different and narrower: only it pins **the whole toolkit as one unit** — `bin/sciagent`, `lib/`,
hooks, `craft.yaml`, and the skills together — against the exact commit an analysis ran. **Earlier
claim that the plugin channel "cannot pin a version": wrong.**

**The standalone-skill problem is smaller than feared — but the count was wrong.** 83 active skills
(`skills/*/SKILL.md`). The first scan reported **8**. A
wider re-scan on 2026-08-11, reproduced independently twice, finds **15**. The original pattern set
keyed on toolkit-ish paths and missed plain `02_analysis/config|modules|scripts` references. **The
number in the first draft of this plan was too low; treat 15 as the count.**

```
# the original 8, all reproduced
anndatar-seurat-scanpy-conversion   bulk-rnaseq-gsea
decision-gate-notebook              figure-style
interactive-breakpoint-explorer     mllmcelltype-consensus-annotation
reasoning-trace                     scrna-pipeline-conventions
# 7 the first scan missed
annotate-bulk-rnaseq-data           bulk-rnaseq-activity-inference
coresh-signature-search             peak-atlas-framework
peak-atlas-multiome                 peak-atlas-unpaired
te-geneset-gsea
```

**The 15 are two different problems, and only the first is about the toolkit:**

| Flavour | Needs | Skills |
|---|---|---|
| **Toolkit coupling** — needs toolkit code reachable at runtime | `02_analysis/helpers/*` shims that source the contract libs `symlink_create_helper_lib` mounts (`symlinks.sh:547-570`) | `figure-style`, `anndatar-seurat-scanpy-conversion`, `bulk-rnaseq-gsea` (+ `mllmcelltype-consensus-annotation`, see below) |
| **Project-scaffold coupling** — needs the analysis-repo *layout*, not toolkit code | `02_analysis/config/analysis_config.yaml`, `02_analysis/stages/`, `02_analysis/notebooks/` | the remaining ~11 |

This distinction matters for ADR-D6: most of the 15 need a project-layout prerequisite in their
body prose rather than toolkit code. Single-skill export is coherent for ~68 of 83 outright and
for most of the rest once the reader satisfies the stated requirement.

`mllmcelltype-consensus-annotation`'s `.claude/skills/` reference (`bin/mllmct:5`) is arguably **not**
coupling at all: `packaged-skills.md` contract #3 *requires* the launcher to dereference its own
`BASH_SOURCE` precisely so it survives being symlinked. Flagged by the pattern, correct by design.

> **Superseded 2026-08-11 by the semantic audit — the two-flavour table above is wrong in three
> ways.** It is kept because this file deliberately preserves corrected claims. What the audit found
> by reading all 15 skills end-to-end rather than grepping them:
>
> 1. **`bulk-rnaseq-gsea` is not toolkit-coupled.** Its "toolkit" is `RNAseq-toolkit`, a *different*
>    submodule under `01_modules/`. `symlink_create_helper_lib` mounts exactly two directories —
>    `lib/figure-style` and `lib/interactive-style` (`symlinks.sh:596-604`; the earlier
>    `symlinks.sh:547-570` citation in the table is also stale) — and `RNAseq-toolkit` is neither.
>    Filing it under toolkit coupling asserts a dependency on code SciAgent does not ship.
> 2. **`anndatar-seurat-scanpy-conversion` is not coupled at all.** It is one of six audited false
>    positives, with `mllmcelltype-consensus-annotation` (already suspected just above),
>    `peak-atlas-framework`, `peak-atlas-unpaired`, `coresh-signature-search`, and
>    `scrna-pipeline-conventions`. They *mention* the layout; they do not *require* it.
> 3. **The set spans four kinds of requirement.** Toolkit code is needed by **2** skills —
>    `figure-style` and `interactive-breakpoint-explorer`; that is the ceiling by the loop body
>    cited above. `peak-atlas-multiome` needs the sibling `peak-atlas-framework` skill, four skills
>    need `RNAseq-toolkit`, `TE-RNAseq-toolkit`, or `pathway-explorer`, and the remainder rely on
>    project-layout paths. Six skills this scan never flagged are genuinely coupled:
>    `cellranger-multi-to-anndata`, `consensus-nmf-multirun`, `scrna-cxg-host`,
>    `bulk-rnaseq-pathway-explorer`, `iterative-peak-merging`, `delegate-cli`.
>
> Net: still 15 skills with reader-relevant coupling, but only 9 of them appear in the list above.
> The corrected breakdown is in `50_ADRs.md` ADR-D6. The 2026-08-12 owner ruling places those facts
> in skill-body prose and accepts the loss of the machine drift guard.
>
> *How it was determined:* each flagged skill was read in full — SKILL.md plus `scripts/`,
> `references/`, `assets/`, `checks/` — and its requirements derived from what it actually reads and
> writes at runtime, not from what its prose names. The toolkit-coupling ceiling of 2 was then
> confirmed against `symlink_create_helper_lib`'s literal `for libdir in ...` loop.

**Correction to ADR-D6's wording:** its phrase "or the toolkit hooks" has **zero** hits. A grep for
`hooks/*.sh`, `PreToolUse`, `PostToolUse` across all 83 skills returns nothing — no skill references
the enforcement hooks. The clause describes a coupling that does not exist in the corpus.

**Artifact size makes `git archive` mandatory, not merely tidy.** The working checkout is 180 MB;
tracked content is 5.6 MB uncompressed. A single ignored skill virtualenv
(`skills/mllmcelltype-consensus-annotation/.venv`) is 171 MB of that. Any tar of the working tree
ships development residue; `git archive` cannot.

**Byte-identical rebuilds are the default, and the `gzip -n` finding here was wrong.** `git archive
HEAD` is deterministic (same SHA256 across runs). So is `| gzip -n`, so is a bare `| gzip`, and so
is `--format=tar.gz` — re-measured 2026-08-11 during implementation, all three identical across
runs, then confirmed independently.

*Superseded claim, kept visible:* this entry previously read "Without `-n`, gzip embeds an mtime and
the checksum changes every run — so `--format=tar.gz` is not usable for a reproducibility claim."
That is false. **gzip embeds an MTIME only when compressing a NAMED FILE**; from a pipe there is no
name and it writes zero — precisely what `-n` forces — and git's `--format=tar.gz` shells out to
`gzip -cn` anyway. `build-release.sh` still passes `-n`, now as insurance against a differing gzip
rather than as a fix for anything observed here. Full reasoning in `10_packaging_contracts.md` §1.

This is the second "verified finding" in this file to fall on re-measurement (the coupled-skill
count was the first). Both were measurements of the right thing done the wrong way, which is the
failure mode to watch for here.

## 6. Files in this plan

| File | Contents |
|---|---|
| `00_INDEX.md` | this file — decision, layers, verified findings |
| `10_packaging_contracts.md` | `scripts/build-release.sh` and `install.sh` contracts |
| `50_ADRs.md` | ADR-D1…D7, house format (`[Arch]`/`[Eng]`/`[Taste]`) |

Sequencing lives in the plan of record: `sciagent-rna/docs/10_sciagent-plan-of-record.md` §3.

## 7. Note on where this lives, and a fragility it exposes

This plan follows the house planning format (`00_INDEX` / `10_` / `50_ADRs`, the `[Arch]`/`[Eng]`/
`[Taste]` taxonomy from `kickoff.md §2`) but **not** the house *location*. Sibling plans live under
`docs/_internal/plans/`, and `.gitignore:11` ignores `_internal/`: **zero files under
`docs/_internal/` are tracked, and `git clean -xdn` reports it as removable.** The entire July 2026
provider-agnostic plan — architecture, tiers, ADRs, blast radius, progress — exists only in this one
working tree and is one `git clean -xdf` from gone.

This plan is therefore tracked under `docs/proposals/` instead, for two reasons: a decision record
that can vanish is not durable state, and this one needs to **travel with the submodule pin** so a
consumer that re-pins receives the rationale for how the toolkit is distributed.

**Recommendation (owner's call, not taken here):** track the July plan too, either by moving it
alongside this one or by narrowing `.gitignore:11` from `_internal/` to the paths that genuinely need
ignoring. Left undone deliberately — moving another plan's files is outside this plan's scope.
