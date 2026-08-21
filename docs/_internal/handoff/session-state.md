# Session state — cold-start handoff

Living file, overwritten rather than accumulated. Read it first, then
**re-derive the numbers** with §2 before acting on them: every figure here was
measured, and every figure here decays.

Authoritative design: `docs/architecture.md`, `docs/propagation.md`. Decisions:
`docs/proposals/2026-08-11-offline-distribution/50_ADRs.md`. This file carries
only what those cannot — where the work stopped, what waits on the owner, and
the facts that cost time to rediscover.

Last verified: **2026-08-21**, toolkit at `beb933b`. 63 tests passing.

---

## 1. The shape of the thing

- **THE ASSET** — `skills/` 85, `agents/` 21, `commands/` 7, `craft.yaml`.
  Non-regenerable: it encodes traps (the CORESH Entrez-integer trap,
  species-mismatch silent failures, iterative peak merging) that took real
  analyses to learn.
- **THE MECHANISM** — `lib/scio/` 2,619 LOC + `bin/scio` 65. Three verbs, one
  concern each: `link` = availability (convergent), `craft` = standing text
  (idempotent), `lint` = enforcement (read-only).
- **THE AUTHORITY** — Git. The submodule pin **is** the lock (ADR-D1). No
  lockfile, and there should not be one.
- **THE RECORD** — `docs/` + 9 ADRs (D1–D9; D10 planned).

Three facts that constrain every future change:

- **Six directory links, not 579 file mounts.** Ownership tracking *dissolved*
  rather than being solved. The price: a directory symlink can neither filter
  nor flatten, so **source shape must equal mount shape**. That is why
  `templates/skill/` exists, why `_attic/` is top-level, why `agents/` is flat.
  Putting non-mountable content in `skills/`, `agents/` or `commands/` breaks
  every consumer.
- **`lint.sh` owns every convention it can measure.** A claim stated twice
  drifts, so CRAFT names the standing boundary, a router skill selects the
  guide, `references/` carries qualitative doctrine, `lint` owns the predicates.
- **The defect class this project keeps producing** (see §5): an instruction
  naming something the mechanism does not guarantee.

---

## 2. Re-derive the state

```bash
cd <toolkit>                       # scbio-docker/toolkits/SciAgent-toolkit
git log --oneline -1 && git rev-list --count hub/dev..dev
bash tests/run-all.sh | tail -3
bin/scio lint --check toolkit --strict && echo clean

# Fleet census. Each directory must BE its own repo root: `git -C` on an empty
# submodule dir silently walks up to the superproject and reports ITS head,
# which reads as a diverged copy that does not exist.
TK=$PWD
find /data1/users/antonz /scratch/current/antonz /data2/users/JCRLab \
     -maxdepth 6 -type d -name SciAgent-toolkit 2>/dev/null | grep -v "/\.git/" |
while read -r p; do
  top=$(git -C "$p" rev-parse --show-toplevel 2>/dev/null)
  [[ "$(readlink -f "$top")" == "$(readlink -f "$p")" ]] || { echo "NOT-A-REPO $p"; continue; }
  s=$(git -C "$p" rev-parse --short HEAD)
  printf "%-9s behind:%-4s %s\n" "$s" "$(git -C "$TK" rev-list --count "$s..dev" 2>/dev/null || echo ?)" "$p"
done | sort

cd /data1/users/antonz/pipeline/module-vendor && ./module-vendor status
```

### What that showed on 2026-08-21

`dev` = `beb933b`, **12 ahead** of `hub/dev` = `origin/dev` = `106f59f`.
`main` = `c83dfe2` both remotes. **63** tests passing, 0 failing. `toolkit` lint
clean. Tree clean.

Unpushed, oldest first:

```
7e9e72a  Route analysis code through convention guides
6cfe500  Anchor analysis code in durable conventions
42ae8a1  Prune an unused release without reporting damage
287e8e4  Bind the CRAFT block to its declared line budget
c7e8c2e  Record where the work stopped for a cold session
f234b2c  Survey how agent memory actually lands in 24 projects
e315de2  Plan the run that makes the toolkit's claims true
1c7599f  Keep delegation prompts disposable; contract, not a path
e33b635  Bring the cold-start handoff up to the session's end state
63d80a7  Record that the sweep would break every container it touches
2a1c702  Let the assets spell the flags the skill kept getting wrong
beb933b  Point the handoff at the shipped delegation assets
```

Only `42ae8a1` (uninstall exit-3 fix), `287e8e4` (CRAFT budget check) and
`2a1c702` (delegation assets) change behaviour. The rest is documentation,
planning and evidence.

Fleet — **25 real copies, 19 tracked** (`JR-MC-Tonsill/JR-MC` appeared
2026-08-21):

| Count | Commit | Behind `dev` | What |
|---|---|---|---|
| 1 | `beb933b` | 0 | scbio-docker — the canonical dev checkout |
| 1 | `e33b635` | 2 | JR-MC — new, bound in-container |
| 1 | `106f59f` | 12 | 14616-DM — the swept pilot |
| 20 | `5e5347e` | **60** | unswept (15 fleet-managed + 5 frozen by decision) |
| 1 | `cf19c6d` | 126 | PanSci — PINNED, deliberately |
| 1 | `fb6012a` | not in canonical history | Gama_Vivian — excluded by config |

`scbio-docker` pin drift: records `3ab3768`, checked out `beb933b`.

---

## 3. The pain points this design is answering

Not abstractions — each was measured on 2026-08-20/21.

**Memory dies with the container.** `craft.yaml` tells every agent durable
memory lives in tracked files under `docs/_internal/`; `link.sh:278` writes that
path into `SCIO:GITIGNORE`; `docs-layout` warns when it is *not* ignored. Across
five sampled projects: 8,317 files, 2 tracked. Two projects independently made
`docs/_internal/` a nested git repo (126 and 210 commits) — with **no remote**,
so it dies with the disk instead of the session.

**Availability failed silently.** 14761-DM: **86 skills in its own pinned
checkout, 60 reachable, 26 invisible** — `delegate-cli` among them. The
Meta-Aging umbrella is worse: `.claude/skills` does not exist, so **86 of 86**
are unreachable. An in-container agent burned a session rediscovering a
container-sandbox finding that `delegate-cli` had recorded on 2026-07-28, and
the umbrella `AGENTS.md` grew a hand-written copy of that skill because nothing
was mounted. This class is already deleted by design — membership *is* directory
membership under `link` — so it needs **delivery, not repair**. That is the
sweep.

**Structure spreads by example, not instruction.** `handoffs/` ships containing
only `.gitkeep`; agents ignored it in **every** project and wrote handoffs to
`sessions/`. Ten spellings of "handoff" exist across the fleet. Meanwhile
`03_results/<stage>/{figures,tables}/` — the part that always holds content — is
followed nearly everywhere. An empty directory is not an example, it is an
unsupported claim.

**`docs/_internal/` became a dumping ground.** DC-nexus 568 MB, 91.5% of it a
Python `.venv`. 14782-DM 231 MB around a 196 MB model checkpoint plus 210 cache
JSONs. Retention is not a later refinement.

**Scratch escapes the project.** `craft.yaml` calls root `_scratch/` "the only
sanctioned throwaway zone"; **no project has one**. Agents write to
`/tmp/claude-<pid>/<hashed-workspace>/<uuid>/scratchpad/`. `delegate-cli`
itself prescribed `/tmp` three times — prescribed behaviour, not drift. The
skill half is **fixed** by `2a1c702`; the `craft.yaml` `_scratch/` claim is
phase 01's to delete.

---

## 4. Blocked on the owner — six decisions

1. **Review `106f59f..dev`** (12 commits; 3 change behaviour). The owner intends to tweak wording
   during review. Use `nvim -c 'DiffviewOpen origin/dev'` — **not** `A..dev`: a
   commit-to-commit range makes both panes read-only because neither side is a
   file on disk. One revision diffs against the working tree, so the right pane
   is editable. From any diffview entry, `gf` opens the real file
   (`goto_file_edit`), `<C-w><C-f>` in a split, `<C-w>gf` in a tab — verified
   against the local install. `devenv/autoreload.lua` calls `DiffviewRefresh`
   automatically after saves.
2. **Plan as a third durable form** — `06_adr-d10.md` marks this OPEN. The
   skeleton says two forms (`reasoning/` + `session.md`), which contradicts the
   fleet's best-working record (14782-DM's
   `plans/2026-08-14_consensus-migration/00_STATE.md`, 93 lines, plus six phase
   files) *and* this toolkit's own shipped `templates/plan/`. Proposed:
   `_project/plans/<date-slug>/`, `00_STATE.md` **is** that plan's `session.md`
   so pick one spelling, created with content and never scaffolded.
3. **ADR-D9 rename** — repo → `scio`, vendor path → `01_modules/scio/`. Four
   vendor-dir spellings are live and each resolves individually: `01_modules`,
   `01_scripts`, `01_Modules`, `01_Scripts`. Decide before the sweep; it makes
   phase 05's wording wrong twice otherwise.
4. **Task #36** — CRAFT token-weight budget. The cap binds shape (17 of 25
   lines) but not cost (4,658 chars; longest bullet 693). A per-bullet cap fails
   today.
5. **The `seccomp=unconfined` posture** (`07_delegation-seam.md`).
   scbio-docker's compose template already ships the line, commented, on both
   services. Leaving it off means codex delegation runs unsandboxed and safety
   rests on prompt prohibitions plus a post-hoc grep — the exact defect shape
   §5 describes. Recommendation: enable it. It is a real loosening on a shared
   machine, so it is the owner's call.

Also carried: commit 14761-DM (was blocked on a live `codex exec` writing
`02_analysis/helpers/grn_viz/`; re-check).

---

## 4b. NEW 2026-08-21 — #37 blocks the sweep, decide it first

`_link_category` (`link.sh:140-172`) writes the six category links **absolute**:
`src="$SCIO_TOOLKIT/$category"`, `ln -s "$src" "$dst"`, no relativization. The
helper-shim path (`link.sh:230-231`) *does* relativize. So relativization
survived the demolition for `02_analysis/helpers/` and was lost for
skills/agents/commands in `3ab3768`. Verified live:

```
.claude/skills                   -> /data1/.../SciAgent-toolkit/skills        absolute
02_analysis/helpers/figure-style -> ../../../../../data1/.../figure-style      relative
```

`CHANGELOG.md:59` still claims all three namespaces and both mirrors are
relative. False since `3ab3768` — a **fifth** instance of §5's defect class,
this time in the changelog.

Observed in the field: `/data2/users/JCRLab/JR-MC-Tonsill/JR-MC` (new
2026-08-21, pinned at `e33b635` = current dev) has all six links pointing at
`/workspaces/JR-MC/01_modules/...` — they resolve in the container and dangle on
the host.

**Why it blocks #30:** `link` from the host rewrites all six to host paths and
breaks the container; in-container it breaks host tooling. Last writer wins. The
sweep runs `link` in 15 repos — from the host that silently breaks every
container-bound project. Fix, or fix the invocation location, before sweeping.

Also: the helper path's absolute fallback is dead — `realpath --relative-to`
always succeeds, so an external toolkit gets a `../../../../../` chain instead
of the intended absolute. Same change should fix it.

**The fleet also grew while blocked.** 19 tracked SciAgent copies now, not 18 —
`JR-MC-Tonsill/JR-MC` appeared overnight. Each day the sweep waits adds a copy
that needs it. Unswept copies are now **57** behind, not 56.

## 5. The organizing finding

Every defect surfaced on 2026-08-20 was the same defect: **an instruction naming
something the mechanism does not guarantee.**

| The toolkit says | The mechanism does |
|---|---|
| durable memory lives in tracked files, in `docs/_internal/` | `link` gitignores that path (`link.sh:278`) |
| `_scratch/` is the only sanctioned throwaway zone | no project has one |
| "`docs/_internal/` missing — run: `scio link`" | `link` never creates it (`lint.sh:604`) |
| AGENTS.md: read the skill | activation mounted 60 of 86, silently |

**The review question that falls out, worth applying to anything this toolkit
asserts: what enforces this claim, and can that thing actually deliver?**

---

## 6. The next implementation run — planned, ready to fan out

`docs/_internal/plans/2026-08-20-memory-and-seams/` — `00_INDEX.md` plus seven
phase briefs, each one bounded implementer. **Owner approved the design.**

Read `00_INDEX.md` first. §4 is a conflict map: **phases 01 and 02 both edit
`lib/scio/lint.sh` and must serialize** — do not put them in parallel worktrees.
§6 lists what is deliberately out of scope.

| # | Ships | Repo |
|---|---|---|
| 01 | delete the three unbacked claims (lint predicate, `_scratch/` claim, CRAFT routes, hook text) | scio |
| 02 | the opt-in `internal-memory` lint check that replaces them | scio |
| 03 | stop shipping `.gitkeep` and empty category dirs | scio |
| 04 | strip scio's path grammar from the user-global dev-env template | scio |
| 05 | 7 docs citing the deleted `sciagent new project` verb | scbio-docker |
| 06 | ADR-D10 — the record for all of it | scio |
| 07 | delegation seam: retention contract, cite scbio-docker for the bwrap cause | scio |
| 08 | **DONE `2a1c702`** — `probe.sh`/`launch.sh` assets; SKILL.md 416 → 150 | scio |

Fan-out for what remains: `{01→02} ∥ 03 ∥ 04 ∥ 05 ∥ 06 ∥ 07`.

**Phase 08 shipped 2026-08-21** (fanned out to two codex workers on disjoint
files, reviewed by re-running every gate rather than trusting their reports).
`delegate-cli` prescribed `--search`, which exists in **neither** codex 0.147.0
nor 0.149.0 — the skill *caused* a field failure. Its backgrounding advice was
a Claude Code tool parameter, absent when the owner runs `!` himself. It said
nothing about concurrency, and two overlapping runs of one unit wrote the same
paths for ~2.5 min until the orchestrator stopped trusting the worker's report.

Now: `skills/delegate-cli/assets/probe.sh` reads `codex exec --help` and emits
capability as `KEY=VALUE` (`HAS_SEARCH_FLAG`, `WEB_MODE`), proving tool use with
a nonce file; `launch.sh` keeps the typed command short (the wrap fix), feeds the
prompt via stdin, locks per unit so a duplicate refuses, requires
`--parallel-ok` to run beside a *different* unit, backgrounds on `--bg`, and
fails when an `--expect` artifact is absent. Verified live: `HAS_SEARCH_FLAG=0`,
`WEB_MODE=config_enable`, and a real end-to-end codex run captured.

**The generalisable lesson: a skill must not assert the flags of a CLI that
version-drifts underneath it.** Two codex versions are live in the fleet
(0.147.0 host, 0.149.0 container, installed by an unpinned `setup_ai_env.sh`).
Prose dated "verified 2026-07-28" already failed once. Probe, do not assert.

Field evidence: `docs/_internal/research/2026-08-21-delegation-invocation/` —
three per-project investigations plus a gpt-5.5 design consult. It also records
that **four** conventions for delegation artifacts existed (JR-MC
`docs/_internal/codex/`, Meta-Aging `_scratch/codex_handoff/`, 14782-DM
`/tmp/…/scratchpad/codex/`, and the skill's own `/tmp/<unit>`), that Meta-Aging
worked on **operator habit rather than mechanism**, and that JR-MC's good
wrapper was itself untracked inside a gitignored tree — the best answer in the
fleet was the least durable.

Evidence: `docs/_internal/research/2026-08-20-internal-skeleton/` — four field
inspections across 24 projects, a digest (`FINDINGS_field.md`), and three
consults (layout, enforcement, delegation). **Cite it; do not re-derive it.**

### Decided, so do not reopen

- **No new hook** — not Claude, not Codex, not git `pre-commit`. Surface
  ranking: populated filesystem > `scio lint` > git hook > harness hooks. Codex
  *does* have a full hook surface (`PreToolUse`/`Stop`/`SessionStart`,
  `~/.codex/hooks.json`, trust by definition hash) and it was **rejected**: it
  would need a project `.codex/` where none exists in any of 24 projects.
- **Never scaffold empty.** A directory arrives with its first real file.
- **Do not name a location nothing creates.** That is why the `_scratch/` claim
  is being deleted, not relocated.
- **No hand-maintained intent metadata.** The owner rejected config-knob
  anchors (`decisions.<stage>`, `status: APPROVED`) from lived experience: they
  produced churn and variables nobody could justify. Intent lives with the
  stage and the figure it produced.
- **No mirror of `02_analysis/helpers/`** — a helper serves several stages, so
  its rationale files with the consuming stage.
- **The seam.** Follows the repo → scio. Follows the human → dev-env.
  Corollary: **dev-env may name a path for a tool, never for an agent.**
- **Delegation prompts stay disposable.** They are an execution projection of a
  plan, not the plan; the durable record already exists in-tree. What was
  missing is a retention contract, not a path.

---

## 7. Traps — expensive to rediscover

**`module-vendor`'s state column says nothing about currency.** Of the 15
fleet-managed unswept copies, 5 report `BEHIND v48` and 10 read `OK clean`, yet
all 20 sit at the same commit — the 10 track an `origin/dev` that is an
*ancestor* of their own HEAD. Verify by SHA.

**`v48` and `56` are both right.** `v48` is measured against the copy's upstream
(`106f59f`); `56` is the true distance to `dev`. The gap is the unpushed
commits.

**Six copies are absent from `module-vendor status` on purpose.** `config.sh`
`MODULE_VENDOR_EXCLUDE` names five **frozen by owner decision 2026-08-18**
(12868-EH, Dakota_Black_NPC_DRP1, 13403-YD_Christina, JBader_scHFD, 13036-DM)
plus Gama_Vivian. Exclusion, not branch-pinning, is what holds them, because
`sync align` advances any copy it can see on a shared branch. **Do not "fix"
their absence.**

**`AdaW-eWAT-WL-bulkRNAseq` is benign** — a gitlink at pin `108da7a` with an
empty directory. Uninitialized submodule, nothing bound, nothing dangling.

**Never run `link` between `align` and `rename` in the sweep.** Align dangles 8
agent and 2 command mounts per repo.

**`docs/_internal/` is gitignored** (`.gitignore:11`). Everything tracked there
was force-added, so this file and its siblings need `git add -f`.

**A `craft.yaml` edit changes every consumer's CRAFT hash.** Verify the rendered
body against `HEAD:craft.yaml` before committing unless the change is *intended*
to alter it. `body:` must stay the last top-level key — `block.sh` reads it as a
trailing block scalar. `max_lines: 25` is enforced by `lint --check toolkit`.

**Release tests build throwaway fixture repos, never this one** — the builder
runs the suite, so a test building this repo would re-enter it. Hence no
`--skip-checks` flag by design.

**`lint` has 11 selectable names; `all` runs 10.** `toolkit` is opt-in and
validates this repo's own assets, so it must never gain a finding that fires in
a consumer. Phase 02 makes it 12/10 on the same principle.

**Bash only on the core path** — no `yq`, `jq`, `python3` (`AGENTS.md` rule 1).

---

## 8. Owner-only, no deadline

- `docs/guidelines/` remainder — `data_processing.md` 281, `gsea_analysis.md`
  342, `master_tables.md` 380 lines await reconciliation against the assay
  skills. Scientific, not structural.
- Decision-gates lint (#34): `notebooks/` is invisible to lint. Deferred.
- **ADR-D5 is closed, not pending.** Earlier revisions of this file said
  "decided, unimplemented" — wrong. The ADR carries its own demolition note of
  2026-08-13: the adapter split and `bind` surface were superseded by the one
  convergent `link`. No `lib/scio/harness/` layer exists to build. Task #6
  closed on that basis.
- `v-last-orchestrated` (= `5e5347e`, where the unswept copies sit) is hub-only.
  Reachable as an ancestor; the named anchor is not on origin.
- Stale `.sciagent/manifest.json` (mode 0600) in the 14616-DM and
  13403-YD_Christina toolkit checkouts — what makes 14616-DM read
  `UNTRACKED ?1`. Untouched: `rm` inside a consumer is the owner's call.

## 9. After the owner unblocks

1. Merge `dev --no-ff` → `main`; push hub **and** origin.
2. **The sweep, 15 repos** (#30) — align → rename → link → craft → commit →
   bump; repoint stale `origin/dev` upstreams to `hub/dev` in the same pass
   (#31). This is what actually delivers the 26 unreachable skills.
3. Clear the `scbio-docker` pin drift off `3ab3768`.
4. dev-env adoption of the 8-file bundle (#33) — the habit layer, live nowhere.
   dev-env clean at `0b5a4b0`.
5. Consumer `scripts/` → `stages/` (#8). Last, after the rename.

## 10. Distribution channels

| Channel | State |
|---|---|
| Vendored submodule (canonical) | shipped, proven end-to-end |
| Global install — `install.sh` + `scripts/build-release.sh` | shipped, exercised end-to-end |
| npx / marketplace | not built, secondary by ADR-D2 |

Two known rough edges in the global channel, both design calls: pruning the
link-*owning* version leaves an installed version with no `bin/scio`; and
uninstalling a version dangles the six category links of any project bound from
it, which `propagation.md` records because the receipt lists installed files
rather than the projects bound from them.
