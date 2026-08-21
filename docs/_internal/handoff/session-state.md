# Session state — cold-start handoff

Living file, overwritten rather than accumulated. Read it first, then
**re-derive the numbers** with §2 before acting on them: every figure here was
measured, and every figure here decays.

Authoritative design: `docs/architecture.md`, `docs/propagation.md`. Decisions:
`docs/proposals/2026-08-11-offline-distribution/50_ADRs.md`. This file carries
only what those cannot — where the work stopped, what waits on the owner, and
the facts that cost time to rediscover.

Last verified: **2026-08-21**, toolkit at `686c408`. 65 tests passing.

---

## 1. The shape of the thing

- **THE ASSET** — `skills/` 85, `agents/` 21, `commands/` 7, `craft.yaml`.
  Non-regenerable: it encodes traps (the CORESH Entrez-integer trap,
  species-mismatch silent failures, iterative peak merging) that took real
  analyses to learn.
- **THE MECHANISM** — `lib/scio/` 2,896 LOC + `bin/scio` 65. Three verbs, one
  concern each: `link` = availability (convergent), `craft` = standing text
  (idempotent), `lint` = enforcement (read-only).
- **THE AUTHORITY** — Git. The submodule pin **is** the lock (ADR-D1). No
  lockfile, and there should not be one.
- **THE RECORD** — `docs/` + 10 ADRs (D1–D10).

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

`dev` = `686c408`, **27 ahead** of `hub/dev` = `origin/dev` = `106f59f`.
`main` = `c83dfe2` both remotes. **65** tests passing, 0 failing. `toolkit` lint
clean. Tree clean.

A file cannot name the commit that adds it, so the tip is this handoff's own
commit, one past `686c408`. Re-derive rather than trusting either number.

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
03d6b0a  Reconcile the handoff's figures with the tip it describes
60c1f2d  Stop naming memory locations the mechanism cannot deliver
3a6870e  Ship no directory that has nothing to copy
432e9cf  Keep project paths out of the user-global habit layer
779eb6f  Record ADR-D10 — memory mirrors the analysis, lint is its surface
710094f  Bring the handoff up to a fully implemented plan
4503fbe  Write category links that resolve from both sides of the container
d1a8dc4  Track the plan directories a clone could not see
4b150fd  Record #37 fixed and the sweep unblocked
4e67fe9  Record the consult on whether the memory tree earns its own depth
3f6d4ef  Let scope be the only structure in the memory tree
a4ba98e  Trace where the delegation time actually went
63f27dc  Bring nine catalog assets onto the one memory grammar
72dff47  Record the publication-topology and delegation-redesign consults
686c408  Count both forms of history the memory tree can already have
```

Six change behaviour: `42ae8a1` (uninstall exit-3), `287e8e4` (CRAFT budget
check), `2a1c702` (delegation assets), `60c1f2d` (a deleted lint predicate plus
the new `internal-memory` check, and a CRAFT body change every consumer sees as
drift on re-pin), `3a6870e` (five `.gitkeep` and six READMEs leave the
scaffold), `4503fbe` (**relative category links** plus `harness-links`; every
bound copy's mounts get rewritten by its next `link`), `3f6d4ef` (**the memory
grammar flattens**; `internal-memory` rewritten, and a second CRAFT body change),
`63f27dc` (nine catalog assets re-routed), `686c408` (a history predicate that
counts worktrees and parent-tracked trees). The rest is documentation, planning
and evidence.

In **scbio-docker**, on branch `feat/bulkirna-v0.5.0`: `5dd9cbd` (seven docs off
the deleted verb) and `127c8a5` (the `si` alias — see §7). The submodule pin
drift is deliberately untouched.

Fleet — **25 real copies, 19 tracked** (`JR-MC-Tonsill/JR-MC` appeared
2026-08-21):

| Count | Commit | Behind `dev` | What |
|---|---|---|---|
| 1 | `686c408` | 0 | scbio-docker — the canonical dev checkout |
| 1 | `e33b635` | 18 | JR-MC — new, bound in-container |
| 1 | `106f59f` | 27 | 14616-DM — the swept pilot |
| 20 | `5e5347e` | **75** | unswept (15 fleet-managed + 5 frozen by decision) |
| 1 | `cf19c6d` | 141 | PanSci — PINNED, deliberately |
| 1 | `fb6012a` | not in canonical history | Gama_Vivian — excluded by config |

`scbio-docker` pin drift: records `3ab3768`, checked out `686c408`.

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

**Scratch escapes the project.** `craft.yaml` called root `_scratch/` "the only
sanctioned throwaway zone"; **no project has one**. Agents write to
`/tmp/claude-<pid>/<hashed-workspace>/<uuid>/scratchpad/`. `delegate-cli`
itself prescribed `/tmp` three times — prescribed behaviour, not drift. Both
halves are now closed: the skill by `2a1c702`, the claim by `60c1f2d`, which
deletes it rather than relocating it. ADR-D10 decision 3 records why.

**What is still open of the five.** The memory tree is still gitignored in every
consumer, so the tracked-files instruction is still unbacked there — ADR-D10
records the nested-repo-plus-pointer topology and deliberately does not
implement it. Availability still needs the sweep. The other three are closed.

---

## 4. Blocked on the owner — six decisions

1. **Review `106f59f..dev`** (20 commits; 6 change behaviour). The owner intends to tweak wording
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
   lines) but not cost (4,742 chars; longest bullet 693). A per-bullet cap fails
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

## 4b. #37 — FIXED 2026-08-21 in `4503fbe`. The sweep is unblocked.

Owner ruling: relativize **and** add the lint predicate, so it cannot come back
silently. `_link_symlink_target` now serves both the six category links and the
two helper libs — relative while the source is inside the project, absolute for
a global install, where a relative chain breaks the moment the project moves.

The part that mattered most: the idempotency test compares the **written** target
rather than where it resolves. A resolution test reads an absolute link as
current, so `link` would have been a silent no-op in all 25 bound copies — the
fix meant to reach them would never have arrived. A re-run now reports each
repair. `lint --check harness-links` (in `all`) fails a mount that does not
resolve or that names an absolute path inside the project, and stays quiet on an
absolute target outside it. `tests/test_lint_harness_links.sh` covers all nine
cases, including the repair path and the global-install exemption.

`CHANGELOG:59` is true again: its historical claim that all three namespaces and
both mirrors are relative now matches behaviour.

**What the sweep must still do about this:** every one of the 25 copies is
carrying absolute links today, so the sweep's `link` step is now also the repair
step — and it must run *after* the re-pin, or it writes relative links against
the old pin. Verify with `readlink .claude/skills` per repo, not by resolution.

### The original finding, kept for the record

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

| The toolkit says | The mechanism does | State |
|---|---|---|
| durable memory lives in tracked files, in `docs/_internal/` | `link` gitignores that path (`link.sh:278`) | open — ADR-D10 topology unimplemented |
| `_scratch/` is the only sanctioned throwaway zone | no project has one | claim deleted `60c1f2d` |
| "`docs/_internal/` missing — run: `scio link`" | `link` never creates it | predicate deleted `60c1f2d` |
| AGENTS.md: read the skill | activation mounted 60 of 86, silently | fixed by design; needs the sweep |
| pass `--search` to codex | the flag exists in no live version | fixed `2a1c702` |
| category links are relative (`CHANGELOG:59`) | they were absolute | fixed `4503fbe` |
| `alias si` runs the toolkit CLI | it named `bin/sciagent`, deleted | fixed `127c8a5` |

**The review question that falls out, worth applying to anything this toolkit
asserts: what enforces this claim, and can that thing actually deliver?**

Seven instances now, six closed. The one that remains is the one where the
mechanism, not the wording, has to change: the memory tree is still gitignored

---

## 6. The implementation run — shipped 2026-08-21

`docs/_internal/plans/2026-08-20-memory-and-seams/` — `00_INDEX.md` plus eight
phase briefs. **Owner approved the design; all eight phases are now landed.**
Read `00_INDEX.md` for the reasoning; §6 lists what stayed out of scope, and it
still holds.

| # | Ships | Repo | State |
|---|---|---|---|
| 01 | delete the three unbacked claims (lint predicate, `_scratch/` claim, CRAFT routes, hook text) | scio | **DONE `60c1f2d`** |
| 02 | the opt-in `internal-memory` lint check that replaces them | scio | **DONE `60c1f2d`** |
| 03 | stop shipping `.gitkeep` and empty category dirs | scio | **DONE `3a6870e`** |
| 04 | strip scio's path grammar from the user-global dev-env template | scio | **DONE `432e9cf`** |
| 05 | 7 docs citing the deleted `sciagent new project` verb | scbio-docker | **DONE `5dd9cbd`** |
| 06 | ADR-D10 — the record for all of it | scio | **DONE `779eb6f`** |
| 07 | delegation seam: retention contract, cite scbio-docker for the bwrap cause | scio | **absorbed by 08** |
| 08 | `probe.sh`/`launch.sh` assets; SKILL.md 416 → 150 | scio | **DONE `2a1c702`** |

**The plan is fully implemented.** 01→02 ran in this session's own lane; 03, 04,
05 and 06 ran as four backgrounded `gpt-5.6-sol` workers launched through the
phase-08 assets themselves — four short invocations, 5.3–6.1 KB prompts through
stdin, one lock per unit, `--parallel-ok` making disjointness an explicit claim.
Every gate was re-run here rather than taken from a worker's report.

Phase 07's two remaining items were already shipped inside 08: `SKILL.md:33`
carries the retention contract and `:65` cites
`scbio-docker/docs/ai-integration.md` for the seccomp cause. Only its owner
decision survives (§4 item 5).

Two things the workers found that the plan had not scoped, both the same defect
class:

- **`setup_ai_env.sh` installed `alias si=".../bin/sciagent"`** — a path with no
  file behind it since the CLI rename, so `si` was broken in every container it
  provisioned. This is why the owner types full vendor paths. Fixed in
  `127c8a5`, and the marker moved to `_scio_si_alias` so a shell already holding
  the stale block receives the corrected one and the later `alias` wins; keeping
  the old marker would have left every live container broken.
- **`init-container.sh` printed the deleted verb as its own "Next steps"** —
  a stale instruction at the moment the reader acts on it.

Phase 03 also went past its brief in the right direction: having deleted the
directories, it rewrote the four templates that pointed at them, so no
instruction outlived its target. Two things it and 04 missed, fixed here: the
dev-env template still cited `SCIAGENT:ROLES`, a block this toolkit deleted, and
ADR-D10 said nothing about repos with no stages, which the software templates now
key as `docs/_internal/<work-stem>/`.

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

### 6b. After the plan — the grammar got simpler and the delegation gap got named

**The memory tree lost a level** (`3f6d4ef`, `63f27dc`). A consult asked to
default every directory to *delete* and argue each back in. Scope survived;
document category did not. A stage holds `session.md` and flat topic notes beside
it — `30_grn/network-selection.md`, not `30_grn/reasoning/network-selection.md`.
Software repos lose the invented `<work-stem>` key and use `_project/` alone,
because there is no observable work key for lint to check. The `_common` README
template is gone: it stated the grammar a third time. Owner rulings, all four:
flatten, `_project/` for software, delete the README, keep Decision 6 with lint
observing it.

**"Git history is its archive" is deleted.** It is false in an ignored,
non-nested tree, where an in-place `session.md` update destroys what it replaced.
ADR-D10 still recommends the nested repository, **no verb creates one**, and
`internal-memory` reports a tree holding continuity records with no history. A
`.git` *file* (worktree, submodule) and a parent that force-tracked the tree both
count as history — the first predicate missed both (`686c408`).

**The grammar was stated in twelve places.** Changing it once changed nothing:
nine catalog assets still routed to `docs/_internal/reasoning/`, including
`agents/handoff.md`, which emitted *"MISSING TRACE — add to
docs/_internal/reasoning/"* into real projects. Fixed. **Five more remain and
need a decision, not a substitution — task #45.** `/decompose` and `/implement`
carry a whole second plan vocabulary (`README.md`, `phase-NN.md`, `_campaign.md`)
that the `00_INDEX.md` predicate will flag on every plan they produce.

**Delegation: the failure was delivery, decisively** — the trace is
`docs/_internal/research/2026-08-21-delegation-monitoring/TRACE_jr-mc-session.md`,
measured from a 5.2 MB transcript. JR-MC ran the pre-asset 416-line skill;
`assets/` was absent from its pinned checkout. `2a1c702` existed in this
repository for **13 h 25 min** before the owner complained and never arrived,
because nothing re-pins a vendored submodule mid-session.

Dead time is *entirely* launch-to-observation: **349.8 min** worst case, **219.7
min** with the owner demonstrably present, against 4 s from observation to action.
The split is mechanical — every gap over six minutes was a `nohup &` with no
watcher; every gap under two minutes was a run the harness knew about. `| tail -5`
was inert, measured: **0 bytes** observable after ten minutes while the run log
held **425,956**.

Three findings that outlive that project:

- **Prose warnings do not work.** The old skill documented the `pgrep -f`
  self-match and the `pkill -f` exit-144 trap, both marked verified. The session
  reproduced both — once with the truth and the falsehood in the *same command
  output*, and the falsehood won for four hours. This is the argument for assets
  over prose, made by experiment.
- **The orchestrator cannot wake itself.** All nine idle gaps end with a user
  message or a harness notification; there is no third mechanism. A detached run
  with no harness task id is unobservable *in principle* — so
  **`launch.sh --bg` reproduces the `nohup` failure** and must go (task #46).
- **The owner's foreground `!` runs were the well-behaved ones**, because timeout
  promotion gave them a task id. The deliberate, correct-looking backgrounding was
  the invisible half.

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

**`docs/_internal/` is gitignored** (`.gitignore:11`, a bare `_internal/`).
Everything tracked there was force-added, so this file and its siblings need
`git add -f`. Measured 2026-08-21: **45 tracked, 88 present-but-untracked** in
this repo's own memory tree — the dev-env bundle, this plan and its two research
directories are tracked; eight earlier plan directories, `refactor-20260602-plan`,
`tracer-20260602-checks` and 30 `codex-step*.{log,report.md}` are not. The logs
are correctly disposable under ADR-D10; the plan directories are the third
durable form, invisible to a clone. **The toolkit carries the same memory defect
it is fixing in its consumers**, and that bears directly on §4 item 2.

**`git check-ignore` lies about tracked paths** unless you pass `--no-index`: it
skips anything in the index, so it reported `docs/_internal/` as *not* ignored
while `git add` refused the directory. Two commands, opposite answers, both
right.

**Two doors reach every skill, and the wrong one is more discoverable.** In the
vendored channel a project holds the whole catalog twice: once behind
`.claude/skills/` (the symlink the harness knows) and once at
`01_modules/SciAgent-toolkit/skills/` (a real directory `ls`, glob and `rg` all
find first). Observed in JR-MC 2026-08-21: an agent asked for a skill by name and
read `01_modules/.../louper-seurat-conversion/SKILL.md` with two `sed` ranges —
the whole 6,610-byte file into context, paying the skill's full cost while
bypassing the Skill tool's routing and its `assets/`. The bytes were identical,
so nothing broke; what it cost was progressive disclosure, and what it teaches
the agent is a path that is wrong under a global install and wrong after ADR-D9.
The mechanism can guarantee that the link is relative and resolves; it cannot
guarantee which of two valid paths an agent types. **The global-install channel
is the only one with a single door** — the link points outside the repo, so no
second copy is visible. That is an argument for it beyond offline use.

**A `craft.yaml` edit changes every consumer's CRAFT hash.** Verify the rendered
body against `HEAD:craft.yaml` before committing unless the change is *intended*
to alter it. `body:` must stay the last top-level key — `block.sh` reads it as a
trailing block scalar. `max_lines: 25` is enforced by `lint --check toolkit`.

**Release tests build throwaway fixture repos, never this one** — the builder
runs the suite, so a test building this repo would re-enter it. Hence no
`--skip-checks` flag by design.

**`lint` has 13 selectable names; `all` runs 11.** (`harness-links` joined `all`; `internal-memory` and `toolkit` are opt-in.) `toolkit` and
`internal-memory` are both opt-in: `toolkit` validates this repo's own assets,
and `internal-memory` would fire in every consumer before any project has
adopted the skeleton. Neither may ever gain a finding that fires from `all`.

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
