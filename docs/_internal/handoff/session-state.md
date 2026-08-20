# Session state — cold-start handoff

Living file, overwritten rather than accumulated. Read it first in a cold
session, then **re-derive the numbers** with the commands in §2 before acting on
them: every figure here was measured, and every figure here decays.

Authoritative design lives in `docs/architecture.md` and
`docs/propagation.md`; decisions in
`docs/proposals/2026-08-11-offline-distribution/50_ADRs.md`. This file carries
only what those cannot: where the work stopped, what is waiting on the owner,
and the facts that cost time to rediscover.

Last verified: **2026-08-20**, toolkit at `287e8e4`.

---

## 1. The shape of the thing

Four layers, and the split between them is the design:

- **THE ASSET** — `skills/` 85, `agents/` 21, `commands/` 7, `craft.yaml`.
  Non-regenerable: it encodes traps (the CORESH Entrez-integer trap,
  species-mismatch silent failures, iterative-peak-merging) that took real
  analyses to learn.
- **THE MECHANISM** — `lib/scio/` 2,619 LOC + `bin/scio` 65. Three verbs, one
  concern each: `link` = availability (convergent), `craft` = standing text
  (idempotent), `lint` = enforcement (read-only).
- **THE AUTHORITY** — Git. The submodule pin **is** the lock (ADR-D1). There is
  no `scio.lock` and there should not be one.
- **THE RECORD** — `docs/` + 9 ADRs (D1–D9).

Two structural facts that constrain every future change:

- **Six directory links, not 579 file mounts.** Ownership tracking *dissolved*
  rather than being solved — no manifest, no `symlinks.sh`. The price: a
  directory symlink can neither filter nor flatten, so **source shape must equal
  mount shape**. That is why `templates/skill/` exists, why `_attic/` is
  top-level, and why `agents/` is flat. Any change that puts non-mountable
  content inside `skills/`, `agents/`, or `commands/` breaks a consumer.
- **`lint.sh` owns every convention it can measure.** A claim stated twice
  drifts, so CRAFT names the standing boundary, a router skill selects the
  guide, `references/` owns qualitative doctrine, and `lint` owns the
  predicates. If you find a convention asserted in prose with no predicate
  behind it, that is a gap, not a style choice — `287e8e4` closed the last known
  one (the CRAFT line budget).

---

## 2. Re-derive the state before trusting it

```bash
cd <toolkit>                       # scbio-docker/toolkits/SciAgent-toolkit
git log --oneline -1 && git rev-list --count hub/dev..dev
bash tests/run-all.sh | tail -3
bin/scio lint --check toolkit --strict && echo clean

# Fleet census. Requires each directory to BE its own repo root: `git -C` on an
# empty submodule directory silently walks up to the superproject and reports
# ITS HEAD, which reads as a diverged toolkit copy that does not exist.
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

### What those commands showed on 2026-08-20

| | Value |
|---|---|
| `dev` | `287e8e4`, **4 ahead** of `hub/dev` = `origin/dev` = `106f59f` |
| `main` | `c83dfe2` on both remotes |
| Tests | 62 passing, 0 failing |
| `toolkit` lint | clean |
| Tree | clean |

Unpushed, oldest first: `7e9e72a` (route analysis code through convention
guides), `6cfe500` (anchor analysis code in durable conventions), `42ae8a1`
(prune an unused release without reporting damage), `287e8e4` (bind the CRAFT
block to its declared line budget).

Fleet — **24 real copies**:

| Count | Commit | Behind `dev` | What |
|---|---|---|---|
| 1 | `287e8e4` | 0 | scbio-docker — the canonical dev checkout |
| 1 | `106f59f` | 4 | 14616-DM — the swept pilot |
| 20 | `5e5347e` | 52 | unswept (15 fleet-managed + 5 frozen) |
| 1 | `cf19c6d` | 118 | PanSci — PINNED, deliberately |
| 1 | `fb6012a` | not in canonical history | Gama_Vivian — excluded by config |

`scbio-docker` pin drift: records `3ab3768`, checked out `287e8e4` (`^4 PIN`).

---

## 3. Traps — the expensive-to-rediscover facts

**The `module-vendor` state column says nothing about currency.** Of the 15
fleet-managed unswept copies, 5 report `BEHIND v48` and 10 report `OK clean`,
yet all 20 sit at the identical commit. The 10 look clean because they track an
`origin/dev` whose ref is an *ancestor* of their own HEAD. Always verify by SHA.

**`v48` and `52` are both correct.** The behind-count is measured against the
copy's upstream (`origin/dev` at `106f59f`); the true distance to `dev` is 52.
The difference is exactly the four unpublished commits.

**Six copies are missing from `module-vendor status` on purpose.** `config.sh`
`MODULE_VENDOR_EXCLUDE` names five copies **frozen by owner decision
2026-08-18** (finished single-project analyses keeping the pre-demolition
11-verb toolkit: 12868-EH, Dakota_Black_NPC_DRP1, 13403-YD_Christina,
JBader_scHFD, 13036-DM) plus `Gama_Vivian_DRP1_bulkRNAseq`. Exclusion rather
than branch-pinning is what holds them, because `sync align` advances any copy
it can see that is on a shared branch. **Do not "fix" their absence.**

**`AdaW-eWAT-WL-bulkRNAseq` is benign.** It records a gitlink at pin `108da7a`
with an empty `01_modules/SciAgent-toolkit` — an uninitialized submodule. No
mounts, no `AGENTS.md`, nothing dangling. Invisible to `module-vendor` by
construction, since there is no repo to discover.

**Never run `link` between `align` and `rename` during the sweep.** Align
dangles 8 agent and 2 command mounts per repo; linking in that window
materializes the damage instead of converging it.

**`docs/_internal/` is gitignored** (`.gitignore:11`). The files tracked there
were force-added, so this handoff and its siblings need `git add -f`. A plain
`git add -A` silently skips them.

**Adding a key to `craft.yaml` must not perturb the rendered body.** Every
consumer stores a SHA1 of the CRAFT block; a changed body shows up as drift in
all of them on re-pin. Verify byte-identity against `HEAD:craft.yaml` before
committing a `craft.yaml` change. Also: `body:` **must** stay the last top-level
key — `block.sh` reads it as a trailing block scalar.

**Release tests build throwaway fixture repos, never this one.** The builder
runs the suite, so a test that built this repo would re-enter it. This is also
why there is deliberately **no `--skip-checks` flag**: an escape hatch on a
release gate eventually gets used for a real release.

**`lint` has 11 selectable names but `all` runs 10.** `toolkit` is opt-in, for
maintainers and the release gate; it validates this repo's own assets, so it
must never gain a finding that fires in a consumer.

**Bash only on the core path** — no `yq`, `jq`, or `python3`. Non-negotiable
(`AGENTS.md` rule 1).

---

## 4. Blocked on the owner

1. **Review `106f59f..dev`** (4 commits). The last two are independent of the
   craft/skill work. `nvim -c 'DiffviewOpen 106f59f..dev'`. The one line to read
   closely is the config invariant in `restartability-and-dataflow.md:72–76`.
2. **Commit 14761-DM.** Was blocked on a live `codex exec` writing
   `02_analysis/helpers/grn_viz/` — re-check whether it is still running. It
   stays out of the sweep until committed, and it is the first real instance of
   the earned-by-cohesion helper family, worth citing in the reference once
   settled.
3. **ADR-D9 rename**, owner's call whether it runs before the sweep. GitHub repo
   → `scio`, vendor path → `01_modules/scio/`. Four vendor-dir spellings exist
   live and each must be resolved individually: `01_modules`, `01_scripts`,
   `01_Modules`, `01_Scripts`.
4. **Task #36** — whether the CRAFT block gets a token-weight budget. The cap
   now binds shape (17 lines of 25) but not cost: those 17 lines are 4,658
   characters and the longest bullet is 693. A per-bullet cap would fail today.
   Picking the number means editing the owner's standing text.

Nothing below is worth starting before item 3 is decided: the sweep touches
vendor paths in 15 repos, and `scbio-docker/AGENTS.md` would need fixing twice.

---

## 5. Next steps, in order

5. Merge `dev --no-ff` → `main`; push hub **and** origin.
6. **The sweep, 15 repos** (#30) — align → rename → link → craft → commit →
   bump. Repoint the stale `origin/dev` upstreams to `hub/dev` in the same pass
   (#31).
7. Clear the `scbio-docker` pin drift off `3ab3768`.
8. dev-env adoption of the 8-file handoff bundle under
   `docs/_internal/handoff/dev-env/` (#33) — this is also what makes the
   interaction stance live; it is active nowhere today. dev-env is clean at
   `0b5a4b0`.
9. `scbio-docker/AGENTS.md` documents a two-step `sciagent new project --type
   analysis <dir>`. That verb is deleted, so the file currently sends an agent
   to a command that does not exist. Owed to `init-project.sh` + the three
   verbs.
10. Consumer `scripts/` → `stages/` (#8). Last, after the repos rename.

## 6. Owner-only, no deadline

- `docs/guidelines/` remainder — `data_processing.md` 281, `gsea_analysis.md`
  342, `master_tables.md` 380 lines await reconciliation against the assay
  skills. That call is scientific, not structural.
- A decision-gates lint check (#34) — `notebooks/` is invisible to lint today.
  Deferred pending a stable contract.
- ADR-D5 (#6) decided, unimplemented: extract project binding behind a common
  layer plus harness adapters. ADR-D10 unwritten.
- `v-last-orchestrated` (= `5e5347e`, exactly where the unswept copies sit) is
  hub-only. Reachable as an ancestor, so nothing is lost, but the named anchor
  is not on origin.
- Stale `.sciagent/manifest.json` (mode 0600) inside the 14616-DM and
  13403-YD_Christina toolkit checkouts — what makes 14616-DM read
  `UNTRACKED ?1`. Deliberately untouched: `rm` inside a consumer is the owner's
  call.

## 7. Distribution channels, for orientation

| Channel | State |
|---|---|
| Vendored submodule (canonical) | shipped, proven end-to-end |
| Global install — `install.sh` + `scripts/build-release.sh` | shipped, exercised end-to-end |
| npx / marketplace | not built, secondary by ADR-D2 |

Two known rough edges in the global channel, both design calls rather than bugs:
pruning the link-*owning* version leaves an installed version with no
`bin/scio` (fixing it means uninstall acquires a version-promotion policy); and
uninstalling a version dangles the six category links of any project bound from
it, which `propagation.md` now records because the receipt lists installed files
rather than the projects bound from them.
