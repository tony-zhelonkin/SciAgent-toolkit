# Tracer-Bullet: Inject/Eject Path — 2026-06-02

**Runner:** Anton (power-user tinkerer persona)
**Model:** claude-sonnet-4-6
**Toolkit:** `/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit`
**Scratch dir:** `/tmp/claude-788715489/tmp.ui2Rbul3WI` (throwaway; also `/tmp/claude-788715489/tmp.N71u1Pkza6` for clean reproductions)

---

## A. Commands Run (in order)

```bash
# Setup
SCIAGENT=.../bin/sciagent
SCRATCH=$(mktemp -d)
$SCIAGENT new project --type analysis $SCRATCH
cd $SCRATCH && $SCIAGENT activate base

# Discovery
$SCIAGENT list skills
$SCIAGENT list agents
$SCIAGENT list commands
$SCIAGENT list deps muon-multimodal-analysis
$SCIAGENT list deps scenic-r-python-interop
$SCIAGENT list dependents anndatar-seurat-scanpy-conversion
$SCIAGENT list dependents scanpy

# Inject - simple auto-detect
$SCIAGENT inject snapatac2-atac-preprocessing
$SCIAGENT status
# -> injected count in manifest, overlay _injected created

# Inject - orchestrator with requires chain
$SCIAGENT inject muon-multimodal-analysis
# -> INJECTED but deps NOT pulled in (Finding #1)

# Eject
$SCIAGENT eject muon-multimodal-analysis
$SCIAGENT eject snapatac2-atac-preprocessing

# Inject by kind flag
$SCIAGENT inject --skill snapatac2-atac-preprocessing
$SCIAGENT inject --agent architect
$SCIAGENT inject --command architect
$SCIAGENT inject --command map
$SCIAGENT inject --command design

# Inject by tag
$SCIAGENT inject --tag qc            # single skill; already in base, nothing injected
$SCIAGENT inject --tag preprocessing # 2 new injected; 4 already-mounted skipped
$SCIAGENT inject --tag bogus-nonexistent-tag-xyz  # expected hard-fail

# Ambiguity & errors
$SCIAGENT inject architect           # hard-fail: exists as both agent and command
$SCIAGENT inject nonexistent-skill-xyz  # hard-fail: not found
$SCIAGENT inject scanpy              # already stack-mounted, silent no-op (exit 0)

# Eject tests
$SCIAGENT eject --skill cellranger-arc-multiome
$SCIAGENT eject --tag preprocessing   # ejects only via:tag:preprocessing entries
$SCIAGENT eject scanpy               # refuses: stack-mounted
$SCIAGENT eject harmonypy-batch-integration  # not injected → no-op

# Roundtrip: inject → status → eject → status → verify no orphans
$SCIAGENT inject louper-seurat-conversion
$SCIAGENT status  # shows as "injected"
$SCIAGENT eject louper-seurat-conversion
# manifest clean, symlinks gone → no orphans

# BLOCKER: inject then eject of inherited-via-requires skill
$SCIAGENT inject --skill signac-footprint-visualization  # silently accepted
$SCIAGENT eject --skill signac-footprint-visualization   # DESTROYS dep symlink
$SCIAGENT status  # signac is GONE (not even "inherited via requires:")

# Eject ambiguity (injected as two kinds)
$SCIAGENT inject --agent architect
$SCIAGENT inject --command architect
$SCIAGENT eject architect  # hard-fail: "ambiguous — injected as both agent and command"
$SCIAGENT eject --agent architect   # succeeds
$SCIAGENT eject --command architect # succeeds

# Deactivate with injected items
$SCIAGENT inject snapatac2-atac-preprocessing
$SCIAGENT inject louper-seurat-conversion
$SCIAGENT deactivate   # cleans up both injected and base-mounted items → no orphans

# Deactivate overlay by name
$SCIAGENT activate base
$SCIAGENT inject snapatac2-atac-preprocessing
$SCIAGENT deactivate _injected
# -> calls cmd_activate base (clean-slate re-activate, drops injected with warning)

# status --source for different kinds
$SCIAGENT status --source snapatac2-atac-preprocessing  # "injected (into _injected)"
$SCIAGENT status --source scanpy              # "base"
$SCIAGENT status --source signac-footprint-visualization  # BUG: "not in active stack"

# status --effective
$SCIAGENT status --effective  # omits inherited-via-requires skills

# status --json
$SCIAGENT status --json  # inherited-via-requires skills missing from JSON output

# Validate
$SCIAGENT validate  # passes; warnings on allowlisted collisions
```

---

## B. Findings

### [BLOCKER] #1 — `eject` of a skill injected-over-inherited-dep destroys the requires symlink

**File:** `lib/sciagent/eject.sh:_eject_remove_symlinks` + `_eject_drop_injected_entry`

**Root cause:** When `tf-footprint-differential-analysis` (in `base.yaml`) lists `signac-footprint-visualization` under `requires:`, `activate.sh` creates a symlink for it and records it in `manifest.symlinks`. The `_inject_entry_is_stack_mounted` guard in `inject.sh` only checks role YAML (not `requires:` closure), so `inject --skill signac-footprint-visualization` succeeds silently. When `eject` is called, `_eject_drop_injected_entry` unconditionally removes the `.claude/skills/signac-footprint-visualization` path from `manifest.symlinks` and `_eject_remove_symlinks` deletes the file. No check is made whether the symlink was also created by `activate.sh`'s requires-resolution phase. Result: `signac-footprint-visualization` vanishes from the active stack mid-session; `tf-footprint-differential-analysis` is silently broken.

**Reproduction (clean):**
```
$ sciagent activate base
Activated stack: base  skills: 39 …
$ ls .claude/skills/ | grep signac
signac-footprint-visualization       ← present as "inherited via requires:"
$ sciagent inject --skill signac-footprint-visualization
injected: signac-footprint-visualization (into _injected)
$ sciagent eject --skill signac-footprint-visualization
ejected: signac-footprint-visualization (skill)
$ ls .claude/skills/ | grep signac
(empty)                              ← GONE; tf-footprint-differential-analysis broken
$ sciagent status | grep signac
(no output)
```

**Fix sketch:** Before deleting a symlink in `_eject_remove_symlinks`, check whether the path appears in `manifest.symlinks` for a reason OTHER than the inject entry being ejected (i.e., whether activate.sh also created it). If yes, preserve the symlink and only drop the injected manifest entry.

---

### [BUG] #2 — `inject` does NOT pull the transitive `requires:` closure

**File:** `lib/sciagent/inject.sh:_inject_named_entry` (no call to `skill_resolve_transitive`)

`skill_deps.sh` exports `skill_resolve_transitive`, and `activate.sh` uses it to resolve and mount all transitive deps. But `inject.sh` injects exactly one entry and never calls `skill_resolve_transitive`. Result: injecting an orchestrator skill leaves its leaf skills absent.

**Evidence:**
```
$ sciagent list deps muon-multimodal-analysis
multimodal-anndata-mudata
scanpy
snapatac2-atac-preprocessing
harmonypy-batch-integration
atac-differential-accessibility
pyranges-peak-gene-linkage
pygenometracks-coverage-plots
scvi-multivi

$ sciagent inject muon-multimodal-analysis
injected: muon-multimodal-analysis (into _injected)

$ ls .claude/skills/ | grep -E "harmonypy|atac-diff|pyranges|pygenome|multimodal-anndata"
(empty — 5 of 8 deps are missing from the active stack)
```

No warning is emitted. `sciagent validate` does not detect this. The AI sees the orchestrator but lacks the leaf skills it references.

**Fix sketch:** In `_inject_named_entry` for `kind=skill`, call `skill_resolve_transitive` and inject any dep not already mounted (same guard as `activate.sh` Phase B). Or at minimum warn: "Note: `muon-multimodal-analysis` requires 5 skills not in the active stack; run `sciagent inject <dep>` for each."

---

### [BUG] #3 — `status --source` returns "not in active stack" for inherited-via-requires skills

**File:** `lib/sciagent/status.sh:_status_render_source` (lines 421–464)

`_status_render_source` only searches `SKILL_ORDER[]` (declared-in-role skills) and `INJECTED_SKILLS[]`. It does not search the `INHERITED_SKILLS[]` array that `_status_render_text` computes from `MANIFEST_SKILLS[]`.

**Evidence:**
```
$ sciagent status | grep signac
  signac-footprint-visualization  inherited via requires:    ← visible in text

$ sciagent status --source signac-footprint-visualization
not in active stack: signac-footprint-visualization          ← contradicts text
exit code: 1
```

Same for `tobias-footprint-bindetect` and `hint-atac-differential-footprint`.

---

### [BUG] #4 — `status --effective` omits inherited-via-requires skills

**File:** `lib/sciagent/status.sh:_status_render_effective` (lines 411–419)

`_status_render_effective` iterates `SKILL_ORDER[]` and `INJECTED_SKILLS[]` only. The `INHERITED_SKILLS[]` array is computed locally inside `_status_render_text` and never passed to `_status_render_effective`.

**Evidence:**
```
$ sciagent status | grep "Skills ("
Skills (39 effective):                      ← includes 3 inherited skills in count

$ sciagent status --effective | wc -l
36                                          ← 3 inherited skills missing

$ sciagent status --effective | grep tobias
(empty)
```

---

### [BUG] #5 — `status --json` omits inherited-via-requires skills

**File:** `lib/sciagent/status.sh` (JSON rendering section, ~lines 496–565)

The JSON `"skills"` array is built from `SKILL_ORDER[]` and `INJECTED_SKILLS[]`, same as `--effective`. `INHERITED_SKILLS[]` is never added. Text mode says 39 effective; JSON has 36.

**Evidence:**
```
$ sciagent status --json | python3 -c "
import json,sys; d=json.load(sys.stdin)
print(len(d['skills']), 'tobias?' , 'tobias-footprint-bindetect' in [s['name'] for s in d['skills']])
"
36 tobias? False

$ sciagent status | grep "Skills ("
Skills (39 effective):
```

The discrepancy (39 text vs 36 JSON) will silently mislead any tooling that parses JSON to determine effective skills.

---

### [BUG] #6 — `inject` silently accepts an inherited-via-requires skill with no warning

**File:** `lib/sciagent/inject.sh:_inject_entry_is_stack_mounted`

`_inject_entry_is_stack_mounted` calls `role_load` and looks for `SKILL <name>` lines from the role YAML. Skills added via `requires:` are NOT in the role YAML, so the guard returns false and inject proceeds. There is no check against the `MANIFEST_SKILLS` set (the skills actually symlinked).

**Effect:** User injects a skill that is already accessible via requires, creating a redundant manifest entry. Later eject removes the shared symlink (Finding #1). The inject succeeds silently even though the skill is visible and functional already.

**Fix sketch:** In `_inject_entry_is_stack_mounted`, also check whether `manifest_symlinks` contains `.claude/skills/$name`; if yes, output a warning: "already accessible via `requires:` inheritance from <parent>; injecting will allow you to eject it, which WILL break the inheritance chain."

---

### [CONFUSING] #7 — Eject error message appends `/_injected` to all stack-mounted errors

**File:** `lib/sciagent/eject.sh:_eject_named_entry` line 123

```bash
echo "sciagent eject: '$name'$kdesc is part of stack-mounted role '$base'${overlay:+/$overlay}; use 'sciagent deactivate' instead" >&2
```

When the stack is `base _injected` (i.e., any time the user has injected anything), the error reads:

```
sciagent eject: 'scanpy' is part of stack-mounted role 'base'/_injected; use 'sciagent deactivate' instead
```

`scanpy` is in `base`, not in `_injected`. The `/_injected` suffix is always appended when `overlay` is non-empty, even when the skill is only in `base`. This is misleading — it implies the skill is split across both layers.

**Fix:** Track which role the skill was found in during `_eject_entry_is_stack_mounted` and report only that role name.

---

### [CONFUSING] #8 — `deactivate _injected` does a clean-slate re-activate, not a targeted overlay teardown

**File:** `lib/sciagent/deactivate.sh:cmd_deactivate` line 47

```bash
if [[ "$target" == "$overlay" ]]; then
    # Re-activate solo-base. cmd_activate handles teardown-first.
    cmd_activate "$base"
    return 0
fi
```

When you run `sciagent deactivate _injected`, you expect it to collapse the `_injected` overlay and eject all injected items, leaving `base` intact. Instead, it calls `cmd_activate base`, which:
1. Tears down everything (symlinks, manifest, block)
2. Re-activates from scratch
3. Emits the "dropping injected entries" warning

**Actual output:**
```
$ sciagent deactivate _injected
sciagent: warning — activate is a clean-slate operation; dropping injected entries:
  - skill muon-multimodal-analysis
  ...
Activated stack: base
  skills:   39
  agents:   7
  commands: 1
```

This is unexpected — it looks like an activation, not a deactivation. It also discards injected items rather than cleanly ejecting them.

---

### [NIT] #9 — `inject --tag` of a tag where all matching skills are stack-mounted exits 0 with only stderr messages

**File:** `lib/sciagent/inject.sh:_inject_by_tag`

```
$ sciagent inject --tag qc
already mounted via stack-role 'base': single-cell-rna-qc — nothing to inject
```

This exits 0, which is correct (idempotent). But the message goes to stderr with no stdout acknowledgement. There's no summary like "0 new injections (1 already mounted)". A script caller relying on stdout would see nothing.

---

### [NIT] #10 — Companion-skill note only fires for agent/command inject, not the reverse

**File:** `lib/sciagent/inject.sh:_inject_named_entry` (lines 228–235)

When you `inject --agent architect`, you get:
```
note: companion skill 'architect' available — `inject --skill architect` to add
```

But there is no equivalent when injecting a skill that has a same-name agent (e.g., if a skill and agent share a name). The discoverability is one-directional only.

---

### [ORPHAN] #11 — Deactivate with stale-symlink state emits warning, not error

**File:** `lib/sciagent/symlinks.sh:symlink_teardown_all` line 261

When a symlink listed in `manifest.symlinks` is already absent (e.g., from the Finding #1 eject bug), `symlink_teardown_all` emits:
```
warning: .claude/skills/signac-footprint-visualization already gone, skipping
```

This is benign in isolation but indicates the manifest and filesystem diverged. There is no way to detect or repair this short of a full deactivate/reactivate cycle. A `sciagent validate` or `sciagent status` could check for manifest↔filesystem drift, but currently do not.

---

## C. Verdict

| Area | Rating | Notes |
|------|--------|-------|
| **Core inject/eject roundtrip** | GOOD | No orphan symlinks or manifest entries for normal usage |
| **Dep-closure on inject** | BROKEN | Orchestrator skills injected without their leaf deps; no warning |
| **Eject of inherited-dep** | BROKEN | Destroys shared symlink; breaks requires chain mid-session |
| **status --source for inherited skills** | BROKEN | Returns "not in active stack" even though they are mounted |
| **status --effective / --json** | BROKEN | Omits inherited-via-requires skills; count mismatch vs text |
| **Ambiguity handling** | GOOD | Hard-fails with clear messages for all cases |
| **Unknown name / bad tag** | GOOD | Clear error messages, exit 1 |
| **Stack-mounted guard** | PARTIAL | Correct for YAML-declared items; misses requires-inherited items |
| **Deactivate cleanup** | GOOD | Injected items fully cleaned up on full deactivate |
| **deactivate <overlay-name>** | CONFUSING | Calls clean-slate re-activate instead of targeted teardown |
| **Tag inject/eject roundtrip** | GOOD | `via:tag:<name>` correctly tracked; eject --tag only removes tag-via entries |
| **Error messages** | CONFUSING | Stack-mounted error appends `/_injected` even when skill is only in base |
| **Overlay disclosure** | GOOD | "injected: X (into _injected)" is clear |
| **Usability overall** | 6/10 | The inject/eject UX is solid for simple cases; breaks badly when orchestrator skills or inherited deps are involved |

### Priority fixes
1. **[BLOCKER #1]** Guard `_eject_remove_symlinks` against destroying requires-inherited symlinks
2. **[BUG #2]** `inject` should resolve and mount the transitive `requires:` closure (or at minimum warn)
3. **[BUG #3/4/5]** Include `INHERITED_SKILLS[]` in `--source`, `--effective`, `--json`
4. **[BUG #6]** Warn (or refuse) when injecting an already-requires-inherited skill
