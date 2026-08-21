# Tracer-bullet check 02 — Role stacking (depth-2 lifecycle)
**Date:** 2026-06-02  
**Model:** claude-sonnet-4-6  
**Persona:** Anton (chromatin/regulatory specialist)  
**Scratch project:** `/tmp/claude-788715489/tmp.n8N67nrXfn/chromatin-analysis`  
**Toolkit:** `/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit`

---

## Commands run (in order)

```bash
# 1. Scaffold
sciagent new project /tmp/.../chromatin-analysis

# 2. Inspect overlays
sciagent list role scatac-regulatory
sciagent list role pathway-signature
# Read roles/scatac-regulatory.yaml, roles/pathway-signature.yaml
# Verified all skill directories exist on disk (both roles: 100% OK)

# 3. Activate depth-2 stack
sciagent activate base scatac-regulatory
# Inspected: .claude/{skills,agents,commands}, .agents/{skills,agents,commands}
# Read .sciagent/manifest.json
sciagent status
sciagent status --effective      # Exit code: 1 (BUG - see Finding 1)
sciagent status --source anndata
sciagent status --source scenic-grn-inference
sciagent status --source tobias-footprint-bindetect

# 4. Depth cap tests
sciagent activate base scatac-regulatory pathway-signature  # 3-role cap test
sciagent activate pathway-signature                          # single-arg on 2-deep stack (Finding 2)

# 5. Overlay switch
sciagent activate base scatac-regulatory   # restore
sciagent activate base pathway-signature   # direct 2-arg switch

# 6. Deactivation
sciagent deactivate pathway-signature      # partial (overlay only)
sciagent deactivate                        # full teardown
sciagent deactivate                        # when nothing active
sciagent deactivate scatac-regulatory      # specific role when nothing active

# 7. Re-activate + idempotency
sciagent activate base scatac-regulatory   # after full teardown
sciagent activate base scatac-regulatory   # second time (idempotency)

# Additional edge cases
sciagent activate base nonexistent-role    # nonexistent role
sciagent activate scatac-regulatory        # overlay-only as sole role
sciagent activate scatac-regulatory pathway-signature  # non-base in first position
sciagent deactivate base                   # deactivating base while overlay active
```

---

## Findings

### FINDING 1 [BUG] `status --effective` exits with code 1 when no injected commands exist

**Location:** `lib/sciagent/status.sh:418`

**Reproducer:**
```bash
sciagent activate base scatac-regulatory
sciagent status --effective; echo $?
# → prints all skill/agent/command names, then: 1
```

**Root cause:** `_status_render_effective()` iterates six arrays in sequence. The last loop is:
```bash
for n in "${INJECTED_COMMANDS[@]:-}"; do [[ -n "$n" ]] && echo "$n"; done
```
When `INJECTED_COMMANDS` is empty, the `:-` expansion yields `""`, the loop body executes once with `n=""`, `[[ -n "" ]]` evaluates false (exit 1), and that becomes the function's — and hence the process's — exit code.

**Bash confirmation:**
```bash
bash -c 'declare -a EMPTY=(); for n in "${EMPTY[@]:-}"; do [[ -n "$n" ]] && echo "$n"; done; echo "Exit: $?"'
# → Exit: 1
```

**Impact:** Any script using `sciagent status --effective` in a pipeline (`if sciagent status --effective | grep -q foo`) will see spurious failure. This is the machine-readable output mode most likely to be used in CI/hooks.

**Fix:** Add `return 0` at end of `_status_render_effective`, or use `|| true` on the last loop, or change the loop to avoid the false-exit pattern.

---

### FINDING 2 [CONFUSING] Single-arg `activate` on an existing 2-deep stack silently nukes the base

**Observed behavior:**
```bash
sciagent activate base scatac-regulatory   # Stack: base + scatac-regulatory (43 skills)
sciagent activate pathway-signature        # User intends: "switch overlay"
# → Activated stack: pathway-signature (8 skills)
# → base is gone
```

**Expected (user mental model):** When a 2-deep stack exists, a 1-arg activate should either (a) refuse and tell you to run `activate base pathway-signature`, or (b) infer "replace the overlay, keep the base." The current behavior — silently wiping the base — is the worst outcome because the user gets a stripped-down harness with no warning.

**What the README says:** `activate <base> [overlay] — Activate role(s); replaces current stack` — technically accurate but the word "replaces" is easy to miss, especially since the 2-arg form *also* replaces.

**Contrast with 2-arg switch (which works cleanly):**
```bash
sciagent activate base pathway-signature   # base+scatac → base+pathway-signature: clean
```

**Missing:** No `activate --overlay <name>` shorthand, no warning when a 1-arg activate implicitly drops an existing base.

---

### FINDING 3 [BUG] Inherited (`requires:`) skills absent from `status --effective` and `status --json`, but present on disk and in text `status`

**Three inherited skills (tobias-footprint-bindetect, hint-atac-differential-footprint, signac-footprint-visualization) are activated by `tf-footprint-differential-analysis`'s `metadata.requires:` chain.**

| Source | Count | Includes inherited? |
|--------|-------|---------------------|
| `status` (text) | 43 | YES — labeled "inherited via requires:" |
| `status --effective` | 40 names | NO |
| `status --json` `.skills` array | 40 | NO |
| disk `.claude/skills/` | 43 symlinks | YES |
| manifest.json `.symlinks` | 43 | YES |

**Root cause:** `_status_render_effective()` only iterates `SKILL_ORDER` and `INJECTED_SKILLS`. `INHERITED_SKILLS` (computed in `_status_load_state` for the text renderer) is a local array not accessible outside that scope, and `_status_render_effective` never builds its own inherited list.

**Impact:** The `--effective` flag is the programmatic interface. An external tool consuming it would see 40 skills but 43 symlinks on disk — unexplained gap. The `status --source tobias-footprint-bindetect` also returns "not in active stack" (exit 1), even though the skill *is* mounted:

```bash
sciagent status --source tobias-footprint-bindetect
# → not in active stack: tobias-footprint-bindetect   (exit 1)
```
This is wrong: the skill IS active (symlink exists, manifest entry exists, AGENTS.md lists it).

---

### FINDING 4 [CONFUSING] `deactivate <overlay>` message is misleading when deactivating base with overlay present

**Observed:**
```bash
sciagent activate base scatac-regulatory
sciagent deactivate base
# → deactivated (removed base implies overlay too)
```

The parenthetical is helpful but easy to miss in scripted use. It would be clearer as a warning line on stderr: `Warning: removing base role also removes overlay scatac-regulatory`.

**Exit code:** 0 (correct).

---

### FINDING 5 [CONFUSING] No guidance on overlay-swap workflow in `--help` or `status` output

**Observed:** After `activate base scatac-regulatory`, there's no hint in the help text or status output about how to swap the overlay. The correct incantation is `activate base pathway-signature` (full 2-arg form), but a user who doesn't know this will try `activate pathway-signature` (1-arg) and lose their base (Finding 2).

**In `sciagent --help`:** `activate <base> [overlay]   Activate role(s) (max stack depth 2)` — terse, doesn't mention that it's a full replace, not an additive operation.

**Suggested addition:** A "Tip" in `status` output when an overlay is active: `To switch overlay: sciagent activate base <new-overlay>`.

---

### FINDING 6 [NIT] `list role <name>` skill count does not match activate-time count

```bash
sciagent list role scatac-regulatory
# → Skills (13): anndata,anndatar-seurat-scanpy-conversion,...
sciagent activate base scatac-regulatory
# → skills: 43
```

**Why:** `list role` counts only entries in the role's own `skills:` YAML array (13 for scatac-regulatory). `activate base scatac-regulatory` merges base (36 skills) + overlay (13, some shadowing), then adds 3 inherited via `requires:` = 43. The numbers refer to completely different things, with no label explaining the discrepancy.

**Impact:** First-time user sees "13 skills" and expects a lightweight overlay; gets 43 skills on disk. Not a bug, but needs a note like: "Skills in this role (additional base skills not shown; run `sciagent status` after activating to see effective count)."

---

### FINDING 7 [NIT] `status --effective` mixes skills, agents, and commands with no separator

```
anndata
scanpy
...
skill-creator
cellranger-arc-multiome
...
docs-librarian      ← agent (no marker)
bio-interpreter
...
commit              ← command (no marker)
```

A machine reader can't distinguish which names are skills vs agents vs commands. The `--json` output correctly separates them into `.skills[]`, `.agents[]`, `.commands[]`. The `--effective` flat list is useful for simple grep but ambiguous for any structured consumer.

---

### FINDING 8 [ORPHAN] `.claude/` directory is completely removed on full deactivate

```bash
sciagent deactivate
ls .claude/   # → ls: cannot access '.claude/': No such file or directory
```

The entire `.claude/` directory (including subdirs) is deleted, not just cleared. If a user has manually placed files in `.claude/` (e.g., a handwritten `settings.json`), they're gone. The manifest at `.sciagent/manifest.json` is also removed (`.sciagent/` directory disappears). Same for `.agents/`.

**This is probably intentional** (clean teardown) but worth documenting explicitly. If any user files are expected to survive in `.claude/`, there's no protection.

---

### FINDING 9 [CONFUSING] Overlay-only activation has no `[overlay]` marker in status

```bash
sciagent activate scatac-regulatory        # solo, no base
sciagent status
# Stack:
#   1. scatac-regulatory roles/scatac-regulatory.yaml ...   ← no [overlay] tag
```

When scatac-regulatory is the only role, it shows without a tag. When it's in position 2 with a base, it shows `[overlay]`. The tag is structurally correct but could mislead a user who activates a role that was designed as an overlay into thinking they have a properly configured stack.

---

### FINDING 10 [NIT] `status --source <name>` exits 1 for valid non-source queries

```bash
sciagent status --source tobias-footprint-bindetect
# → not in active stack: tobias-footprint-bindetect (exit 1)
```

`tobias-footprint-bindetect` IS mounted (it's in `.claude/skills/` and `manifest.json`) but it arrived via `requires:` inheritance, not as an explicit role entry. `--source` has no path for inherited skills — it only queries `SKILL_ORDER` (role-declared) and `INJECTED_SKILLS`. The error message is actively wrong: the skill IS in the active stack, it's just inherited.

---

## Cross-checks: status vs disk

### base + scatac-regulatory active

| What | Status claims | Disk reality | Match? |
|------|--------------|--------------|--------|
| Total skills | 43 (text) / 40 (JSON) | 43 symlinks | Mismatch (3 inherited missing from JSON) |
| Agents | 7 | 7 symlinks | OK |
| Commands | 1 | 1 symlink | OK |
| Inherited skills in text | 3 (labeled) | 3 on disk | OK |
| Inherited skills in JSON | 0 | 3 on disk | MISMATCH |
| AGENTS.md hash | OK | — | OK |

### After full deactivate

| What | Status | Disk |
|------|--------|------|
| Skills | 0 | 0 |
| Agents | 0 | 0 |
| Commands | 0 | 0 |
| `.claude/` dir | gone | gone |
| `.agents/` dir | gone | gone |
| `.sciagent/` dir | gone | gone |
| AGENTS.md block | removed | — |

**Clean teardown confirmed: no orphaned symlinks, no dangling manifest, no block drift.**

---

## Depth cap enforcement

| Command | Expected | Actual | Correct? |
|---------|----------|--------|----------|
| `activate base scatac-regulatory pathway-signature` | Refuse (depth 3) | `Maximum stack depth is 2 (base + overlay)` exit 1 | YES |
| `activate pathway-signature` when base+scatac active | Unclear | Silently replaces stack with solo pathway-signature | CONFUSING (Finding 2) |
| `activate base pathway-signature` when base+scatac active | Replace overlay | Works cleanly | YES |

---

## Overlay switch path

**Clean path exists but is not obvious:**
```bash
# From base+scatac-regulatory to base+pathway-signature:
sciagent activate base pathway-signature   # WORKS — replaces stack entirely
```

**Friction:** You must repeat `base` in the command. There is no `sciagent activate --overlay pathway-signature` shorthand. The user's natural instinct ("just add pathway-signature") silently nukes the base (Finding 2).

---

## Idempotency

```bash
sciagent activate base scatac-regulatory   # first time
sciagent activate base scatac-regulatory   # second time
# → Activated stack: base scatac-regulatory  skills: 43  (same)
```
**Idempotency confirmed.** Double-activate produces identical state: same symlink count, same manifest, no duplicate entries.

---

## Role clarity: is scatac-regulatory obviously the right overlay for a chromatin person?

**`list role scatac-regulatory`:**
> scATAC-seq regulatory analysis — CREscendo, ChromVAR, TF footprinting, pycisTopic, SCENIC+

**From the YAML header:**
```yaml
# Use when
#   - scATAC-seq differential accessibility analysis
#   - TF motif activity scoring (ChromVAR)
#   - TF footprinting (TOBIAS, HINT-ATAC, Signac)
#   - Sub-peak CRE analysis (CREscendo)
#   - Building regulatory networks from ATAC data
```

**Verdict: YES.** For a chromatin/regulatory specialist the role description is unambiguous. The skill list in `list role` output directly names the tools (ChromVAR, CREscendo, pycisTopic, SCENIC+). No YAML reading was needed to confirm fit — the `list role` output is sufficient.

**One caveat:** The `list role` output shows 13 skills, but `activate base scatac-regulatory` gives 43 effective (Finding 6). This gap could cause confusion about what you're getting.

---

## Verdict on stacking UX

| Dimension | Score | Notes |
|-----------|-------|-------|
| Activation (2-arg) | 9/10 | Clean, fast, correct output. |
| Depth cap enforcement | 9/10 | Correct error, clear message. |
| Overlay switch (2-arg) | 8/10 | Works; requires repeating "base". No shorthand. |
| Overlay switch (1-arg) | 3/10 | **Silently nukes base.** No warning. (Finding 2) |
| Deactivation (full) | 10/10 | Clean, no orphans, block removed, dirs deleted. |
| Deactivation (partial) | 8/10 | Works. Deactivating base also removes overlay (announced). |
| Idempotency | 10/10 | Activate twice = same state. |
| `status` (text) | 9/10 | Rich, shadow-annotated, inherited skills labeled. |
| `status --effective` | 4/10 | **Wrong exit code (exit 1), omits inherited skills** (Findings 1, 3). |
| `status --json` | 7/10 | Correct structure, but omits inherited skills (Finding 3). |
| `status --source` | 6/10 | Wrong answer for inherited skills (Finding 10). |
| Role discoverability | 9/10 | `list role` description + Use-when YAML comments are clear. |
| Skill count coherence | 6/10 | `list role` shows 13; `status` shows 43 — same word "skills" (Finding 6). |

**Overall:** The core stacking lifecycle (activate → inspect → switch → teardown → re-activate) is solid. The mechanism works. Two structural bugs (`status --effective` exit code, inherited skills absent from programmatic outputs) and one serious UX trap (1-arg activate nukes the base) are the main issues to fix before this is production-ready for scripted use.
