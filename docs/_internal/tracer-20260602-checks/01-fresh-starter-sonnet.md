# Tracer-Bullet Report — Fresh Starter (RNA Project) — 2026-06-02

**Model:** claude-sonnet-4-6  
**Scratch project:** `/tmp/claude-788715489/tmp.npMWCYDf0a/my-scrna-project`  
**Toolkit:** `/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit`

---

## (a) Path Summary — Commands Run in Order

```bash
# Step 1 — Discovery
/data1/.../sciagent --help

# Step 2 — Scaffold project
/data1/.../sciagent new project /tmp/.../my-scrna-project \
    --type analysis --species mouse --genome mm10 \
    --title "My First scRNA-seq Analysis"

# Step 3 — Discover roles
/data1/.../sciagent list roles
/data1/.../sciagent list role base

# Step 4 — Activate role
cd /tmp/.../my-scrna-project
/data1/.../sciagent activate base
/data1/.../sciagent status
/data1/.../sciagent status --json

# Step 5 — Edge cases
/data1/.../sciagent activate          # no args
/data1/.../sciagent activat base      # typo
/data1/.../sciagent activate base     # second time (idempotency)
/data1/.../sciagent activate scrna-starter  # nonexistent role
/data1/.../sciagent new project       # no dir arg (from /tmp)
/data1/.../sciagent roster            # in README Quick start but not in CLI

# Extra probes
/data1/.../sciagent new project /tmp/.../my-tool --type software-tool
/data1/.../sciagent status --effective
/data1/.../sciagent status --source anndata
/data1/.../sciagent list deps tf-footprint-differential-analysis
/data1/.../sciagent validate
```

---

## (b) Numbered Findings

### 1. [BUG] `software-tool` type template is silently skipped — AGENTS.md, README.md, tool_config.yaml never rendered

**Evidence:**

```
$ sciagent new project /tmp/.../my-tool --type software-tool
wrote: .../docs/_internal/README.md
wrote: .../CLAUDE.md
wrote: .../gitignore
# — no AGENTS.md, no README.md, no tool_config.yaml
```

`lib/sciagent/new.sh:170` calls:
```bash
_render_tree "$proj_tpl/$type" "$dir" "$force"
```
where `$type` is `"software-tool"`. The template directory is
`templates/project/software/` (not `software-tool/`). `_render_tree` guards with:
```bash
[[ -d "$tpl_root" ]] || return 0   # lib/sciagent/new.sh:83
```
So when `templates/project/software-tool/` doesn't exist it silently returns 0.
`templates/project/software/` contains `AGENTS.md.template`, `README.md.template`,
and `tool_config.yaml.template` that are **never rendered**.

**Impact:** A software-tool project has no `AGENTS.md`, no project README, no
`tool_config.yaml` — the three most important files for that type.

**Fix:** Either rename `templates/project/software/` → `templates/project/software-tool/`
or add a type→dirname mapping in `_new_project`.

---

### 2. [BUG] `status --json` omits inherited (via `requires:`) skills; human text counts them

**Evidence:**

```
$ sciagent status
Skills (39 effective):
  ...
  tobias-footprint-bindetect       inherited via requires:
  hint-atac-differential-footprint inherited via requires:
  signac-footprint-visualization   inherited via requires:

$ sciagent status --json | python3 -c "import json,sys; d=json.load(sys.stdin); print(len(d['skills']))"
36
```

The JSON payload has 36 skills; the human output says 39. The three inherited skills
(`tobias-footprint-bindetect`, `hint-atac-differential-footprint`,
`signac-footprint-visualization`) appear in `.claude/skills/` symlinks and in
`AGENTS.md` managed block but are absent from `status --json`'s `skills[]` array.

**Impact:** Any consumer of the JSON manifest (CI, tooling, the Pi harness) sees
an incomplete skill list and will not know about the inherited skills that are
actually in `.claude/skills/`.

---

### 3. [BUG/CONFUSING] `AGENTS.md` "Active role" section stays stale after activation

**Template produces** (`templates/project/analysis/AGENTS.md.template:16-19`):
```markdown
## Active role

_No role activated yet. Run `sciagent activate <role>` to populate this section._
```

After `sciagent activate base` **this text is never updated**. The managed block
is appended below it but the "Active role" section still reads "No role activated yet"
every time an AI agent reads the file — directly contradicting the managed block.

**Impact:** An AI agent reading AGENTS.md top-to-bottom hits "no role activated"
first, then finds the managed block further down. Confusing and misleading.

**Fix:** Either delete the "Active role" section from the template (the managed
block is authoritative), or have `activate` splice-update it.

---

### 4. [BUG] `sciagent new project` (no dir arg) runs from CWD and sprays files into `/tmp` when used as shown in README

**README Quick start (`README.md:53-54`):**
```bash
# Bootstrap a new project directory with AGENTS.md, CLAUDE.md, context.md
sciagent new project
```

This is a real command that works (defaults to `.`). Running it from `/tmp`
scattered `AGENTS.md`, `CLAUDE.md`, `docs/`, `.gitignore` etc. across `/tmp` and
also spewed `find: Permission denied` noise to stdout for every protected systemd dir:

```
find: './systemd-private-208a44f2f955499591e7d204d3bf2d6b-...': Permission denied
...
Project scaffolded at . (type: analysis).
```

The `find` in `_seed_gitkeep` (`new.sh:74`) searches `"$dir"` which is `/tmp`
(all of `/tmp`). The errors go to stdout mixed with the normal output.

**Fix 1:** Change `find "$root"` to `find "$root" 2>/dev/null` (suppress permission errors).  
**Fix 2:** README Quick Start should show `sciagent new project <dir>` with an explicit path.

---

### 5. [ORPHAN] `README.md` Quick Start references `sciagent roster` — a verb that does not exist

**Evidence:**
```
$ sciagent roster
sciagent: unknown verb 'roster'
```

**README line 66:**
```bash
# List the active agents (add --json for a machine-readable manifest)
sciagent roster
```

No `roster` verb in the CLI (`--help`) nor in `lib/sciagent/`. The command for
listing effective agents is `sciagent status`. `sciagent list agents` lists all
available agents (not the active stack).

**Impact:** First-time reader runs the documented example and gets an error.

---

### 6. [CONFUSING] Skill count is inconsistent: `list roles` shows 36, `list role base` shows 36, `activate` reports 39

After activation:
```
Activated stack: base
  skills:   39
```
But `list roles` and `list role base` both display `36` because they only count
explicit role YAML entries. The 3 additional skills are pulled in transitively via
`tf-footprint-differential-analysis`'s `requires:` graph. This is not mentioned
anywhere in `list role base` output.

**Impact:** A newcomer reading `list role base` thinks they'll get 36 skills;
activation says 39; there's no explanation of where the extra 3 come from at the
`list` stage. The `list deps` command exists but is never mentioned in the role
description.

**Fix:** Add `(+N via requires)` or a `Requires expansion:` line to `list role` output.

---

### 7. [CONFUSING/NIT] `.agents/` directory is created silently with no documentation in new user flow

After `sciagent activate base`:
```
.agents/
├── agents/   (7 symlinks)
├── skills/   (39 symlinks)
└── commands/ (1 symlink)
```

This mirrors `.claude/` completely. The README mentions `.agents/` is for Pi harness
support (line 183) but the `status` output says:
```
Harness:  Claude Code detected (.claude/ present)   Pi: not detected (.pi/ absent)
```

A new user sees `.agents/` created in their project with no explanation of what it
is, and `status` tells them "Pi: not detected" — without ever explaining what Pi is
or whether they should care. "Pi" appears nowhere in `--help`.

Additionally, neither `.agents/` nor `.sciagent/` are in the `.gitignore` seed.
If the user commits the project the symlinks will break for anyone who doesn't
have the toolkit at the same absolute path.

---

### 8. [CONFUSING] Step 2 of Next Steps hard-codes `<scbio-docker>` placeholder — opaque to standalone toolkit users

After `sciagent new project`:
```
Next steps:
  1. cd /tmp/.../my-scrna-project
  2. Add container substrate (from scbio-docker):
       <scbio-docker>/scripts/init-container.sh /tmp/.../my-scrna-project --type analysis
  3. Open in VS Code → Reopen in Container
  4. sciagent activate base
```

A new user who found `SciAgent-toolkit` independently (not via `scbio-docker`)
has no idea what `<scbio-docker>` is or where to find it. The placeholder `<scbio-docker>`
is not resolved at print time. Steps 3 and 4 depend on step 2 being done, but
step 2 is optional for non-container workflows. There's no "skip if not using
Docker" note.

Source: `lib/sciagent/new.sh:238-239`.

---

### 9. [CONFUSING] README Quick Start comment says `context.md` but template creates `docs/_internal/scientific-context.md`

**README line 53:**
```bash
# Bootstrap a new project directory with AGENTS.md, CLAUDE.md, context.md
```

The actual file created is `docs/_internal/scientific-context.md`. There is no
`context.md` anywhere in the scaffolded tree. A user looking for `context.md`
won't find it.

---

### 10. [NIT] `list agents` output has inconsistent indentation across agent groups

```
$ sciagent list agents
bio-interpreter          ← no leading spaces (analysis-base group)
  captions
  ...
  insight-explorer
  architect              ← indented (architect group)
  bioinf
  ...
```

The `analysis-base/` agents are listed without leading spaces; the `architect/`
sub-agents are indented. Both groups are mixed into one flat list with no separator,
header, or group label. A newcomer cannot tell which agents are grouped, which are
for which role, or which to use.

---

### 11. [NIT] `templates/project/AGENTS.md.template` at project root is orphaned

File exists: `templates/project/AGENTS.md.template`

It is never picked up by `_render_tree` because `_new_project` only calls:
```bash
_render_tree "$proj_tpl/_common" ...
_render_tree "$proj_tpl/$type"   ...    # analysis/ or (missing) software-tool/
```

The root-level `AGENTS.md.template` is shadowed by the type-specific ones but
sits one directory up where `_render_tree` never walks.

---

### 12. [NIT] Double activation is silently idempotent but provides no feedback that it was already active

```
$ sciagent activate base   # second time
Activated stack: base
  skills:   39
  agents:   7
  commands: 1
```

No "already active, re-synced" or "replaced existing stack" message. For a new
user who accidentally runs activate twice, the identical output gives no signal
about whether anything changed.

---

### 13. [NIT] `docs/plan/` directory exists and is referenced in AGENTS.md documentation table but has no README

The documentation namespace table in `AGENTS.md` routes `public phased plan` to
`docs/plan/`. The directory is created with a `.gitkeep`, but there is no
`README.md` in it (unlike every other `docs/_internal/` subdirectory which has one).
An AI agent writing a plan has no naming convention guidance.

---

## (c) Newcomer Verdict

**Overall intuitiveness: 5/10**

**What works well:**
- `--help` is concise and shows correct usage; all verbs are listed.
- `new project` with explicit args is clean and fast; the scaffolded tree is
  well-documented and the `analysis_config.yaml` is thorough.
- `activate base` is one command and symlinks resolve correctly.
- `status` output is dense but readable; the managed block in AGENTS.md is clearly delimited.
- `validate` runs cleanly and the 2 warnings are expected/benign.
- `list deps` / `list dependents` are useful power-user commands once you know they exist.

**Where I got stuck:**
1. The Quick Start in README led me to run `sciagent roster` → immediate error.
2. `sciagent new project` (no dir) from `/tmp` polluted `/tmp` with project files
   and spewed permission errors — the recommended pattern in the docs.
3. After `activate base`, AGENTS.md still says "No role activated yet" — I had to
   scroll to the bottom to find the managed block.
4. `.agents/` appeared with no explanation; `status` reported "Pi: not detected"
   with no hint of what Pi is.
5. `software-tool` type produces an almost-empty project — no AGENTS.md, no tool
   config — because the template directory name is wrong.

**Blank spots a real new user needs:**
- Where to put the scientific question → `docs/_internal/scientific-context.md` is
  the answer but the README comment says `context.md`.
- What `.agents/` is and whether to commit it.
- How to add new analysis stages (mention of `analysis_config.yaml:stages` is in
  AGENTS.md but not in the Quick Start or README).
- What to do after step 4 (`sciagent activate base`) — there's no "now open Claude
  and start working" or pointer to any workflow docs.
