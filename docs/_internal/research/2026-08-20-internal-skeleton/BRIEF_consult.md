# Consultation — the stable skeleton for agent working memory

You are consulting on **architecture**, not implementing. Produce a written
recommendation. Do not modify any file.

The owner's summary of the problem, in their words: *"agents extensively use
scratchpad and tmp which then basically ships nowhere and the intent of any such
scratch and exploration stays unobservable, forgotten and lost."* And: *"It's
been constantly drifting for me and that's another thing that I would want to
stabilize across projects and have agents be aware of stable locations, no matter
the harness, or the specific agent running within the container, be it a Claude
Code, Codex or whatever."*

## 1. The concrete defect

`scio` is a three-verb toolkit (`link`, `craft`, `lint`) vendored as a git
submodule into ~24 analysis projects. It renders an always-on managed block into
each project's `AGENTS.md` containing this instruction:

> Durable project memory lives in tracked files alone — `AGENTS.md`, committed
> stages, decisions in `docs/_internal/reasoning/`, handoffs in
> `docs/_internal/sessions/`. Harness auto-memory is off and untracked harness
> state counts as no record.

The same toolkit writes `docs/_internal/` into the `SCIO:GITIGNORE` managed block
of every project's `.gitignore` (`lib/scio/link.sh:278`), and its `docs-layout`
lint check warns when `docs/_internal/` is **not** gitignored (fails under
`--strict`). The mandated memory location is therefore guaranteed untracked, and
by the instruction's own definition counts as no record.

Measured: 14761-DM holds 1,161 files / 13 MB under `docs/_internal/` (122 plans,
57 reasoning traces, 9 sessions, 3 reports), zero tracked. DC-nexus 4,169 files,
2 tracked. 14782-DM 2,881 files, 0 tracked. Across five sampled projects: 8,317
files, 2 tracked.

Note the ignore rule was not stupid — publishing candid agent reasoning to a
public GitHub repo is a genuine hazard. It chose safety and created amnesia. Any
recommendation must address both, not trade one for the other.

Also relevant: `.claude/hooks/no_ephemeral.sh` already detects agents writing
`/tmp/*.py|R|sh` probe scripts and tells them to *"capture throwaway probes in
`_scratch/` or a reasoning trace (`docs/_internal/`)"* — pointing at the location
the toolkit then ignores. Doctrine, detection, and destination all exist; the
ignore defeats all three.

## 2. What the owner has already rejected — do not re-propose it

An earlier suggestion was to anchor intent traceability on scio's **decision
gates** (`analysis_config.yaml` under `decisions.<stage>`, `status: APPROVED`)
and plan phases, requiring every reasoning trace to cite one. The owner rejected
this from lived experience:

> "never worked, they only created churn for agents to figure out what those
> knobs mean, and keeping such intent knobs within config files that are supposed
> to have their values readable by some scripts referencing and calling for
> parameter definitions in the knob — meh, I don't really like the idea. I've
> lived with this for a while now, and usually I only see agents churning on it,
> and the code instead of being clean code, having variables that are hard to
> define why are they in the first place. I would rather simplify things and not
> introduce such things."

**Their actual model of where intent lives:**

> "The intent rather lives, in my case at least, with the staged analysis
> artifact producing some viz I can essentially reason about, and based on that
> delve into branches of analysis, into different ideas of different viz ideas
> I'd want to explore etc."

So the anchor is **a stage and the figure it produced**. The human looks at a
figure, reasons, and branches. Any traceability scheme must hang off that natural
artifact, not off metadata a human must maintain by hand. Treat "the owner will
maintain an intent field" as a failed hypothesis.

**Hard constraint: `I really don't want this make more complex.`** A
recommendation that adds moving parts needs to justify each one against this. A
recommendation that *removes* things is more likely to be right. Deleting a
convention counts as progress.

## 3. Context you need about scio's design discipline

Violating these is a defect, not a trade-off:

- **Three verbs, one concern each.** `link` = availability (convergent), `craft`
  = standing text (idempotent), `lint` = enforcement (read-only). A fourth verb
  needs a very strong argument; prefer extending a concern.
- **Bash only on the core path.** No `yq`, `jq`, `python3`.
- **One directory symlink per harness category** (`.claude/{skills,agents,commands}`,
  `.agents/…` → toolkit dirs). Consequence: a directory symlink cannot filter or
  flatten, so **source shape must equal mount shape**.
- **The submodule pin IS the lock.** No lockfile.
- **`lint.sh` owns every convention it can measure.** A claim stated twice
  drifts, so each claim has exactly one home. If a convention exists only as
  prose with no predicate, that is a gap.
- **Always-on text is scarce.** The `craft.yaml` block is capped at 25 rendered
  lines (17 used) because it is preloaded in every session. Depth belongs in
  skills, loaded on demand.
- Managed blocks (`SCIO:CRAFT`, `SCIO:GITIGNORE`) are SHA1 drift-guarded, so
  scio can own a region of a file the human also edits.

**Fleet reality that constrains any design:** 24 real toolkit copies exist. 20
sit at one commit 52 behind `dev`; 10 of those report `OK clean` because they
track an upstream that is an ancestor of their own HEAD. `scbio-docker` shows pin
drift right now. **Any design whose correctness depends on pin discipline or
on humans bumping a reference will inherit this failure mode.** Weigh this
heavily.

## 4. Questions to answer

1. **Is `docs/_internal/` one thing or several?** A working hypothesis (accept,
   reject, or replace): three classes are currently conflated — **scratch** (the
   `/tmp` probe script; dies with the session; the artifact is near-worthless but
   the *intent sentence* behind it is the value), **reasoning** (why this cutoff,
   what was rejected; permanent), **continuity** (handoffs, session state; lives
   until superseded). They differ in retention, value, and publish policy. Does
   splitting them simplify or complicate?

2. **Retention.** 8,317 files across five projects is already most of the way to
   a landfill. What is kept, what is pruned, and what enforces the pruning? A
   memory system that ships everything is worse than none, because nothing in it
   can be found. Be concrete about the mechanism.

3. **Topology.** Evaluate, with worked mechanics rather than a preference:
   - (a) **Status quo** — untracked, dies with the container.
   - (b) **Track in the analysis repo.** Simplest. Cost: agent churn floods the
     human commit log; the publish hazard returns.
   - (c) **Sibling repo with a recorded backlink** — a separate repo, each of its
     commits naming the analysis SHA it describes in a trailer. No pin bump, so
     it cannot go stale, but no atomic guarantee either. **The owner explicitly
     could not envision how this works** — if you recommend it, show the concrete
     mechanics: where the directory lives, how an agent finds it, who commits and
     when, how a reviewer walks from a figure back to the reasoning, and what
     breaks when someone clones only one of the two repos.
   - (d) **Submodule of the analysis repo.** Atomic (the analysis commit pins the
     memory state) but makes durability depend on pin bumps — see §3.
   - (e) Something else.

   **A specific hypothesis to attack, derived from §2.** `docs/_internal/`
   becomes a git repo *in place* — a nested repo, not a submodule. The parent
   already gitignores `docs/_internal/`, so the nesting needs no new wiring: the
   ignore rule that causes today's amnesia is exactly what keeps the agent log
   out of the human log. No pin, no bump, nothing that can go stale. Backlink =
   a commit trailer naming the analysis SHA.

   Paired with it, an organising key that needs no maintained metadata: file
   memory **by stage number**, mirroring structure that already exists —
   `02_analysis/stages/30_grn.R` → `03_results/30_grn/figures/` →
   `docs/_internal/30_grn/`. The human already reasons in those numbers, so
   there is no field to fill and nothing for an agent to churn on, which is
   precisely what the rejected decision-gate scheme got wrong.

   Attack this. Known hazards: a nested repo is easy to forget to commit and
   easy to destroy with the project; `git status` in the parent reveals nothing
   about it; tooling (IDEs, `provc`, `module-vendor`) may not expect it; and
   stage-numbered directories may not fit memory that spans stages or precedes
   any stage. Say whether those are fatal, mitigable, or acceptable, and if
   fatal, what replaces it.

4. **Harness neutrality.** The skeleton must work for Claude Code, Codex, and
   whatever comes next. Today `.claude/`, `.agents/`, `.gemini/` all exist and
   `docs/_internal/` is the only harness-neutral location. Where does harness
   state stop and project memory begin? Note the observed failure: an agent wrote
   its scratch harness to
   `/tmp/claude-<pid>/<hashed-workspace>/<uuid>/scratchpad/` — a
   harness-chosen path outside the project, per-session, which no project-level
   convention can reach after the fact. Can that be redirected, or must the
   design accept it and capture only intent?

5. **The publish knob.** Cases, roughly in order of how soon each bites:
   personnel turnover (the record is the handoff); abandoned analysis (the record
   is the *only* output of value); published paper (ship intent alongside code);
   mid-flight collaboration (traces contain candid assessments of other people's
   work); audit (evidentiary value depends on being contemporaneous and
   tamper-evident); sensitive content (patient data, embargo — default closed).
   Note cases 5 and 6 pull opposite ways: audit wants an unbroken history,
   sensitivity wants selective withholding, and a published git history cannot be
   retroactively redacted. Does that argue for private-by-default plus an export
   step rather than a repo that flips public?

6. **What enforces it.** Which parts are `link` (materialize/bind), which
   `craft` (always-on standing text, ≤8 spare lines), which `lint` (measurable
   predicate), and which are simply documentation? Name anything that cannot be
   enforced and should therefore not be claimed.

## 5. Deliverable

1. **The skeleton** — a concrete directory tree with real names, that you would
   put in every one of the 24 projects. Show it as a tree. Justify each entry;
   if you cannot justify one, drop it.
2. **A topology recommendation** with the worked mechanics asked for in Q3.
3. **The enforcement map** — entry → mechanism (`link` / `craft` / `lint` /
   docs-only), per Q6.
4. **What to delete.** Existing conventions, directories, or claims that should
   go. Be specific and expect this section to be acted on.
5. **The migration.** 8,317 existing files across 24 projects, mid-flight. What
   happens to them, and what is the smallest first step that improves the
   situation without waiting on the full design?
6. **Where you disagree with this brief.** If the framing is wrong, say so and
   say why. If the three-class hypothesis is wrong, replace it.

Ground every recommendation in what you can verify by reading the repositories.
State plainly what you could not determine. Prefer the answer that removes
complexity; the owner has explicitly asked not to add any.
