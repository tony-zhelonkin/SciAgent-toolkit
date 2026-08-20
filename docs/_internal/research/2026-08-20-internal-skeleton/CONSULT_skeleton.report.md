I read the brief, the field digest, all four survey reports where needed, and spot-checked the two nested repos plus the largest landfill examples. The important correction from `FINDINGS_field.md` is decisive: `14761-DM` and `14782-DM` already use nested repos in `docs/_internal/`; they have 195 and 708 tracked files respectively, but no remotes, so the failure is reachability and retention.

**1. The Skeleton**

```text
docs/
├── internal-memory.md              # tracked by parent repo; remote/policy pointer
└── _internal/                      # ignored by parent; private nested git repo
    ├── README.md                   # short local contract
    ├── .gitignore                  # Markdown-only memory, deny caches/binaries
    ├── _project/                   # repo-wide or pre-stage memory
    │   ├── reasoning/
    │   │   └── scope.md
    │   └── session.md
    └── 30_grn/                     # stage id copied from 02_analysis/stages/30_grn.R
        ├── reasoning/
        │   └── differential-network-selection.md
        └── session.md
```

Only `README.md`, `.gitignore`, and `docs/internal-memory.md` should be scaffolded. Stage directories are created when content exists. No `.gitkeep`, no empty `handoffs/`, no default `plans/`, no default `reports/`, no internal `_scratch/`.

`_project/` handles memory that precedes or spans stages. Normal memory is keyed by the stage stem: `02_analysis/stages/30_grn.R`, `03_results/30_grn/figures/`, and `docs/_internal/30_grn/`.

`session.md` is updated in place. Git history preserves old handoffs without leaving a stale dated pile in the working tree. Reasoning files are topic files updated or deleted when superseded; history remains the archive.

**2. Topology Recommendation**

Use `docs/_internal/` as an in-place nested private git repo, with a required remote and a tracked parent pointer at `docs/internal-memory.md`.

Mechanics:

- `link` creates `docs/_internal/`, initializes a nested repo if absent, writes its local `.gitignore`, and keeps the parent `.gitignore` rule for `docs/_internal/`.
- `docs/internal-memory.md` is tracked by the parent and records the private memory remote or an explicit `ask-owner` placeholder. This fixes the current no-remote failure without adding a submodule pin.
- Agents find memory with `git -C docs/_internal rev-parse --is-inside-work-tree`, then route by stage stem.
- Agents commit `docs/_internal` after updating `session.md` or reasoning, before ending a substantive session.
- Each memory commit includes `Analysis-Commit: <parent HEAD>` as a trailer. If the parent tree is dirty, the session file states that plainly; the backlink is still the last committed analysis SHA.
- A reviewer starting from `03_results/30_grn/figures/foo.png` looks at `docs/_internal/30_grn/`, then uses nested git history for chronology and trailers for the analysis SHA.
- A clone of only the parent repo gets code plus `docs/internal-memory.md`, but not private memory. A clone of only the memory repo gets stage-keyed reasoning plus trailers, but not executable code. That break is acceptable for private-by-default memory; publication should use an export step into tracked public docs.

I would not use a submodule. The fleet already shows pin drift, and submodules make durability depend on the same behavior that is failing. I would not track all memory in the analysis repo because candid reasoning, audit material, and sensitive content have different publish rules.

**3. Enforcement Map**

| Entry | Mechanism |
|---|---|
| Parent ignores `docs/_internal/` | `link` writes `SCIO:GITIGNORE`; `lint` checks `git check-ignore` |
| `docs/internal-memory.md` exists and is tracked | `link` materializes; `lint` checks parent `git ls-files` |
| `docs/_internal/` is a nested git repo | `link` initializes; `lint` verifies |
| Nested repo has remote/upstream | `lint`; fail under `--strict` |
| Nested repo has no uncommitted memory at session end | `lint`; agent workflow in `craft` |
| Memory commits carry `Analysis-Commit:` | `lint` can check recent commits or HEAD |
| Stage directories match `NN_slug` or `_project` | `lint` |
| `session.md` is the continuity file | `craft` names it; `lint` warns on new `handoffs/`, root handoff files, and dated session piles |
| Markdown-only memory | nested `.gitignore` from `link`; `lint` scans tracked and on-disk files for binaries, caches, venvs, checkpoints, parquet, bytecode, logs |
| Scratch stays ephemeral | `craft` and hook text; `lint` can forbid `docs/_internal/_scratch` |
| Semantic quality, honesty, redaction correctness | docs-only; cannot be truthfully enforced |

**4. What To Delete**

Delete the `handoffs/` convention. The survey found it mostly empty while real handoffs went to `sessions/`.

Delete `.gitkeep` scaffolding under internal memory. Empty directories misled agents in the STING child projects.

Delete `reports/` as a default internal sink. In the field it means everything from durable syntheses to arbiter caches.

Delete the root `_scratch/` claim from always-on text. It is absent almost everywhere; current lint already treats any `_scratch/` as ephemeral.

Delete the decision-gate requirement from always-on memory doctrine. The owner rejected it, and stage plus figure is the natural anchor.

Delete generic `_internal/` ignore rules. Use the exact `docs/_internal/` rule so unrelated `_internal` directories are not silently hidden.

Delete locally-authored project memory from harness directories. `.claude/`, `.agents/`, and `.gemini/` should hold harness bindings/settings, while project memory lives in the neutral nested repo.

**5. The Migration**

Do not bulk-add the 8,317 files.

First, preserve what already works: add remotes to the existing nested repos in `14761-DM` and `14782-DM`, push them, and commit their current dirty state or explicitly triage it. That immediately changes “dies with disk” into reachable memory.

For projects with `docs/_internal/` but no nested repo, initialize the nested repo in place, add the Markdown-only `.gitignore`, then commit only curated Markdown records. Leave venvs, checkpoints, caches, PNG previews, parquet, bytecode, and logs out. Extract one summary Markdown file when a bulky artifact contains useful intent.

For projects with no `docs/_internal/`, create only the minimal skeleton when the next agent needs memory. Absence is the dominant fleet state, so migration should not manufacture empty category trees.

For existing stale piles, use git history as the archive: keep only current `session.md` and current reasoning topic files in the working tree. Delete superseded files after their content is folded forward. The old versions remain in nested repo history.

Smallest toolkit step: change the `craft.yaml` memory sentence and the `no_ephemeral.sh` guidance so agents stop writing durable intent to a parent-ignored non-repo location. Point them to `docs/_internal/<stage>/session.md` in the nested repo, and point probes to local `_scratch/` with only intent captured in memory.

**6. Where I Disagree With The Brief**

The brief’s “zero tracked” framing is wrong after the field correction. The parent tracks zero, but two projects already track memory in nested repos. The missing piece is a remote plus enforcement.

The three-class hypothesis is close, but I would simplify it further: scratch is not memory; it is disposable work. The durable split is only `reasoning/` and `session.md`, both stage-keyed. Plans, reports, research caches, and scratch artifacts should not be scaffolded as first-class memory categories.

I also disagree with physically putting a full directory tree in every project. The field shows empty scaffold as active harm. Put the contract everywhere; create stage directories only when they carry content.