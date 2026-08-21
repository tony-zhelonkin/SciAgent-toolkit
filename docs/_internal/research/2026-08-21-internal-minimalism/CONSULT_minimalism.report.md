# Recommendation

The structure is still one layer too deep. Keep directories for scope; remove directories that merely classify documents.

Recommended consumer shape:

```text
docs/_internal/
├── scientific-context.md                 # analysis only, populated when useful
├── _project/
│   ├── session.md
│   ├── <topic>.md
│   └── plans/<date-slug>/
│       ├── 00_INDEX.md
│       ├── session.md                    # optional mutable plan state
│       └── NN_<slug>.md
└── <stage-stem>/
    ├── session.md
    └── <topic>.md
```

Everything arrives with substantive content. A one-file plan is a project topic note; the plan directory is reserved for a genuinely phased plan.

## 1. Verdict on the eight locations

| Location | Verdict | Failure prevented |
|---|---|---|
| `<stage-stem>/session.md` | **Keep** | After interruption, the next worker otherwise has to reconstruct the stopping point, blockers, and next action from code and logs. |
| `<stage-stem>/reasoning/<topic>.md` | **Delete as written**; retain `<stage-stem>/<topic>.md` | Durable rationale prevents model choices, cutoffs, and rejected alternatives from being re-derived. The `reasoning/` directory itself prevents no failure. |
| `_project/session.md` | **Keep** | Pre-stage and cross-stage work otherwise gets falsely assigned to one stage or duplicated across several stages. |
| `_project/reasoning/<topic>.md` | **Delete as written**; retain `_project/<topic>.md` | Cross-stage rationale needs a home, but the category directory adds no information beyond the filename and scope directory. |
| `_project/plans/<date-slug>/` | **Keep for phased plans only** | A real multi-phase plan needs one authoritative phase map plus bounded implementation briefs. Flattening six phases into one growing file weakens handoff, review, and phase ownership. |
| `scientific-context.md` | **Keep, optional and populated** | Stable questions, cohorts, datasets, and hypotheses otherwise drift between stage sessions or are repeated throughout them. Never instantiate the current placeholder body. |
| `README.md` | **Delete** | No failure occurs: `AGENTS.md` routes writers and lint defines the mechanical contract. The field found generic READMEs participating in naming conflicts and empty scaffolds. |
| `<work-stem>/` in software repos | **Delete from the toolkit convention** | There is no existing software work key for lint to validate. Inventing one makes agents create vocabulary. Software repositories should use `_project/`; their issue/branch systems already partition work. |

The field evidence supports this reduction: 10 of 24 projects have no internal tree, empty category scaffolds misdirected every STING child, naming drift is extensive, and larger trees had lower signal density ([field findings](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/docs/_internal/research/2026-08-20-internal-skeleton/FINDINGS_field.md:6), [empty scaffolds](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/docs/_internal/research/2026-08-20-internal-skeleton/FINDINGS_field.md:31), [signal-to-noise](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/docs/_internal/research/2026-08-20-internal-skeleton/FINDINGS_field.md:106)).

## 2. Direct answers

### `reasoning/`

Delete it. Scope is useful structure; document type is not.

A stage directory containing `session.md`, `network-selection.md`, and `peak-floor.md` is unambiguous. The current extra level usually promises a collection and often delivers one file. The field’s conflicting reasoning spellings and README rules are evidence that this category has become a vocabulary problem ([findings](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/docs/_internal/research/2026-08-20-internal-skeleton/FINDINGS_field.md:43)).

### `_project/`

Keep it. `_project` is a scope key symmetrical with the stage key, rather than a document category. It earns the directory even when its first content is one `session.md`.

Root-level cross-stage notes would mix mutable working records with `scientific-context.md` and repository control files. A `00_` stem would claim preliminary-stage scope; cross-stage decisions remain applicable after preliminaries. `_project` states the actual scope.

### Plans

Keep a directory only for phased plans:

```text
_project/plans/YYYY-MM-DD-slug/
├── 00_INDEX.md
├── session.md
└── NN_slug.md
```

`session.md` wins over `00_STATE.md`. It is the same mutable continuity form everywhere. `00_INDEX.md` remains because it performs a different job: stable plan map, scope, phase order, and dependencies.

The 14782-DM record demonstrates that a 93-line state plus six phase briefs can carry real orchestration. The toolkit’s eight unread plan directories show the limiting rule: a directory does not make a plan current, reachable, or valuable. Therefore:

- A phased plan earns a directory.
- A single-file plan is `_project/<topic>.md`.
- No `00_STATE.md`.
- No `00_ARCHITECTURE.md`; fold architecture into `00_INDEX.md` or a project topic note.
- No plan-local `research/`, `reviews/`, or `_handoff/` categories.

### Two grammars

Today they are two conventions, and that is a smell. The analysis grammar has an observable key; `<work-stem>` does not. Current lint explicitly skips stem validation when stages are absent ([lint](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/lib/scio/lint.sh:760)).

Keep stage-keying for analysis. Use `_project/` alone for software. This removes the invented software vocabulary instead of pretending it is an equivalent grammar.

## 3. Lint recommendation

### Remove

- Remove the `handoffs|sessions|plans|reports|research` name blacklist. In analysis repositories, the unmatched-stage predicate already rejects them. In no-stage repositories it currently creates arbitrary false positives for possible work stems ([current predicate](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/lib/scio/lint.sh:775)).
- Remove the dependency on `reasoning/*.md`; require a nonempty immediate Markdown topic file or `session.md` ([current content test](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/lib/scio/lint.sh:806)).
- Remove the substring-based `*session*.md` counting rule. It can flag a valid phase such as `01_session-recovery.md`, while missing a lone legacy session file ([current counter](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/lib/scio/lint.sh:819)).
- Replace the extension blacklist with an allowlist.
- Remove wording that calls nested Git the recommended topology. Continue pruning `.git` so lint does not inspect repository machinery.

### Add

- Immediate child rule:
  - With stages: allow `_project/` and exact canonical stage stems.
  - Without stages: allow `_project/` only.
- Normalize `NN_topic_viz` to `NN_topic` when constructing recognized stage stems, so compute and visualization twins cannot create two memory scopes.
- Flag every empty directory below `_internal/`, excluding `.git`. Current lint catches `.gitkeep`, but misses a truly empty `reasoning/` beneath an otherwise valid stage.
- Permit only Markdown plus repository control files such as `.gitignore` and `.gitattributes`; flag other regular files and symlinks. The current blacklist misses the field’s PNG, PDF, PPTX, JSON, CSV, and similar landfill ([field payloads](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/docs/_internal/research/2026-08-20-internal-skeleton/FINDINGS_field.md:61)).
- Flag each legacy continuity filename individually, even when it is the only one: dated `session_*`, `handoff*`, `HANDOFFS.md`, `STATE.md`, and `00_STATE.md`.
- If `_project/plans/` exists:
  - plan directories match `YYYY-MM-DD-slug`;
  - each has nonempty `00_INDEX.md`;
  - phase files match `NN_<slug>.md`;
  - optional mutable state is exactly `session.md`;
  - reject other files.
- Make `internal-memory` a silent no-op when `--project-dir` is the toolkit checkout itself. It is a consumer grammar, and the toolkit has its own `toolkit` check.

### Keep loose

- An absent `_internal/` remains silent.
- `scientific-context.md` remains optional.
- `session.md` remains optional when a scope contains only durable topic notes.
- No predicate for correctness, completeness, currentness, or scientific quality.

Current misses include malformed plan trees because `_project/` passes after finding any nonempty descendant ([lint](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/lib/scio/lint.sh:783)). Current false findings include the toolkit’s own `handoff/`, `plans/`, and `research/`; I verified this by running the opt-in check against this checkout.

## 4. What to delete from ADR-D10

From [ADR-D10](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/docs/proposals/2026-08-11-offline-distribution/50_ADRs.md:319):

- Correct the Context’s present-tense claim that `docs-layout` warns when `_internal/` is absent. The implementation now treats absence as silence.
- In Decision 1:
  - delete both `reasoning/` layers;
  - delete the software `<work-stem>` convention;
  - retain stage scope and `_project`.
- Rewrite Decision 2 around optional phased plans; choose `session.md`, retain `00_INDEX.md`, and delete `00_STATE.md`.
- Reduce Decision 5 to one sentence: the opt-in tree lint is the enforcement surface. Delete the hook inventory and future pre-commit language; that question is settled.
- Delete Decision 6, the nested-repository-plus-parent-pointer commitment. The actual lint merely tolerates `.git`; it requires neither a nested repository nor `docs/internal-memory.md`. Both observed nested repositories lack remotes, so they still fail reachability ([field findings](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/docs/_internal/research/2026-08-20-internal-skeleton/FINDINGS_field.md:12)). Recording this as architecture while explicitly leaving creation and pushing unimplemented repeats the unsupported-claim defect.
- Consequently, delete “Git history is its archive” everywhere. It is false for the ordinary ignored, non-nested tree.
- Move or delete Decision 7’s ownership seam; it is valid architecture context but does not decide the memory structure.

The nested-repository question can return in a separate ADR when creation, reachability, and recovery are operationally defined together.

## 5. Template deletions

### Common README template

Delete [the entire template](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/templates/project/_common/docs/_internal/README.md.template:1). It duplicates `AGENTS.md`, already uses the vague `<work-stem>` wording, and does no work lint cannot do.

### Analysis `AGENTS.md`

In [the analysis template](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/templates/project/analysis/AGENTS.md.template:10):

- Make `scientific-context.md` conditional: “If present, read…”.
- Flatten both reasoning paths.
- Add the phased-plan route.
- Delete “Git history is its archive.”
- Delete the “never-grow-this-file” rule that routes excess instructions into private reasoning. Stable project instructions deserve an authoritative tracked home.

The placeholder [scientific-context template](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/templates/project/analysis/docs/_internal/scientific-context.md.template:1) should also be deleted or instantiated only after its placeholders are filled; otherwise “never scaffold empty” remains false in substance.

### Software `AGENTS.md`

From [the software template](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/templates/project/software/AGENTS.md.template:10):

- Delete every `<work-stem>` route.
- Route continuity and topic notes to `_project/`.
- Flatten `_project/reasoning/`.
- Delete the Git-history claim.
- Add phased plans only if the software template is intended to support them.

There are three additional live contradictions outside those templates:

- [CRAFT still routes plans to root `plans/`](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/craft.yaml:29).
- [The no-ephemeral hook still routes to `reasoning/`](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/templates/project/_common/.claude/hooks/no_ephemeral.sh.template:110).
- [`docs-layout` still recommends the retired `reports/` namespace](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/lib/scio/lint.sh:620).
- [The plan template itself names root `plans/` and root `reasoning/`](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/templates/plan/00_INDEX.md.template:4).

## Toolkit-only directories

- `handoff/`: **Delete as a category.** Keep one live `session.md` in a project-scope location. The directory currently mixes the live handoff with maps and a copied dev-env code bundle.
- `plans/`: **Keep for substantive multi-file plan bundles.** Single-file plans do not justify it.
- `research/`: **Keep in this toolkit.** Its dated evidence packets are directly cited by ADRs and consultations; removing the grouping would blur measured evidence with decisions. This is a toolkit audit function, outside the consumer grammar.

## Strongest case for “session only”

It is a serious case: most projects either have no tree or have trees whose signal declines as they grow; lint cannot determine whether a reasoning note is useful or stale; code, config, captions, and public documentation already hold reproducible outcomes. Every additional private form creates another retention channel.

I stop one step short of that conclusion. A mutable session cannot safely carry durable rejected alternatives, cohort exclusions, or cross-stage architecture: it either grows indefinitely or overwrites them. The minimum that survives the evidence is therefore:

- one mutable `session.md` per active scope;
- flat topic notes only when rationale has lasting value;
- one optional populated scientific frame;
- phased plan bundles only when there are actual phases.

That is the smallest structure in which every retained form prevents a distinct failure.

## Verification boundary

I verified ADR-D10, the current lint implementation and tests, the requested templates, the plan templates, architecture, handoff, and the toolkit’s present tree. I relied on the supplied 24-project field digest for fleet facts and did not re-derive them. The failure analysis and keep/delete judgments above are my inference from those verified sources.

No files were modified; pre-existing worktree changes were left untouched.