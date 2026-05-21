# Review

Run a composable expert panel over a feature's map + design. Each reviewer writes one file to `docs/{feature}/review/<reviewer>.md`. Reviewers run in parallel.

## Usage

```
/review <feature>                                    — interactive: asks which reviewers
/review <feature> --as all                           — all 6 reviewers, parallel-independent (round-1 regime)
/review <feature> --as all --but stat,divergent      — all 6 minus listed
/review <feature> --as bioinf,ml                     — specific subset
/review <feature> --as wetlab                        — single reviewer

/review <feature> --iterate [--as <spec>]            — cross-informed iterate round
                                                       (reconciliation regime; see MADR-H5)
```

### Two regimes

Round-1 (default) and iterate rounds serve different epistemic purposes:

| Regime | Invocation | Reviewer context | Purpose |
|---|---|---|---|
| **Round 1 — parallel-independent** | `/review <feature> --as ...` (no `--iterate`) | Own lens + `map.md` + scope + any `design/` drafts. **Blind to other reviewers.** | Surface divergence. Each lens gives an honest independent read. |
| **Iterate — cross-informed** | `/review <feature> --iterate --as ...` | All of the above PLUS every other reviewer's prior-round output PLUS `synthesis.md` (incl. `§ User resolutions`). | Reconcile after divergence has been seen. Each reviewer *refines or explicitly stands by* its prior position. |

Round-1 is always parallel-independent; iterate is always cross-informed. Do not mix.

## Reviewer roster (the `all` set)

| Name        | Lens                                                             |
|-------------|------------------------------------------------------------------|
| `bioinf`    | Computational biology / bioinformatics methods correctness       |
| `wetlab`    | Bench utility, experimental prioritisation, tractability         |
| `graphic`   | Tufte information graphics: data-ink, perceptual encoding, UX    |
| `stat`      | Statistical inference, multiple testing, power, robustness       |
| `divergent` | Contrarian: blind spots, hidden assumptions, what-breaks         |
| `ml`        | Manifolds, embeddings, distance metrics, dim-reduction choice    |

## Phase 0: Parse arguments

- `$ARGUMENTS[0]` — feature slug (required)
- Flag `--as <spec>` — `all` OR a comma-separated list of reviewer short-names
- Flag `--but <csv>` — ONLY valid with `--as all`. Comma-separated reviewers to exclude.
- Flag `--iterate` — triggers the cross-informed iterate regime. Requires at least one prior round of reviews (`review/*.md`) on disk; rejects otherwise. Default reviewer set when `--iterate` is passed without `--as` is the set of reviewers from the most-recent round (all of them, iterating together is the common case).

If slug is missing, ask for it.

If `--as` is missing AND `--iterate` is not set, show the roster table and ask which reviewers to run (default suggestion: `all` for features with a design drafted, `bioinf,ml,divergent` for a feature still at map stage).

If `--as` is missing AND `--iterate` IS set, resolve the reviewer set from the filenames under `docs/{feature}/review/*.md` (the round-about-to-be-archived). This is the common case: re-run the same reviewers that produced the prior round so every lens reconciles. Do not prompt — the prior round's filenames are the authoritative default.

If `--but` is used without `--as all`, reject:
```
--but is only valid with --as all. To run a specific subset, use --as a,b,c instead.
```

If `--iterate` is passed but no prior reviews exist, reject:
```
--iterate requires at least one prior round of reviews under docs/{feature}/review/. Run /review {feature} --as <spec> first (round-1 is parallel-independent; iterate reconciles).
```

## Phase 1: Verify prerequisites

Check `docs/{feature}/map.md` exists. If not:
```
docs/{feature}/map.md not found. Run /map {feature} first — reviewers read the map rather than re-exploring code.
```
Then stop.

## Phase 2: Resolve reviewer set

Given `--as` + `--but`, compute the final set:
- `--as all` → `{bioinf, wetlab, graphic, stat, divergent, ml}`
- `--but stat,divergent` on `all` → `{bioinf, wetlab, graphic, ml}`
- `--as bioinf,ml` → `{bioinf, ml}`

Validate each name exists in the roster. If any name is unknown, error with the valid names listed.

**On `--iterate`, classify the resolved set into two buckets** by checking which have prior-round files under `docs/{feature}/review/*.md`:

- **Returning reviewers** — have a prior-round file. They get the iterate dispatch prompt (cross-informed, position-shift tagging required).
- **Fresh reviewers** — named in `--as` but no prior-round file (e.g., adding `divergent` in round-2 of a feature whose round-1 was `--as bioinf,ml`). They get the round-1 dispatch prompt instead (blind to siblings, no position-shift tagging — nothing to reconcile against). Surface this to the user before dispatch:

```
Round {N} dispatching with mixed regime:
  Returning (reconcile): {list}
  Fresh (round-1 regime): {list}
Fresh reviewers skip the iterate position-shift section; they have no prior position to refine.
```

If every reviewer is fresh (zero overlap with prior round), reject — this is a round-1 invocation wearing an `--iterate` hat, and the archive/prompt-regime distinction adds nothing:
```
No reviewer in --as {list} has a prior-round file. Iterate has nothing to reconcile. Run without --iterate to produce a fresh round-1 for these reviewers, or --as a set that overlaps with the prior round.
```

Print:
```
Running {N} reviewers in parallel: {list}
Each will read docs/{feature}/map.md + scope.md + any existing design drafts, and write docs/{feature}/review/<name>.md.
```

## Phase 2.5: Archive prior round (`--iterate` only)

If `--iterate` was passed:

1. Determine round numbers. Let `H = max({i : docs/{feature}/review/.history/round-i/ exists})` (or `0` if no history). The currently-on-disk round is `M = H + 1` (i.e. round-1 if no history, round-2 after the first iterate, etc.). The round we are about to dispatch is `N = M + 1`.
2. Create `docs/{feature}/review/.history/round-{M}/`.
3. Move every `docs/{feature}/review/<name>.md` that exists to `docs/{feature}/review/.history/round-{M}/<name>.md`. Use `mv` via Bash; do NOT overwrite any file that already exists under `.history/` (indicates a prior move was interrupted — abort and ask the user).
4. Archive the prior synthesis: if `docs/{feature}/synthesis.md` exists, create `docs/{feature}/synthesis/.history/` if needed, then `mv docs/{feature}/synthesis.md docs/{feature}/synthesis/.history/round-{M}.md`. Do NOT overwrite an existing `.history/round-{M}.md` — abort and ask the user if collision detected (same discipline as step 3). The next `/synthesize` run will produce the new `synthesis.md`.

After archive, print:
```
Archived round {M} to docs/{feature}/review/.history/round-{M}/ ({K} reviewer files).
Archived prior synthesis to docs/{feature}/synthesis/.history/round-{M}.md.
Round {N} reviewers dispatching with cross-informed context.
```

If `--iterate` was NOT passed: skip this phase. Round-1 writes directly to `review/<name>.md`; there's no prior round to archive.

## Phase 3: Dispatch reviewers in parallel

**CRITICAL:** Use a SINGLE assistant message with MULTIPLE Agent tool calls — this is what makes them run in parallel. Do not dispatch sequentially.

For each reviewer in the resolved set, one Agent call with `subagent_type: "<reviewer-name>"`. The prompt differs by regime:

- **No `--iterate`** → every reviewer gets the round-1 prompt.
- **With `--iterate`** → returning reviewers get the iterate prompt; fresh reviewers (no prior-round file; see Phase 2 classification) get the round-1 prompt so they don't hallucinate a prior position.

### Round-1 prompt (parallel-independent — default regime, or fresh reviewers on an iterate run)

```
Review the feature `{slug}` from your lens ({lens summary}).

Read in this order:
1. docs/{slug}/map.md
2. docs/{slug}/scope.md (if present)
3. docs/{slug}/design/*.md (if present — this is a re-review post-design)
4. docs/{slug}/synthesis.md (if present — prior synthesis you're reacting to)

Do NOT re-grep the codebase. If the map lacks context you need, add an "Open questions for mapper" section instead of investigating yourself.

You are in ROUND-1 regime: you are deliberately blind to what other reviewers are saying in parallel. Give your honest independent read from your lens. Divergence from sibling reviewers is a feature, not a bug — do not pre-reconcile.

Write exactly one file: docs/{slug}/review/{your-name}.md, matching the output format in your agent definition.
```

### Iterate prompt (cross-informed, `--iterate`)

```
You are reviewing the feature `{slug}` in the ITERATE regime (round {N}), from your lens ({lens summary}). You have a prior-round output from your own lens; you are REFINING it in light of what other reviewers said and what the user has locked as a design direction.

Read in this order:
1. docs/{slug}/map.md (current)
2. docs/{slug}/synthesis.md (if present — may have been re-synthesized post round-1; also read `§ User resolutions` for the user's locked decisions)
3. docs/{slug}/synthesis/.history/round-{N-1}.md (prior synthesis you're refining against)
4. docs/{slug}/review/.history/round-{N-1}/{your-name}.md (your OWN prior-round output)
5. docs/{slug}/review/.history/round-{N-1}/*.md (every OTHER reviewer's prior-round output — you are explicitly cross-informed now)
6. docs/{slug}/design/*.md (if present)

Do NOT re-grep the codebase.

Regime discipline — you are REFINING, not re-doing:
  - For each position in your prior round, EXPLICITLY state one of:
      (a) "UPDATED" — the sibling reviews or locked direction changed my mind. Say what changed and why.
      (b) "STOOD BY" — the sibling reviews did not move me. Say why you stood by it.
      (c) "RETIRED" — my prior position turned out to be based on a wrong premise (cite which sibling or which fact).
  - Do NOT silently re-write your round-1 as if you'd always held the round-2 view. Position-shift history is what makes iterate a reconciliation, not a premature-consensus machine.
  - Divergence is still legitimate if the siblings didn't convince you. Stand by explicitly rather than drift.

The user's locked direction (from synthesis.md § User resolutions) is binding; work under it as a constraint, not a position to re-debate.

Write exactly one file: docs/{slug}/review/{your-name}.md, matching your agent's output format. Begin the file with a `## Round {N} iterate summary` section that lists your UPDATED / STOOD-BY / RETIRED tags before the lens-specific body — the user and `synth` should see the position-shift at a glance.
```

For the iterate prompt, **inject the raw bodies** of sibling reviewers' prior-round files into the dispatch prompt (or ensure the reviewer can `Read()` them via the paths above — the agent has Read access). Injecting verbatim is cheaper on token count and avoids agent-miss on filename conventions; use `Read()` only if the combined body would exceed the dispatch prompt size budget.

## Phase 4: Summarize and route

After all parallel reviewers complete:

1. List the files written
2. Read the Summary section of each review
3. Surface a short cross-reviewer headline

```
Reviews complete ({N} files):
  - docs/{feature}/review/bioinf.md    ({M1} lines)
  - docs/{feature}/review/wetlab.md    ({M2} lines)
  ...

Per-reviewer headlines:
  bioinf:    {their summary's first sentence}
  wetlab:    {their summary's first sentence}
  ...

Recommended next steps:
  /synthesize {feature}           — only if ≥2 reviewers ran (otherwise read the single review directly)
  /design {feature}               — skip synthesis if reviews agree and you're ready
```

### Phase 4 addendum — iterate rounds

For iterate rounds, the summary block above is replaced (not supplemented) by a position-shift overview read from each reviewer's `## Round {N} iterate summary` section:

```
Iterate round {N} complete.

Position shifts (from each reviewer's round-{N} iterate summary):
  bioinf:    {C_updated} updated / {C_stood} stood by / {C_retired} retired
  wetlab:    ...
  ...

Prior round archived at docs/{feature}/review/.history/round-{N-1}/.
Prior synthesis archived at docs/{feature}/synthesis/.history/round-{N-1}.md.

Recommended next step: /synthesize {feature}
  (synth will read the new round-{N} reviews + the archived round-{N-1} synthesis to produce the new consensus)
```

## Rules

1. **Parallel dispatch is mandatory.** One message, multiple Agent calls. If you dispatch sequentially you've violated the design.
2. **Reviewers read the map — they do not re-explore.** Enforce this in their prompt.
3. **One file per reviewer.** If a reviewer writes elsewhere, flag it and ask user.
4. **Synthesis is conditional.** With only one reviewer, tell the user to just read the one file — don't auto-invoke `/synthesize`.
5. **Never invent a reviewer name.** The six names are canonical. Adding more requires a new agent file + role YAML update.
6. **Round-1 is always parallel-independent; iterate is always cross-informed.** Do not mix the regimes. Round-1 without `--iterate` forbids injecting sibling reviews (divergence-surfacing would be defeated); iterate with `--iterate` requires injecting them (reconciliation is the whole point). See MADR-H5 for the rationale.
7. **Never overwrite a review file.** On `--iterate`, the prior round is archived to `review/.history/round-N/<name>.md` *before* the new round is dispatched. If a file under `.history/` would be overwritten, abort and ask the user — this indicates a prior move was interrupted or the round-numbering is wrong.
8. **Iterate requires a prior round.** There is no "round-1 iterate." Reject `--iterate` if `review/*.md` is empty; the user should run `/review <feature> --as ...` first.
9. **Position-shift must be explicit.** Iterate reviewers are required to tag each prior position as UPDATED / STOOD BY / RETIRED. Silently rewriting a round-1 position as if it was always held defeats the reconciliation purpose. The `## Round {N} iterate summary` section at the top of each iterate-round review makes this auditable.
