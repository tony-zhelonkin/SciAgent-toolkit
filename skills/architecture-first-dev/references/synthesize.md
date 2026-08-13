# Synthesize

Collapse multiple reviewer outputs into one consensus document. Only useful when ≥2 reviewers ran.

## Phase 0: Parse arguments

- `$ARGUMENTS[0]` — feature slug (required)

If missing, ask for it.

## Phase 1: Count reviewer files

List files matching `docs/{feature}/review/*.md` (exclude the `.history/` subtree).

- **0 files** → tell user to run `/review {feature}` first. Stop.
- **1 file** → tell user synthesis is unnecessary — reading one review directly is cheaper. Show the file path and offer to print its Summary section. Stop.
- **≥2 files** → proceed to Phase 1.5.

## Phase 1.5: Archive prior synthesis on re-run

If `docs/{feature}/synthesis.md` exists AND at least one `review/<name>.md` has `mtime > synthesis.md.mtime` (i.e., a newer review round landed after the last synthesis), this is a re-synthesis — archive the prior synthesis before overwriting:

1. Create `docs/{feature}/synthesis/.history/` if needed.
2. Determine round numbers. Let `H = max({i : docs/{feature}/synthesis/.history/round-i.md exists})` (or `0` if no history). The currently-on-disk synthesis is round `M = H + 1` (i.e. round-1 if no history, round-2 after the first re-synthesis, etc.). The synthesis we are about to produce is round `N = M + 1`.
3. `mv docs/{feature}/synthesis.md docs/{feature}/synthesis/.history/round-{M}.md` — do NOT overwrite a file already under `.history/`; abort and ask the user if collision detected (indicates a prior move was interrupted or round-numbering is wrong).
4. Print: `Archived prior synthesis to synthesis/.history/round-{M}.md; producing round-{N} synthesis from {K} reviewer files.`

If `synthesis.md` does not exist OR all reviewer files pre-date it (i.e., nothing new): this is first-run or no-op, skip archive.

Note: `/review --iterate` already archived the prior synthesis in its own Phase 2.5. This phase is a safety net for the case where the user re-ran `/review` without `--iterate` (manually overwrote reviews) and then re-runs `/synthesize`. Idempotent — if `synthesis.md` is already absent because `/review --iterate` archived it, Phase 1.5 is a no-op.

## Phase 2: Dispatch synth

Use the Agent tool with `subagent_type: "synth"`.

Prompt template:
```
Synthesize the reviews in docs/{slug}/review/ into docs/{slug}/synthesis.md.

Files to read:
  {list each review file — exclude the .history/ subtree; only primary files}

{If `docs/{slug}/synthesis/.history/round-i.md` exists for any i, inject — let M = max such i (post-archive, after Phase 1.5 has moved the prior synthesis), so the new synthesis is round N = M+1:}
You are producing an ITERATE synthesis (round {N}). Also read:
  docs/{slug}/synthesis/.history/round-{M}.md  (prior synthesis — preserve its `§ User resolutions` verbatim in the new synthesis unless a new reviewer explicitly supersedes a scribed decision)
  docs/{slug}/review/.history/round-{M}/*.md    (prior reviewer round — for provenance when the new round says "STOOD BY" a prior position)

Each new reviewer file contains a `## Round N iterate summary` section with UPDATED / STOOD BY / RETIRED tags — carry that shift-structure into the new synthesis so the reader sees what actually moved between rounds.

Output: docs/{slug}/synthesis.md with sections (in this order):
  - Headline (1-2 sentences, the single most important insight)
  - Convergent concerns (with reviewer attribution)
  - Divergent opinions (preserve disagreement, don't flatten)
  - Position shifts — ITERATE ROUNDS ONLY. Mirror the UPDATED / STOOD BY / RETIRED
    tags from each reviewer's `## Round N iterate summary`. Omit the section on
    round-1 syntheses.
  - Open questions
  - Prioritized recommendations (P0/P1/P2)
  - Reader map (index of who flagged what)
  - User resolutions — always-present. On round-1, initialise with the literal
    placeholder `<!-- empty -->`. On iterate runs, carry the prior synthesis's
    `## User resolutions` content forward verbatim (unless a new reviewer
    explicitly supersedes a scribed decision). The main agent appends binding
    chat decisions here between /synthesize and /design (scribe-on-latest,
    MADR-H3).

Rules:
  - Read only the review files. Do not re-explore the codebase.
  - Preserve disagreement — faithfully state each position when reviewers diverged.
  - Every claim traceable to a reviewer. No synthesis-original recommendations.
  - Keep total length ≤ combined reviewer length.
```

## Phase 3: Surface the synthesis

Once `synth` completes:

1. Read `docs/{feature}/synthesis.md`
2. Quote the Headline and P0 list to the user
3. Recommend next step:

```
Synthesis created: docs/{feature}/synthesis.md

Headline:
  {quote the headline}

P0 items:
  {quote the P0 list}

Recommended next step:
  /design {feature}                 — translate synthesis into design docs
  /review {feature} --as <name>     — add more reviewer lenses before designing
```

## Rules

1. **Do not invoke for <2 reviewers.** Tell user to read the single review directly.
2. **Synth reads only review files + (on iterate) the prior synthesis.** If you're tempted to have synth dig into the codebase, you've misunderstood the pipeline — escalate open questions to `/map` regeneration instead.
3. **Preserve disagreement.** Never let synth produce a flattened consensus when reviewers actually disagreed — that's how bad designs sneak past.
4. **Synthesize never destroys.** On re-synthesis (after a new review round landed), the prior `synthesis.md` is archived to `synthesis/.history/round-N.md` before being overwritten. The prior `§ User resolutions` section is preserved in the new synthesis (carry verbatim unless a new reviewer explicitly supersedes a scribed decision). See MADR-H5.
5. **Iterate synthesis reflects position-shift.** When synthesizing an iterate round, the new synthesis must surface the UPDATED / STOOD BY / RETIRED structure from the reviewer files — reconciliation is the point of the round; a round-2 synthesis that reads like a round-1 synthesis has failed.
