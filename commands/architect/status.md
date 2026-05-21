# Status

Produce a one-page portfolio-state summary from the artifacts under `docs/`. Dispatches the `status-reporter` sub-agent (Sonnet, tight-scoped, frontmatter + grep only); surfaces its markdown block to the chat — writes nothing.

Purpose: replace "open 10 files to figure out where we are" with a <30-second table. Use at the start of a session, after a context switch, or any time the architect feels "swirled."

## Phase 0: Parse arguments

- `$ARGUMENTS` — optional comma-separated list of feature slugs.
  - Empty → `scope = all` (sub-agent resolves by globbing `docs/*/map.md` and excludes `_meta`).
  - Explicit list → `scope = <the list>`.

No further argument parsing. Announce resolved scope in one line before dispatch.

## Phase 1: Dispatch `status-reporter`

Use the Agent tool with `subagent_type: "status-reporter"`.

Prompt template:

```
Produce a portfolio status report.

scope: {all | <comma-separated slugs>}
today: {YYYY-MM-DD}

Follow your Phase 1 (state extraction) + Phase 2 (cross-feature conflict scan) + Phase 3 (emit summary) exactly. Frontmatter + grep only; no full file reads; ≤4 Bash calls per concern; return the markdown block as your final message.
```

The sub-agent's tool-set is deliberately narrow (`Read, Grep, Glob, Bash`) and its model is Sonnet — the task is aggregation, not reasoning, and Sonnet is the cost-appropriate ceiling. Full orientation and output-format contract live in `agents/status-reporter.md`.

## Phase 2: Surface the sub-agent's output

The sub-agent returns a single markdown block (portfolio table + meta layer + conflicts + Recommended next step). Relay it to the user verbatim — no additional editorial, no reformatting. Do NOT re-read artifacts yourself to "verify" the sub-agent's numbers; that defeats the purpose of dispatch.

If the sub-agent's output is malformed (missing table, missing legend, missing Recommended next step), re-dispatch once with a reminder to follow the Phase-3 template exactly. If the second attempt also fails, surface the failure to the user and suggest they re-run `/status` after the next stable session — do not patch the output by hand.

## Phase 3: Diagnose `/status` ↔ `/verify` disagreement (cheap consistency check)

After the sub-agent's output is surfaced, the main agent runs ONE additional cheap check (no full reads):

For each feature whose row shows `shipped/total` < `total` (i.e., open phases), check whether `docs/{feature}/verify.md` exists and what its frontmatter `verdict:` says. If `verify.md` says `verdict: CLEAN` but the row shows open phases, the two artifacts disagree — emit a one-line warning beneath the table:

```
⚠ Disagreement: {feature} — /status sees {shipped}/{total} phases shipped, but verify.md (last run YYYY-MM-DD) says CLEAN. One of them is stale: either re-run `/implement {feature} --auto` to ship the open phases, or re-run `/verify {feature}` to update the verdict.
```

Conversely, if `shipped/total == total/total` but `verify.md` says `INCOMPLETE` or `DRIFT`, emit:

```
⚠ Disagreement: {feature} — all phases marked `shipped` but verify.md says {VERDICT}. Re-run `/verify {feature}` to refresh the verdict against current code.
```

This is the **diagnostic disagreement signal** — `/status` (cheap, frontmatter-driven) and `/verify` (expensive, code-driven) are now distinct sources of truth that should converge. When they don't, surface the gap immediately so the user can pick which one to refresh.

One Bash call total:

```bash
grep -H '^verdict:' docs/*/verify.md 2>/dev/null
```

Cross-reference against the sub-agent's per-feature shipped/total counts. Emit warnings only when there's a mismatch. Stay silent when they agree.

## Rules

1. **Read-only.** `/status` writes no files. Ever. If the user wants durable state, they commit the sub-agent's output to `docs/_meta/STATUS.md` manually.
2. **Fast.** Target <30 seconds wall-clock end-to-end (dispatch overhead + sub-agent work). Sub-agent target is <20s; the ~5-10s dispatch overhead is the cost of model-pinning to Sonnet and keeping the main-agent context clean.
3. **Sub-agent, not main-agent.** Dispatch `status-reporter`. The earlier "main-agent only" rule was reversed after observing that main-agent runs on Opus drifted from the frontmatter-only discipline and read full file bodies in loops — a dedicated Sonnet sub-agent with a tight prompt is both cheaper and more disciplined.
4. **No sub-sub-agent.** `status-reporter` is a leaf. It dispatches nothing.
5. **Relay verbatim.** The main agent does not re-narrate, re-format, or "improve" the sub-agent's table. Surface and stop.
6. **Conflict scan is best-effort.** The sub-agent marks suspected conflicts with `?`; the main agent does not upgrade or downgrade those flags.
