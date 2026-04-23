# Map

Produce a single cartography artifact (`docs/{feature}/map.md`) that every downstream stage reads instead of re-exploring the codebase.

This is the token-efficiency foundation of the `architect` role: 60–80% of agent tokens typically go to orientation — we pay that cost once, in this command, and amortize it across all downstream stages.

## Phase 0: Parse arguments

- `$ARGUMENTS[0]` — feature slug OR path to an existing feature directory
- `$ARGUMENTS[1+]` — optional feature description / issue link

If the first argument is missing, ask:
```
Please provide:
1. Feature slug (e.g., "umap") OR path to feature dir (e.g., "docs/umap")
2. Optional: one-paragraph description of the feature / problem to map
```

### Phase 0.1: Slug resolution

If `$ARGUMENTS[0]` contains `/` or starts with `.` or `/`, treat it as a **path** and derive the slug:
- Absolute path ending in `docs/umap` → slug = `umap`, feature dir = `docs/umap`
- Relative path `docs/umap/` → slug = `umap`, feature dir = `docs/umap`
- If the path exists but is outside the project `docs/` tree, warn the user and ask whether to proceed with the feature dir at that exact path or to use `docs/{basename}` instead.

If `$ARGUMENTS[0]` is a plain slug (no `/`), use it directly: slug = the argument, feature dir = `docs/{slug}`.

**Always announce the resolution before proceeding**:
```
Resolved:
  slug = {slug}
  feature dir = {path}
  map file = {path}/map.md
```

If the resolution could have been ambiguous (path given but basename doesn't match a known convention; or the feature dir doesn't exist yet), **WAIT for user confirmation** before Phase 1.

## Phase 1: Check for existing map

Always check if `{feature_dir}/map.md` already exists — this phase always runs, just with different outcomes:

- **Exists:** ask the user to choose — (1) use existing, (2) regenerate/overwrite, (3) append updates. WAIT for choice.
- **Does not exist:** state `No existing map.md — proceeding to mapper dispatch.` and continue to Phase 2 without prompting.

Do NOT describe this phase as "skipped" when the file doesn't exist — the check was performed; there was simply nothing to prompt about.

## Phase 2: Dispatch mapper

Use the Agent tool with `subagent_type: "mapper"`.

Prompt template:
```
Map the feature `{slug}`. {optional description from $ARGUMENTS[1+]}

Output path: docs/{slug}/map.md

Cover:
1. Blast radius — files touched, tests, configs, explicit out-of-scope
2. Entry points and data flow with file:line references and a Mermaid flowchart
3. Existing patterns to reuse
4. Open questions — anything you could not resolve by reading alone

End the file with an always-present section:

```markdown
## User resolutions

<!-- Binding chat decisions the user makes between /map and /review|/design
     are appended here by the main agent (scribe-on-latest). Empty on first
     draft. Every downstream command reads this section as authoritative. -->
```

Facts only. No opinions, no suggestions, no critiques.
```

The mapper agent has tools: Read, Grep, Glob, Bash, Write. It will produce exactly one file.

## Phase 3: Surface the map

After mapper completes:

1. Read `docs/{feature}/map.md`
2. Surface the "Open questions" section to the user — these often need human/domain input before the next stage
3. Recommend the next command:

```
Map created: docs/{feature}/map.md ({N} lines, {M} file:line references)

Open questions for you:
{quote the open questions from the map}

Recommended next steps:
  /review {feature} --as all           — full 6-reviewer panel (large or cross-cutting features)
  /review {feature} --as bioinf,ml     — targeted review
  /design {feature}                    — skip review if the feature is small/unambiguous
```

## Rules

1. **Map is facts, not design.** If the user asks you to make the map suggest fixes, refuse — that's `/review`'s job.
2. **One file.** Mapper writes `map.md` only.
3. **Every claim in the map cites file:line.** Verify this on mapper's output before surfacing to user.
4. **If mapper's output has no Open questions section,** ask mapper to add it — every feature has some irreducible uncertainty; papering over it creates downstream surprises.
