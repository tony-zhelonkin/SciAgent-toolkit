# Commit

Commit pending work as a clean, readable git history: small atomic commits, brief descriptive messages, no AI attribution.

Use this any time you'd otherwise type "please commit cleanly with brief messages and no AI co-author trailer." `$ARGUMENTS` is an optional free-text hint (e.g. a scope like `notebooks` or a directive like `squash into one`); if empty, infer scope from `git status`.

## Rules (non-negotiable)

1. **No AI attribution.** Do NOT add `Co-Authored-By: Claude ...`, "Generated with Claude Code" footers, or any reference to an AI/LLM/Anthropic in the commit message, body, or trailers. Commits authored solely by the human committer, not the LLM Provider.
2. **No `git add -A` / `git add .`.** Stage files explicitly by path. Never sweep in untracked files you didn't inspect — risk of leaking `.env`, credentials, large binaries, or unrelated WIP.
3. **Atomic commits.** One logical change per commit. If `git status` shows changes spanning multiple concerns (e.g. a doc edit + a script fix + a config bump), split them into separate commits in a sensible order. Don't bundle unrelated changes "to save a commit."
4. **Brief, descriptive subject lines.** Imperative mood, ≤72 chars, no trailing period. Match the existing repo style (check `git log --oneline -20`). Body only when the *why* isn't obvious from the diff — and keep it tight.
5. **Never amend or force-push** unless the user explicitly asks. Always create new commits. If a hook fails, fix the underlying issue and create a NEW commit — do not `--amend` or `--no-verify`.
6. **Don't push** unless asked.

## Procedure

1. Run in parallel: `git status`, `git diff`, `git diff --staged`, `git log --oneline -10`.
2. Scan for sensitive files (`.env`, `*credentials*`, `*.key`, `*.pem`, large binaries). If any are staged or in the change set, **stop and ask** before proceeding.
3. Group changes into atomic commits. Announce the planned grouping in one short block before staging — e.g.

   ```
   Plan: 3 commits
     1. docs: <scope> — <what>           [files: README.md, docs/foo.md]
     2. fix: <scope> — <what>            [files: src/bar.py]
     3. chore: <scope> — <what>          [files: config/baz.yaml]
   ```

   If `$ARGUMENTS` says "squash" or "one commit", collapse to a single commit instead.
4. For each commit: stage by explicit path, then commit with a HEREDOC message. Example:

   ```bash
   git add path/to/file1 path/to/file2
   git commit -m "$(cat <<'EOF'
   <type>: <subject ≤72 chars>

   <optional body — only if the why is non-obvious>
   EOF
   )"
   ```

   No `Co-Authored-By` trailer. No tool footer. No emoji unless the repo's existing log uses them.
5. After all commits: `git status` (verify clean) and `git log --oneline -<N>` (show the new history). Report the list of new commit SHAs and subjects to the user — nothing else.

## Failure modes to avoid

- Mixing a real fix with formatting/whitespace churn in the same commit — split them.
- Vague subjects like `update files`, `fix stuff`, `wip` — re-read the diff and write what actually changed.
- Re-committing a hook-rejected change with `--no-verify` — fix the lint/test/format issue instead.
- Staging a whole directory when only two files changed — be precise.