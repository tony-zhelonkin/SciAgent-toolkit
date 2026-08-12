# Architect command migration map

The source inventory contains **17 command files**, although the step brief says 16. The named groups resolve to 8 per-feature + 5 portfolio/meta (including status) + 4 audit commands. Their 2,927 `wc -l` lines plus the 47-line README reconcile to the stated 2,974 lines. This map covers all 17 command files.

`status.md` remains its own reference because the portfolio snapshot is a narrow frontmatter-and-grep route. The four `meta-*` specifications share `portfolio.md`; the four probationary audit specifications share `architecture-treemap-audit.md`.

Counts in the “after” column measure the preserved source section, excluding wrapper headings and `BEGIN/END SOURCE` markers. Grouped destination totals are 803 lines for `portfolio.md` and 597 lines for `architecture-treemap-audit.md`.

| Source command | Destination reference | Before | After |
|---|---|---:|---:|
| `commands/architect/map.md` | `skills/architecture-first-dev/references/map.md` | 103 | 103 |
| `commands/architect/review.md` | `skills/architecture-first-dev/references/review.md` | 231 | 231 |
| `commands/architect/synthesize.md` | `skills/architecture-first-dev/references/synthesize.md` | 102 | 102 |
| `commands/architect/design.md` | `skills/architecture-first-dev/references/design.md` | 378 | 378 |
| `commands/architect/architect.md` | `skills/architecture-first-dev/references/architect.md` | 129 | 129 |
| `commands/architect/plan.md` | `skills/architecture-first-dev/references/plan.md` | 140 | 140 |
| `commands/architect/verify.md` | `skills/architecture-first-dev/references/verify.md` | 291 | 292* |
| `commands/architect/diagram.md` | `skills/architecture-first-dev/references/diagram.md` | 112 | 112 |
| `commands/architect/meta-map.md` | `skills/architecture-first-dev/references/portfolio.md` | 125 | 125 |
| `commands/architect/meta-design.md` | `skills/architecture-first-dev/references/portfolio.md` | 256 | 256 |
| `commands/architect/meta-apply.md` | `skills/architecture-first-dev/references/portfolio.md` | 253 | 253 |
| `commands/architect/meta-plan.md` | `skills/architecture-first-dev/references/portfolio.md` | 154 | 154 |
| `commands/architect/status.md` | `skills/architecture-first-dev/references/status.md` | 71 | 71 |
| `commands/architect/components-extract.md` | `skills/architecture-first-dev/references/architecture-treemap-audit.md` | 109 | 109 |
| `commands/architect/audit-slice.md` | `skills/architecture-first-dev/references/architecture-treemap-audit.md` | 150 | 150 |
| `commands/architect/synthesize-audit.md` | `skills/architecture-first-dev/references/architecture-treemap-audit.md` | 228 | 228 |
| `commands/architect/architecture-treemap.md` | `skills/architecture-first-dev/references/architecture-treemap-audit.md` | 95 | 95 |

\* `verify.md` ended without a newline, so `wc -l` reported 291 for 292 textual lines. The migrated file contains the same 292 lines and a normalized final newline.

## Preservation checks

- Every per-feature and status reference matches its source after final-newline normalization.
- Every grouped section matches its source between the named source markers after the documented migration rewrites.
- Five internal source pointers were rewritten from `commands/design.md` or `commands/review.md` to their new `skills/architecture-first-dev/references/` paths; their line counts are unchanged.
- The audit slice-index template uses an equivalent HTML anchor for `slices/01_{slug}.md`, preventing the skill link-integrity check from treating an output placeholder as a bundled reference. Its rendered link and line count are unchanged.
- `$ARGUMENTS`, `$ARGUMENTS[N]`, and remainder/flag contracts remain inside the complete preserved specifications. The pre-migration repository-wide audit found indexed positional forms in 15 of 24 live commands.
- `commands/implement/implement.md` and `commands/decompose/decompose.md` remain in place and are routed from the skill as external command boundaries.
- `docs/workflows/architect/00-quickstart.md` and `docs/workflows/architect/01-architecture.md` remain canonical and unchanged.
- The audit blob records `status: probationary` where present and repeats the prune-review date `2026-11-20` at blob level.
