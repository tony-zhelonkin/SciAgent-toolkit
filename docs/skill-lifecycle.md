# Skill lifecycle

Skills have a *natural*, judgment-driven lifecycle. It is **not** enforced by
hooks or fail-closed checks — there is no PreToolUse/Stop hook and nothing in
`validate` fails on a skill's maturity. Deciding when a skill graduates or
retires is a human call made during practice.

```
in use  →  _attic  →  (delete if truly dead)
```

- **in use** — the skill lives at `skills/<name>/` and is mounted.
- **_attic** — retired. Moved to `_attic/<name>/`, reference-only. See
  [`_attic/README.md`](../_attic/README.md).
- **delete** — when even the reference value is gone, remove the directory.

The directory a skill sits in *is* its lifecycle state. There is no
`status:` frontmatter field: it went with the rest of the `metadata:` block in
the Phase 4 frontmatter diet, because every field there was preloaded into
context for every skill in the catalog whether or not it was used. A skill that
is in flux, or on its way out, says so in the first line of its body — where the
reader is already looking.

## The attic (`_attic/`)

A skill under `_attic/<name>/` is **reference-only**. Active-skill operations
walk `skills/<name>/`, so the top-level attic stays outside validation,
resolution, mounting, and active-skill counts.

### Retiring a skill

1. `git mv skills/<name> _attic/<name>` (preserve history).
2. In the moved `SKILL.md`, add a one-line `> **Deprecated …**` banner at the
   top of the body.

### Reviving a skill

1. `git mv _attic/<name> skills/<name>`.
2. Drop the banner.
