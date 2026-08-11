# Skill lifecycle

Skills have a *natural*, judgment-driven lifecycle. It is **not** enforced by
hooks or fail-closed checks — there is no PreToolUse/Stop hook and nothing in
`validate` fails on a skill's maturity. Deciding when a skill graduates or
retires is a human call made during practice.

```
in use  →  _attic  →  (delete if truly dead)
```

- **in use** — the skill lives at `skills/<name>/` and is mounted.
- **_attic** — retired. Moved to `skills/_attic/<name>/`, reference-only. See
  [`skills/_attic/README.md`](../skills/_attic/README.md).
- **delete** — when even the reference value is gone, remove the directory.

The directory a skill sits in *is* its lifecycle state. There is no
`status:` frontmatter field: it went with the rest of the `metadata:` block in
the Phase 4 frontmatter diet, because every field there was preloaded into
context for every skill in the catalog whether or not it was used. A skill that
is in flux, or on its way out, says so in the first line of its body — where the
reader is already looking.

## How it surfaces

`sciagent list skills` prints active skills, then a trailing attic section:

```
  scanpy
  scvi-mrvi
  ...

  Attic (retired, reference-only — not installed by any role):
    _attic/shinymultiome-uio-host
```

## The attic (`skills/_attic/`)

A skill under `skills/_attic/<name>/` is **reference-only**: not walked by
`validate`, not resolvable by `activate` (the walker looks at
`skills/<name>/`, not `skills/_attic/<name>/`), and not counted as an active
skill. The `skills/*/` glob does not recurse into it, and every skill walker
also skips underscore-prefixed dirs explicitly.

### Retiring a skill

1. Remove it from any `roles/*.yaml` that install it.
2. `git mv skills/<name> skills/_attic/<name>` (preserve history).
3. In the moved `SKILL.md`, add a one-line `> **Deprecated …**` banner at the
   top of the body.

### Reviving a skill

1. `git mv skills/_attic/<name> skills/<name>`.
2. Drop the banner.
3. Re-add it to the appropriate role(s).
