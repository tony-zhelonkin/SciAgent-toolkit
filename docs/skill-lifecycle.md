# Skill lifecycle

Skills have a *natural*, judgment-driven lifecycle. It is **not** enforced by
hooks or fail-closed checks — there is no PreToolUse/Stop hook, no `validate`
hard-fail on status. The toolkit only *surfaces* the state softly (in
`sciagent list skills`); deciding when a skill graduates or retires is a human
call made during practice.

```
experimental  →  stable  →  deprecated  →  _attic  →  (delete if truly dead)
```

- **experimental** — newly authored or in flux; usable but expect changes.
- **stable** — the default. Battle-tested; safe to rely on.
- **deprecated** — superseded or no longer maintained, but still installed
  somewhere or kept for one more cycle. A nudge to migrate off.
- **_attic** — retired. Moved to `skills/_attic/<name>/`, disconnected from all
  roles, reference-only. See [`skills/_attic/README.md`](../skills/_attic/README.md).
- **delete** — when even the reference value is gone, remove the directory.

## The `metadata.status:` field

Declared in a skill's `SKILL.md` frontmatter, under `metadata:`:

```yaml
metadata:
  scope: implementation
  status: experimental   # experimental | stable | deprecated
  ...
```

Allowed values: `experimental`, `stable`, `deprecated`. **Absent or empty means
`stable`** — untouched skills are stable by convention, so you only ever set the
field when a skill is *not* stable. `validate` does not hard-fail on this field
(or on an unexpected value); it is a soft convention.

## How it surfaces

`sciagent list skills` tags each non-stable skill with its status and prints a
trailing attic section:

```
  scvi-mrvi                  [experimental]
  some-old-skill             [deprecated]
  scanpy
  ...

  Attic (retired, reference-only — not installed by any role):
    _attic/shinymultiome-uio-host
```

Stable skills are shown plain (no tag) to keep the common case tidy.

## The attic (`skills/_attic/`)

A skill under `skills/_attic/<name>/` is **reference-only**: not walked by
`validate`, not resolvable by `activate`/`inject` (the resolver looks at
`skills/<name>/`, not `skills/_attic/<name>/`), and not counted as an active
skill. The `skills/*/` glob does not recurse into it, and every skill walker
also skips underscore-prefixed dirs explicitly.

### Retiring a skill

1. Remove it from any `roles/*.yaml` that install it.
2. `git mv skills/<name> skills/_attic/<name>` (preserve history).
3. In the moved `SKILL.md`, set `metadata.status: deprecated` and add a one-line
   `> **Deprecated …**` banner at the top of the body.

### Reviving a skill

1. `git mv skills/_attic/<name> skills/<name>`.
2. Set `metadata.status: stable` (or `experimental`) and drop the banner.
3. Re-add it to the appropriate role(s).
