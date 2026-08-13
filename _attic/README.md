# Skill attic (`_attic/`)

Retired skills live here, **reference-only**.

A skill under `_attic/<name>/` is:

- **not walked** by `sciagent validate` (and the skill-walking tests skip it),
- **not resolvable** by `sciagent activate` / `sciagent inject` (the resolver
  looks at `skills/<name>/`),
- **not counted** as an active/available skill in `sciagent list skills`
  (it appears only in the trailing "Attic" section, marked retired).

The attic is a soft, judgment-driven stage in the natural skill lifecycle
(`experimental → stable → deprecated → _attic → delete`). It is **not** enforced
by any hook or fail-closed check — see `docs/skill-lifecycle.md`.

## Reviving an attic'd skill

```bash
git mv _attic/<name> skills/<name>
```

Then edit `skills/<name>/SKILL.md`:

- set `metadata.status: stable` (or `experimental` if re-entering as a draft),
- drop the `> **Deprecated …**` banner at the top of the body,

and re-add the skill to whichever `roles/*.yaml` should install it.

## Truly dead?

If a skill is dead beyond reference value, delete its directory outright.
The attic is for skills worth keeping around for reference, not a graveyard.
