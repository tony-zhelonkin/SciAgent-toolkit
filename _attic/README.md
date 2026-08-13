# Skill attic (`_attic/`)

Retired skills live here, **reference-only**.

A skill under `_attic/<name>/` is:

- outside the `sciagent lint --check toolkit` walk,
- outside both harness catalog links,
- excluded from active-skill counts.

The attic is a soft, judgment-driven stage in the natural skill lifecycle
(`experimental → stable → deprecated → _attic → delete`). It is **not** enforced
by any hook or fail-closed check — see `docs/skill-lifecycle.md`.

## Reviving an attic'd skill

```bash
git mv _attic/<name> skills/<name>
```

Then drop the `> **Deprecated …**` banner at the top of
`skills/<name>/SKILL.md` and run `sciagent lint --check toolkit`.

## Truly dead?

If a skill is dead beyond reference value, delete its directory outright.
The attic is for skills worth keeping around for reference, not a graveyard.
