# Phase 01 — remove the three claims nothing backs

**Repo:** scio · **Blocks:** phase 02 (same file) · **Read first:** `00_INDEX.md` §1

## Goal

Delete three instructions that name something the mechanism does not deliver.
Nothing is added here; phase 02 adds the replacement. After this phase the
toolkit says less and every remaining claim is true.

## The three

**(a) `lib/scio/lint.sh:604`** — the `docs-layout` check emits
`"docs/_internal/ missing — run: scio link"` whenever a project has `docs/` but
no `docs/_internal/`. `link` never creates that directory: its only mention of
the path is `lib/scio/link.sh:278`, inside the `SCIO:GITIGNORE` block. So the
remediation cannot remediate, and the warning pressures exactly the empty
scaffolding the field survey shows is harmful (`FINDINGS_field.md` §3).

Delete the predicate and its comment block at `lint.sh:601-605`. Keep the
following `docs/_internal/ is NOT gitignored` predicate — that one is true and
enforceable. Keep the absent-subject guard above it.

**(b) `craft.yaml`** — the Reproducibility bullet calls `_scratch/` "the only
sanctioned throwaway zone". No project has one at its root; one has an empty
one (`FINDINGS_field.md` §6). Remove the claim. Do not replace it with another
location — phase 06 records why scratch is not memory at all.

**(c) `craft.yaml`** — the same bullet routes durable memory to
`docs/_internal/reasoning/` and handoffs to `docs/_internal/sessions/`. Both
are the flat namespaces phase 06 retires in favour of stage-keyed
`docs/_internal/<stage-stem>/`. Rewrite the memory sentence to name the
stage-keyed route and `session.md`, and nothing else.

**(d) `templates/project/_common/.claude/hooks/no_ephemeral.sh.template:110`**
— the message says *"Capture throwaway probes in `_scratch/` or a reasoning
trace (`docs/_internal/`)"*, pointing at a directory the toolkit gitignores and
one that does not exist. Point it at the stage-keyed route.

Changing the hook body changes its hash. `templates/PROVENANCE.sha1` must be
regenerated (`tools/gen-template-provenance.sh`) so a field copy of the old body
is recognised and refreshed rather than ceded — otherwise every already-bound
project warns once and writes a ceded marker.

## Verify

```bash
bash tests/run-all.sh                        # 62 → still 62, 0 failing
bin/scio lint --check toolkit --strict       # clean; CRAFT body still ≤ 25 lines
# (a): a project with docs/ but no docs/_internal/ must now be SILENT
# (b,c): read the rendered CRAFT diff deliberately — this change is INTENDED to
#        alter the block, so every consumer will see drift on re-pin. Confirm the
#        body got SHORTER.
SCIO_TOOLKIT=$PWD bash -c '. lib/scio/block.sh; . lib/scio/craft.sh; _craft_render_body' | wc -l
# (d): the regenerated PROVENANCE.sha1 must contain the OLD hash too
bash tests/test_template_provenance.sh
bash tests/test_hook_body_propagation.sh
bash tests/test_validate_docs_layout.sh
```

Add a test asserting a project with `docs/` and no `docs/_internal/` produces no
`docs-layout` output. That test is the point of (a).

## Do not

- Do not add the new `internal-memory` check — that is phase 02, same file.
- Do not create `docs/_internal/` anywhere, in templates or tests.
- Do not add net lines to `craft.yaml`'s rendered body.
