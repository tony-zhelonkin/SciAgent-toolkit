# Phase 07 — the codex-delegation seam holds; make the knowledge reachable

**Repo:** scio (one skill edit) · **Parallel-safe** · **Plus one owner decision**

## The finding: nothing is misplaced

An in-container agent hit *"Codex can't write files here — the container's
bwrap restriction is live (workspace-write refused both write and read)"* and
investigated it from scratch. Both repos had already answered it:

- **scbio-docker owns the cause and the fix.** `docs/ai-integration.md`
  §"Codex sandbox (bubblewrap) requirement": `bwrap` is baked into the image,
  Docker's default seccomp profile blocks the `clone`/`unshare` calls its
  unprivileged user namespace needs, so run with
  `--security-opt seccomp=unconfined`. The compose template already ships the
  line, commented, on both services
  (`templates/devcontainer/.devcontainer/docker-compose.yml.template:38-39, 89-90`).
- **scio owns the workflow and the fallback.** `skills/delegate-cli/SKILL.md`
  §"Linux devcontainer", *"Verified 2026-07-28 in the Meta-Aging containers.
  Four differences, each found by execution"* — states the same cause, gives the
  exact `--dangerously-bypass-approvals-and-sandbox` invocation, and supplies
  the compensating verification (`grep -cE '\bgit (add|commit|push|…)'` over the
  log plus `rev-parse HEAD`/`reflog` in **every** reachable repo, submodules
  included).
- **dev-env owns nothing here**, correctly.

So the ownership test needs no adjustment. `delegate-cli` was one of the **26
skills unreachable of 86** in that project — the skill documenting the failure
was invisible where the failure was being rediscovered. This is the same
delivery gap as everything else in this plan, in its most expensive form: an
agent spent a session re-deriving a written answer.

## Change (small)

`skills/delegate-cli/SKILL.md` restates scbio-docker's *cause* as well as the
workflow. One claim, one home: keep the invocation and the verification recipe
(workflow — scio's), and replace the cause sentence with a citation of
`scbio-docker/docs/ai-integration.md`, noting the commented compose line as the
alternative posture. Two sentences, not a rewrite.

Sanity-check the rest of §"Linux devcontainer" against reality while there —
it was verified 2026-07-28 and codex has since moved to 0.147.0. Do not
re-verify by execution as part of this phase; flag anything that looks stale.

## Verify

```bash
bin/scio lint --check toolkit --strict     # description ≤ 350 chars, name matches dir
bash tests/run-all.sh
bash tests/test_skill_frontmatter_valid.sh
```

## The owner decision this surfaces

Two valid postures, and the containers currently sit in the second:

**(a) Uncomment `security_opt: ["seccomp=unconfined"]`.** codex's sandbox works;
`-s workspace-write` behaves as documented; the container remains the outer
isolation boundary.

**(b) Leave it commented.** Keep the hardened seccomp profile and rely on
`--dangerously-bypass-approvals-and-sandbox` plus prompt-level prohibitions and
post-hoc grep verification.

**Recommendation: (a), for the delegation workflow.** Posture (b) makes safety
depend on an instruction in a prompt and a grep afterwards — which is precisely
the "instruction the mechanism does not enforce" shape this entire plan exists
to remove, applied to the one case where the consequence is an agent writing to
a real repo. (a) restores a mechanism. Against it: `seccomp=unconfined` is a
genuine loosening on a shared lab machine, so it is the owner's call, not an
implementer's.

If (a) is chosen it is a **scbio-docker** change (uncomment in the template,
note it in `docs/ai-integration.md` as the default rather than the option) and
requires container recreation via `provc` to take effect — out of scope for this
plan either way.
