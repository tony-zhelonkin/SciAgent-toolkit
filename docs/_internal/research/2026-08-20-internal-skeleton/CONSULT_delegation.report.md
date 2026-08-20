1. Verification and corrections

- §2(a) is correct. The skill already requires exact numeric targets with deviations reported, forbids weakening/deleting tests, and makes independent review mandatory ([SKILL.md](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/skills/delegate-cli/SKILL.md:152)). It also records Anton’s `!`-prefix role ([SKILL.md](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/skills/delegate-cli/SKILL.md:211)).
- §2(b) is correct: the codex prompt, final response, run log, and agy prompt/log examples all use `/tmp` ([SKILL.md](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/skills/delegate-cli/SKILL.md:100)). `/tmp` is appropriate for runtime prompt renderings, logs, byte counts, and final messages. It is inappropriate only when the prompt is the sole record of durable decisions.
- §2(c) is correct about reachability: Meta-Aging currently has 86 skill directories in its pinned toolkit, zero mounted skills, no `.agents/`, and no `.claude/skills/`. The umbrella therefore points to a skill it cannot load.
- The supposed catalog gap is no longer real. Both the canonical 416-line skill and Meta-Aging’s pinned 419-line copy already document `agy --model` labels/slugs, rejection of `--effort`, and the `-p`/bare-`--print` parsing trap ([SKILL.md](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/skills/delegate-cli/SKILL.md:70)). The line-count difference is frontmatter/“See also” material, not those operational facts.
- Future generic field findings should be recorded locally as evidence, then promoted into `delegate-cli`. The project copy should be removed once the catalog version carries the finding.

2. Position on the artifact split

| Artifact | Verdict |
|---|---|
| `probe.sh` | Yes: toolkit asset. It is reusable capability machinery. Parameterize the working directory/model and retain the grounded file-read check. |
| `launch.sh` | Split as proposed. The PID-safe waiter, output capture, exit reporting, and byte counts belong in an asset. The worker/path table is per-run input. |
| Worker prompts | Disposable by default. A prompt is an execution projection of a plan, current state, and guardrails—not automatically project memory. Preserve its stable decisions, ownership partition, acceptance baseline, and deviations in the project’s declared session/reasoning record. Retain the raw prompt only when it is itself the authoritative spec or required audit evidence. |
| Logs, byte counts, `-o` messages | Disposable. Promote only findings, deviations, or decisions that changed the work. Stateful research output remains the existing exception; the skill already routes that into reasoning records. |

The worked instance supports this distinction. Its durable `00_STATE.md` already records worker ownership, interpreter/`.pth` findings, the dot-directory glob trap, test counts, the GitHub-only tag, the 18.8 GB exclusion, and unresolved owner decisions ([00_STATE.md](/data1/users/antonz/projects/Meta-Aging/14782-DM/docs/_internal/plans/2026-08-14_consensus-migration/00_STATE.md:3)). Archiving three absolute-path/model-specific prompts would largely duplicate that record.

Therefore, do not change the three `/tmp` examples to an in-tree prompt route. Keep the per-run prompt rendering in temporary storage; ensure its durable source and residue exist in the project memory structure. This avoids adding another Markdown accumulation channel to trees already shown to mix decisions with environments, caches, checkpoints, and logs ([FINDINGS_field.md](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/docs/_internal/research/2026-08-20-internal-skeleton/FINDINGS_field.md:49)).

3. Ownership of normalization

Scio owns it: the `delegate-cli` skill defines the workflow and assets, while the project’s scio-owned memory grammar defines where durable residue goes. The catalog mount links the complete `skills/` tree, so `assets/` becomes reachable without another delivery mechanism ([link.sh](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/lib/scio/link.sh:140)).

The verification stance belongs in `delegate-cli` as part of the delegation workflow. The already-settled delegation seam says scio owns the workflow/fallback and dev-env owns nothing here ([07_delegation-seam.md](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/docs/_internal/plans/2026-08-20-memory-and-seams/07_delegation-seam.md:14)). Scbio-docker owns the container sandbox cause and supported configuration; the skill should cite that source while retaining the invocation and review procedure.

No delegation-specific lint mechanism is warranted. The planned general internal-memory lint can reject obvious logs, caches, checkpoints, and other non-memory payloads. Dev-env should neither duplicate the stance nor name a project path.

4. What is enforceable

Mechanically enforceable within the supplied scripts:

- waiting on captured PIDs;
- unique output paths;
- exit-status and byte-count reporting;
- failure when expected outputs are absent;
- recording pre/post HEADs, reflogs, and worktree state;
- catalog asset existence and mount reachability.

The following are not enforceable by `scio lint` or a shell launcher without introducing semantic metadata and judgment:

- whether a prompt contains every important decision;
- whether its path partition is conceptually disjoint;
- whether an acceptance command is correct;
- whether a changed assertion was legitimately strengthened or improperly weakened;
- whether review was genuinely independent;
- whether the durable residue captured everything worth retaining;
- whether a Markdown prompt is memory or landfill.

Under the current unsandboxed codex workaround, “do not run git” is also not enforced. Log, HEAD, reflog, and status checks provide post-run evidence; prompt prohibitions do not recreate a sandbox. Lint cannot inspect vanished `/tmp` state. No hook should be added.

5. Resolving the duplication

`delegate-cli` is the single home for the CLI flags, container workflow, prompt traps, and verification stance. Meta-Aging’s duplicated operational section exists because its catalog was unreachable, not because it owns newer knowledge ([AGENTS.md](/data1/users/antonz/projects/Meta-Aging/AGENTS.md:203)).

Once the fleet sweep has mounted the catalog and reachability is verified, replace that section with at most one router sentence—“Use the `delegate-cli` skill for codex/agy delegation”—or remove it entirely if skill triggering is reliable. Delete the duplicated facts and stance. Until reachability lands, the current copy remains a necessary temporary bridge.

6. First bounded move

Make one edit to `skills/delegate-cli/SKILL.md` adding an explicit retention contract:

> Runtime prompt renderings, logs, final messages, and byte counts stay in temporary storage. Before launch, put stable decisions, worker ownership, and acceptance baselines in the repository’s declared durable-memory record; after review, add material deviations and outcomes. Preserve raw prompts/output only when they are authoritative specifications, research evidence, or required audit records.

That single-file change corrects the ambiguity behind the loss without moving artifacts, adding a directory, creating metadata, or depending on the asset work or fleet sweep. No files were modified during this consultation.