Recommendation: keep `docs/_internal/` as an independent Git repository during work, preferably backed by a private remote. At publication, leave the memory repository separate and commit a derived `docs/internal-memory.md` into the parent. If publication is approved, that file names an anonymously accessible memory repository at an immutable ref and full commit ID. If publication is declined, it records only that the omission is intentional.

This preserves the owner’s late choice, keeps daily memory commits out of the parent, and avoids the lasting dirty-worktree and clone behavior of a submodule.

## 1. Ranked topologies

| Rank | Topology | Daily cost | Publication cost | Main failure |
|---|---|---|---|---|
| 1 | Nested repository at `docs/_internal/`; public remote and parent pointer added only after approval | Commit memory separately. Parent remains clean. A private remote is advisable for backup. | Review all reachable history, publish an immutable memory ref, then commit the derived pointer file. | A repository without commits or an off-machine remote still loses data. A branch-only pointer can drift or disappear. |
| 2 | Submodule registered and pinned only at publication | Before first publication it is simply an ignored nested repo. After registration, advancing its HEAD makes the parent appear dirty even without pin commits. | Create and populate the public remote first; add `.gitmodules` and the gitlink at the reviewed commit. | The pin is deliberately stale between publications; ordinary clones leave it uninitialized; recursive clones fail when access is unavailable. |
| 3 | Separate repository outside `docs/_internal/`, named by a pointer | Requires a checkout, worktree, symlink, or agent routing convention to connect the required path to the external repository. | Publish and write the same immutable pointer as rank 1. | The pointer can name a repository that was never checked out at the promised path. If its working tree is placed at `docs/_internal/`, this collapses into rank 1. |
| 4 | Snapshot committed into the public parent at release | No daily parent churn. | Copy the reviewed tip into the parent and make one potentially large release commit. | It loses the history needed to archive overwritten `session.md` states and makes the disclosed snapshot permanent in parent history. |
| 5 | `git subtree` | Can defer parent changes until publication if the source stays separate. | Import or squash the memory into the parent. Repeated releases require subtree updates. | Full import floods parent history; squash reduces to the snapshot option. Separation is lost after publication. |
| 6 | Orphan branch in the parent repository | Requires branch switching or a linked worktree to edit memory while preserving the main checkout. Default-branch history remains clean. | Publish an allowlisted set of refs and deliberately include or exclude the memory branch. | It mixes public and private objects in one repository, making `push --all`, mirrors, bundles, migrations, and visibility changes disclosure hazards. |

### Does publication-only pinning dissolve the submodule objection?

Yes. It dissolves the original objection about frequent parent-log pin bumps: the parent changes once per publication, which is the correct coupling point.

Other costs replace it:

- Before the first publication, there can be no valid submodule entry until its repository exists and the pinned commit is reachable.
- After registration, daily memory commits move the submodule HEAD and normally leave the parent working tree reporting a modified gitlink. `ignore=all` can suppress that signal, but then stale state becomes invisible.
- The public pin represents the last published memory state, intentionally lagging current private memory.
- Plain `git clone` does not initialize submodules. `git clone --recurse-submodules` or a later update does.
- A reader without access cannot initialize the submodule. Therefore a private submodule must never be the public repository’s “no” answer.
- The pinned commit must remain reachable from a durable advertised ref; relying on a raw object ID alone is fragile.

For these reasons, publication-only submodules are defensible but remain second choice.

### Orphan branch

The premise needs one correction: making a repository public exposes every branch already present on the server; it does not automatically push every local branch. The hazard remains substantial.

It is salvageable only if publication uses a sterile staging repository created from an explicit allowlist of public refs, and the public remote is never used directly from the repository containing the memory branch. That ensures memory-only objects are absent from the pushed object graph. A branch-name deny rule is insufficient because a tag or renamed branch could expose the same objects.

Without a hook or server-side allowlist, this safety property depends on always using the publication procedure. The separate object database of rank 1 is a stronger boundary with less machinery.

## 2. What Scio should change

Scio should provide:

- The settled documentation describing the daily nested-repository arrangement and the two publication outcomes.
- One corrected `internal-memory` history predicate.
- Tests for real repositories, Git worktree `.git` files, parent-tracked trees, initialized-but-unborn repositories, and unprotected populated trees.

Scio should not add a hook, scaffold, pointer generator, remote manager, or publication verb. `link` remains a binding operation, `craft` renders standing text, and `lint` remains read-only. Repository initialization, remote creation, redaction, and the parent publication commit belong to the owner-operated publication procedure.

The pointer is created by that procedure from Git-derived values. It is release metadata generated at the decision point, rather than daily metadata someone must manually synchronize.

## 3. Publication sequence

### Answer: yes, publish memory

1. Freeze edits to both repositories.
2. Run `scio lint --check internal-memory --strict`.
3. Confirm the memory repository has no untracked or uncommitted content.
4. Select the exact memory ref to publish. Mechanically inspect every object, commit message, tag message, filename, and author identity reachable from that ref.
5. Perform the human redaction review. If full history cannot be cleared confidently, publish a sanitized single-root snapshot repository and explicitly record that its private history was withheld.
6. Push the reviewed memory commit to a public repository under an immutable tag or release ref.
7. From an unauthenticated environment, clone that URL and check out the recorded commit.
8. Create `docs/internal-memory.md` in the parent from the verified URL, ref, and full object ID.
9. Confirm `docs/_internal/` remains absent from the parent index and from every parent ref being published.
10. Publish only the intended parent refs, then perform a fresh recovery test starting from the public parent clone.

### Answer: no, memory remains private

1. Keep the memory repository and any private remote separate.
2. Confirm that all public parent refs contain no `docs/_internal/` files, gitlink, `.gitmodules` entry, memory branch, memory tag, or private remote URL.
3. If an orphan branch ever shared the parent object database, publish from a clean staging repository containing only allowlisted public refs.
4. Create the status-only `docs/internal-memory.md` described below.
5. Publish the parent and verify from an unauthenticated fresh clone that no memory content or retrieval endpoint appears.

A previous “yes” cannot later be undone reliably. Future updates may stop, but disclosed commits can remain in clones, archives, and caches.

## 4. Redaction map

| Mechanically checkable | Requires human review |
|---|---|
| Allowed file types, blob-size limits, binaries, symlinks, submodules, and unexpected object types | Whether scientific reasoning, hypotheses, results, or rejected alternatives are unpublished or confidential |
| Known credential prefixes, private-key headers, URLs containing credentials, and other fixed high-confidence secret patterns | Credentials or access instructions that do not match known syntax |
| Absolute Unix, lab-storage, container, home-directory, Windows-drive, and UNC path syntax | Whether a syntactically ordinary path reveals a lab, project, dataset, or collaborator |
| Email-address and author-identity inventory | Whether names, emails, and collaborator identities may be disclosed |
| Every reachable historical file, deleted file, commit message, tag message, filename, and ref—not merely the current checkout | Patient-adjacent identifiers, cohort descriptions, re-identification combinations, and context-dependent IDs |
| Parent refs contain no ignored memory, gitlink, private URL, or unintended branch | Whether the overall narrative is safe and appropriate to publish |
| Anonymous clone and exact-revision recovery | Whether publication consent and licensing are adequate |

Mechanical scanning should report paths and revisions without printing matched secret values. Any discovered credential must be removed from all published history and rotated.

No generic Bash predicate can determine that a person’s name, numeric identifier, cohort description, or scientific statement is safe. Those remain human review.

## 5. Pointer-file contract

For a public memory repository, `docs/internal-memory.md` must contain:

- Status: `published`.
- Canonical anonymously cloneable URL.
- Full Git object ID.
- Immutable tag or release ref that keeps that object reachable.
- Publication date and corresponding parent release or commit.
- Intended checkout path: `docs/_internal/`.
- Exact clone and detached-checkout commands.
- Whether the publication includes complete history or a sanitized snapshot.
- Applicable license or access statement.
- Optional archival mirror or DOI.

A moving branch name is useful context but cannot replace the immutable ref and object ID.

For private memory:

- Status: `not published`.
- Decision date and corresponding parent release.
- A clear statement that the omission is intentional and no public retrieval endpoint is declared.
- Optionally, a public contact route such as the repository issue tracker.

It must contain no private URL, local absolute path, repository name, commit ID, access token, or dead placeholder link. This is a status record, not an unresolvable pointer.

## 6. Lint predicate

The current in-progress predicate in [lint.sh](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/lib/scio/lint.sh:735) is insufficient:

- `-d "$root/.git"` misses valid Git worktrees and submodules where `.git` is a file.
- A directory merely named `.git` is accepted without verifying Git.
- The current test accepts `git init` without checking that the path is an actual repository root.
- It only considers `session*.md`, although the settled ruling says a populated memory tree.
- It does not recognize memory genuinely tracked by the parent.

An implementable Bash predicate is:

1. An absent `docs/_internal/` is silent.
2. “Populated” means at least one nonempty allowed memory Markdown file exists below it, pruning repository machinery. Repository-control files alone do not count.
3. Resolve physical paths with `pwd -P`.
4. Ask Git, rather than inspecting `.git`:

   - Run `git -C docs/_internal rev-parse --show-toplevel`.
   - The returned physical path must equal the physical memory root.
   - This recognizes nested repositories, linked worktrees, submodules, and separate Git directories.
   - An unborn but valid dedicated repository may count as the brand-new grace state.

5. Otherwise, test parent coverage:

   - Resolve the parent Git root.
   - Require every substantive memory Markdown file to have a committed parent path. A handful of force-tracked files among thousands does not qualify.
   - An ignored, untracked tree therefore fails.

6. If neither condition holds, emit:

   `docs/_internal/ is populated and has no detectable Git history mechanism; initialize or attach history before updating mutable memory in place`

7. Do not test for a remote here. A remote concerns reachability and publication; local history is sufficient for this specific archive predicate.

There is no exact way to distinguish a newly created unversioned directory from an old unversioned directory using only its current tree and no stored state. Filesystem timestamps provide only a heuristic and are reset or preserved by copy, restore, and `touch`. The clean deterministic answer is for the documented first-content procedure to initialize the nested repository at the same time it creates the first authentic file. That makes a valid unborn repository the objectively detectable “brand new” case.

## 7. Where I disagree

- “Publishing the repository publishes every branch” is too broad. Publication exposes server-side refs; the practical warning about orphan branches remains correct.
- Waiting until publication to create any remote leaves daily memory vulnerable to disk loss. A private remote should exist earlier when feasible; visibility is the late decision.
- Presence of a `.git` directory is not evidence of a history mechanism. Git must verify the worktree root.
- A redaction scan cannot certify publication safety.
- Exact reproducibility requires updating either a pointer commit or submodule pin at each publication. Zero parent changes across publications and an exact parent-to-memory binding cannot both hold.

Verified by reading: ADR-D10, current and in-progress lint implementation/tests, the field digest, the accepted minimalism consultation, architecture, handoff, the managed ignore block, and installed Git clone/worktree help. The ranking, procedures, and publication-risk analysis are my inference. The supplied fleet measurements were accepted from the cited field report rather than re-derived.

No files were modified.