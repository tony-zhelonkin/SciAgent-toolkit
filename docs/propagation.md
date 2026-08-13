# How a change reaches a project

Toolkit changes reach a project through versioning and binding.

```text
toolkit commit
    |
    | hop 1: re-pin the project's toolkit checkout
    v
<project>/01_modules/SciAgent-toolkit/
    |
    | hop 2: run link or craft for materialized project state
    v
<project>/.claude/  <project>/.agents/  <project>/AGENTS.md
```

## Hop 1: version

Consumer projects hold this toolkit as a Git submodule pinned to a commit.
Re-pinning changes the bytes available in that project's own toolkit checkout.
The project runs the CLI at its pin.

When `$SCIAGENT_TOOLKIT` points directly at the checkout being edited, each
save has already completed hop 1 for that working tree.

## Hop 2: binding

Two commands materialize project state from the pinned checkout:

- `sciagent link` creates the six category links and refreshes guardrail hook
  bodies, registrations, and the `SCIAGENT:GITIGNORE` block.
- `sciagent craft` renders `SCIAGENT:CRAFT` into `AGENTS.md`.

Both commands accept `--project-dir D` and are idempotent.

## Artifact classes

### Linked catalog content

The project has one directory link for each harness and category:

| Project link | Toolkit directory |
|---|---|
| `.claude/skills`, `.agents/skills` | `skills/` |
| `.claude/agents`, `.agents/agents` | `agents/` |
| `.claude/commands`, `.agents/commands` | `commands/` |

Edits, additions, removals, and renames inside these directories propagate on
hop 1. Directory membership stays live through the category link, so catalog
changes need no relinking after the initial binding.

### Materialized project content

The project owns copies or rendered text:

| Project path | Producer |
|---|---|
| `.claude/hooks/*.sh` | `sciagent link` |
| `.claude/settings.json` hook registrations | `sciagent link` |
| `AGENTS.md` `SCIAGENT:CRAFT` block | `sciagent craft` |
| `.gitignore` `SCIAGENT:GITIGNORE` block | `sciagent link` |

These require hop 1 followed by the producing command. A hook-template change
needs `link`; a `craft.yaml` change needs `craft`.

### Toolkit code

`bin/sciagent` and `lib/sciagent/*.sh` execute from the pinned checkout. Their
behavior changes on hop 1 and takes effect at the next invocation.

## Performing propagation

For a vendored project:

```bash
git submodule update --init 01_modules/SciAgent-toolkit
./01_modules/SciAgent-toolkit/bin/sciagent link
./01_modules/SciAgent-toolkit/bin/sciagent craft
```

The project owner chooses the new submodule commit before these commands. The
toolkit has no command that re-pins a consumer automatically.

For a release install, `install.sh` places an immutable toolkit version on the
machine. Running that version's `sciagent link` performs the project binding.
A project that vendors its own toolkit must use its in-repo binary; `link`
refuses a global copy in that situation.

## Hook ownership

Hook bodies are copied byte-for-byte from templates. Each managed body carries
a content-hash record under `.sciagent/hook_state/`.

| Evidence | `link` action |
|---|---|
| destination absent | write current template and record its hash |
| content matches the current template | adopt silently |
| content matches a shipped template version or its ownership record | refresh and report |
| content matches neither | preserve, warn once, and write a ceded marker |

`templates/PROVENANCE.sha1` records shipped hook hashes, allowing an older
unmodified body to be recognized in a shallow clone or release archive. A
ceded marker suppresses repeated warnings for the same user-owned content and
allows a later, different edit to be reported once.

## Legacy mount convergence

Before creating category links, `link` scans the old flat mount locations. It
removes symlinks that resolve inside the active toolkit and dangling legacy
links whose path identifies a `SciAgent-toolkit` source. This covers historical
absolute `/workspaces/...` mounts and the retired `.claude/output-styles/`
tree. Regular files and outside-pointing links remain.

After toolkit-owned child links are swept, an empty category directory is
replaced by the category link. A directory with any remaining entry is
preserved and refused with a listing and relocation instructions.

## Manual removal

Project teardown is explicit:

```bash
rm .claude/{skills,agents,commands} .agents/{skills,agents,commands}
# Remove complete SCIAGENT managed blocks from AGENTS.md.
```

Hook bodies and registrations are materialized guardrails. Retiring those is a
separate manual edit so user-customized hooks and unrelated settings remain
under project control.
