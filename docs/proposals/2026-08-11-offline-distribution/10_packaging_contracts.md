# Packaging contracts — `build-release.sh` and `install.sh`

Two small Bash scripts replace Phase 6's `plugin.json` + release-workflow + `npx` proof. Neither
touches the network. Neither mutates a project.

The division of labour is absolute and is the whole point of the design:

| Script | Owns | Never does |
|---|---|---|
| `scripts/build-release.sh` | producing bytes from a Git ref | fetch, publish, install, touch a project |
| `install.sh` | placing bytes on this machine | resolve URLs, fetch, choose a harness, touch a project |
| `sciagent` itself | mutating projects (`activate`, later `bind`) | install itself |

---

## 1. `scripts/build-release.sh`

### Contract

- **Requires an explicit clean commit or ref.** Refuse on a dirty tree; refuse on an implicit
  "whatever is checked out." The ref is an argument, not an inference.
- **Uses `git archive`**, so ignored content — `.venv`, caches, development residue — cannot enter
  the artifact. This is mandatory, not stylistic: the working checkout is 180 MB against 5.6 MB of
  tracked content, and 171 MB of the difference is one skill's virtualenv.
- **Produces exactly three things:**
  ```
  sciagent-<tag>-<short-sha>.tar.gz
  sciagent-<tag>-<short-sha>.tar.gz.sha256
  release metadata recording the COMPLETE git SHA (not abbreviated)
  ```
- **Runs `./bin/sciagent validate` and `bash tests/run-all.sh` before producing the artifact.** A
  release that fails its own suite must not exist as a file.
- **Deterministic:** a test builds the same ref twice and asserts identical checksums.
- **No network operations whatsoever.**

### Determinism, verified

`git archive HEAD` is already byte-stable across runs. The compression step is not, unless told:

```bash
git archive --format=tar "$ref" | gzip -n > "$out"    # -n omits mtime → stable
git archive --format=tar.gz "$ref"                    # NOT stable: embeds an mtime
```

Both forms were measured on this checkout; only the first repeats its SHA256. `--format=tar.gz` is
therefore unusable for any reproducibility claim.

### Why the full SHA in metadata

The short SHA in the filename is for humans. The full SHA in the metadata is what a reader of the
paper resolves back to a commit, and what `install.sh` uses as the content-addressed directory name.
Abbreviated SHAs collide as history grows; a supplement that ships one is a supplement that rots.

---

## 2. `install.sh`

### Invocation — local inputs only

```bash
./install.sh --archive  sciagent-0.1.0-<sha>.tar.gz \
             --checksum sciagent-0.1.0-<sha>.tar.gz.sha256 \
             --prefix   "$HOME/.local"
```

A local checkout is equally acceptable as input. Nothing else is.

### Contract

- **No URL parsing.** Not "URLs are discouraged" — there is no code path that accepts one.
- **No `curl`, no `npm`, no telemetry, no registry, no automatic update check.**
- **Verifies the checksum before extracting anything.**
- **Installs atomically into a content-addressed directory** keyed by the full Git SHA:
  ```
  ~/.local/share/sciagent/versions/<full-git-sha>/
  ```
  Two versions coexist by construction; an interrupted install leaves no half-tree.
- **Links only the executable** into `~/.local/bin`. Nothing else escapes the version directory.
- **Writes an installation receipt** sufficient for an exact uninstall — same ownership-record
  discipline the toolkit now uses everywhere else (symlinks carry their target, managed blocks carry
  hash-framed markers, settings artifacts carry state tags). An artifact without a record is an
  artifact that cannot be reversed, which is the class of bug Phase 5c and the `deactivate` work
  existed to eliminate.
- **Supports `--dry-run`.**
- **Never mutates a project.** Installation and project binding are different verbs on different
  layers.

### What is deliberately absent

No harness detection. No `settings.json`. No `AGENTS.md`. No skill mounting. `install.sh` puts a
program on the machine; `sciagent` decides what a project gets. Merging those is exactly the mistake
§2 of `00_INDEX.md` names.

---

## 3. Project binding stays in `sciagent`

Eventually the project-mutating verb becomes `bind` (or `link`), with the harness set explicit:

```bash
sciagent bind --project . --harness portable
sciagent bind --project . --harness portable,claude
```

`portable` = the common layer only: `AGENTS.md` + `.agents/skills` + the ownership receipt. Adding
`claude` adds `CLAUDE.md` + `.claude/*` + hooks/settings. This is ADR-D5's adapter split surfaced as
a flag, and it is what makes rule 1 ("packaging must not select or configure harnesses") checkable:
the harness set is always something the user typed.

**Fleet precedence is unchanged and must stay that way:** a project-local
`01_modules/SciAgent-toolkit` overrides any global installation. The locality guard already enforces
this distinction — `bin/sciagent` sources its modules from `$SCIAGENT_TOOLKIT`, so pointing that
variable elsewhere runs *that* checkout's code, and `_guard_toolkit_locality` plus the inline guard
in `deactivate.sh` refuse the mismatch. Exact-commit reproducibility depends on that precedence, not
on convention.

---

## 4. Distribution methods under this model

| Method | Network / telemetry | Reproducibility | Best use |
|---|---|---|---|
| Existing submodule + local hub | none | exact SHA | **the analysis fleet** |
| `git bundle` | none during transfer or install | full history and tags | air-gapped / lab transfer |
| Release tarball + SHA256 | none during install | immutable snapshot | paper supplement, collaborators |
| Direct git clone from chosen host | only the chosen Git host | tag / SHA pin | technical external users |
| `npx skills add` with opt-out | npm + Git host; telemetry disabled by env | client lock semantics | optional compatibility test only |
| Claude marketplace | Claude-specific | vendor lifecycle | no reason to adopt now |

For publication, the tarball and its checksum can be deposited with the paper or in an archival
repository. That gives readers an immutable artifact **without making any marketplace the package
authority** — which is the property the whole design is protecting.
