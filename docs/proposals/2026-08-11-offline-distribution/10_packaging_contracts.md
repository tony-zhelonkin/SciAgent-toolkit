# Packaging contracts — `build-release.sh` and `install.sh`

Two small Bash scripts replace Phase 6's `plugin.json` + release-workflow + `npx` proof. Neither
touches the network. Neither mutates a project.

The division of labour is absolute and is the whole point of the design:

| Script | Owns | Never does |
|---|---|---|
| `scripts/build-release.sh` | producing bytes from a Git ref | fetch, publish, install, touch a project |
| `install.sh` | placing bytes on this machine | resolve URLs, fetch, choose a harness, touch a project |
| `sciagent` itself | project binding, managed context, and linting | install itself |

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
  scio-<version>-<short-sha>.tar.gz
  scio-<version>-<short-sha>.tar.gz.sha256
  scio-<version>-<short-sha>.metadata.json recording the complete Git SHA
  ```
- **Runs `./bin/sciagent lint --check toolkit` and `bash tests/run-all.sh` before producing the artifact.** A
  release that fails its own suite must not exist as a file.
- **Deterministic:** a test builds the same ref twice and asserts identical checksums.
- **No network operations whatsoever.**

### Determinism, verified — **and this section's first draft was wrong**

`git archive HEAD` is byte-stable across runs. So, it turns out, is every compression form below.

```bash
git archive --format=tar "$ref" | gzip -n > "$out"    # stable  — what we ship
git archive --format=tar "$ref" | gzip    > "$out"    # ALSO stable
git archive --format=tar.gz "$ref"        > "$out"    # ALSO stable
```

**Superseded claim, kept visible:** this section previously stated that only the `gzip -n` form
repeats its SHA256, and that `--format=tar.gz` is "NOT stable: embeds an mtime" and "therefore
unusable for any reproducibility claim." Re-measured 2026-08-11 during implementation and confirmed
independently: all three forms produce identical SHA256 across runs, one second apart.

The mechanism explains why, and it is not host luck. **gzip embeds an MTIME only when it compresses
a NAMED FILE.** Reading from a pipe there is no name and no mtime, so it writes zero into that
header field — which is exactly what `-n` forces. Git's internal `--format=tar.gz` shells out to
`gzip -cn`, so it never embedded one either. The original measurement most likely compared
`gzip < file` against a named-file invocation.

`gzip -n` is still what the script uses. It is now belt-and-braces against a gzip implementation
that behaves differently rather than a fix for a defect that was ever observed here, and the script
says so. Anything built on the retracted claim — for instance rejecting `--format=tar.gz` on
reproducibility grounds — should be revisited.

### Why the full SHA in metadata

The short SHA in the filename is for humans. The full SHA in the metadata is what a reader of the
paper resolves back to a commit, and what `install.sh` uses as the content-addressed directory name.
Abbreviated SHAs collide as history grows; a supplement that ships one is a supplement that rots.

---

## 2. `install.sh`

### Invocation — local inputs only

```bash
./install.sh --archive  scio-0.1.0-<sha>.tar.gz \
             --checksum scio-0.1.0-<sha>.tar.gz.sha256 \
             --prefix   "$HOME/.local"
```

A local checkout is equally acceptable as input. Nothing else is.

### Contract

- **No URL parsing.** Not "URLs are discouraged" — there is no code path that accepts one.
- **No `curl`, no `npm`, no telemetry, no registry, no automatic update check.**
- **Verifies the checksum before extracting anything.**
- **Installs atomically into a content-addressed directory** keyed by the full Git SHA:
  ```
  ~/.local/share/scio/versions/<full-git-sha>/
  ```
  Two versions coexist by construction; an interrupted install leaves no half-tree.
- **Links only the executable** into `~/.local/bin`. Nothing else escapes the version directory.
- **Writes an installation receipt** sufficient for an exact uninstall. The
  installer removes only paths the receipt proves it wrote; project catalog
  links carry ownership in their targets, and managed bodies carry recognized
  hashes or markers.
- **Supports `--dry-run`.**
- **Never mutates a project.** Installation and project binding are different verbs on different
  layers.

### What is deliberately absent

No harness detection. No `settings.json`. No `AGENTS.md`. No skill mounting. `install.sh` puts a
program on the machine; `sciagent` decides what a project gets. Merging those is exactly the mistake
§2 of `00_INDEX.md` names.

---

## 3. Project binding stays in `sciagent link`

Project binding is one convergent operation:

```bash
sciagent link --project-dir .
```

`link` creates six whole-tree category links under `.agents/` and `.claude/`,
materializes the two Claude guardrail hooks, merges their settings
registrations, and refreshes the managed gitignore block. `craft` separately
renders the `SCIAGENT:CRAFT` block into `AGENTS.md`. The local installer still
performs neither operation, preserving the packaging/project boundary from
ADR-D3.

**Fleet precedence is unchanged and must stay that way:** a project-local
`01_modules/SciAgent-toolkit` overrides any global installation. The locality guard already enforces
this distinction: `link` detects the in-project toolkit and refuses an active
external copy. Exact-commit reproducibility depends on that precedence.

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
