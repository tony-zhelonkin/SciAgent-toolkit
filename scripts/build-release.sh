#!/usr/bin/env bash
# scripts/build-release.sh — produce a deterministic, offline release artifact
# from an explicit Git ref.
#
# Implements `docs/proposals/2026-08-11-offline-distribution/10_packaging_contracts.md` §1.
# Owns exactly one thing: turning a Git ref into bytes. It never fetches, never
# publishes, never installs, and never touches a project.
#
# ---------------------------------------------------------------------------
# Usage
# ---------------------------------------------------------------------------
#   scripts/build-release.sh <ref> [--out DIR] [--version V] [--quiet]
#
# The ref is a REQUIRED positional argument. "Whatever is checked out" is
# refused on purpose: a release is a statement about a commit, and inferring
# which commit makes that statement unfalsifiable.
#
# ---------------------------------------------------------------------------
# What it produces (exactly three files, per the contract)
# ---------------------------------------------------------------------------
#   <out>/scio-<version>-<short-sha>.tar.gz
#   <out>/scio-<version>-<short-sha>.tar.gz.sha256      (sha256sum -c format)
#   <out>/scio-<version>-<short-sha>.metadata.json      (carries the FULL sha)
#
# A fourth copy of the metadata — minus the self-referential digest fields —
# rides INSIDE the tarball as `.scio-release.json`, because `install.sh` takes
# only a local archive and needs the full commit SHA to name its
# content-addressed install directory. Deriving that SHA from the filename's
# abbreviation would reintroduce exactly the collision the contract's "full SHA
# in metadata" rule exists to avoid.
#
# ---------------------------------------------------------------------------
# Naming — the one place this file departs from a literal reading of the docs
# ---------------------------------------------------------------------------
# `10_packaging_contracts.md` §1 writes the filename as
# `sciagent-<tag>-<short-sha>.tar.gz`; ADR-D7 writes `scio-0.1.0.tar.gz`. The two
# disagree in two independent ways, and both are resolved here deliberately:
#
#   stem      — `scio` wins. ADR-D7 is the later and explicitly DECIDED record,
#               and it decided the stem for a hard reason (`sciagent` is taken).
#   short sha — KEPT, against ADR-D7's example. §1's rationale is the sound one:
#               the abbreviation is for humans reading a filename, the full SHA
#               in the metadata is the identity. Dropping it makes two builds of
#               different commits at the same version indistinguishable on disk,
#               which is precisely the failure mode a paper supplement cannot
#               afford. ADR-D7's `scio-0.1.0.tar.gz` is an illustration of the
#               NAME question it was deciding, not a decision about the SHA.
#
# SCOPE LIMIT, load-bearing: ADR-D7 renames the release ARTIFACT only. The CLI
# stays `sciagent` — the installed executable, `$SCIAGENT_TOOLKIT`, the `si`
# alias and project-locality contract are all untouched by this script. Only the
# tarball and its metadata carry the `scio` name.
#
# ---------------------------------------------------------------------------
# Determinism
# ---------------------------------------------------------------------------
#   git archive --format=tar "$sha" | gzip -n > out
#
# `-n` is mandatory by contract and kept even though it is belt-and-braces on
# this host. Measured here (git 2.34.1, GNU gzip 1.10), and this CORRECTS
# `00_INDEX.md` §5's recorded finding:
#
#   * `git archive --format=tar HEAD | gzip -n`  → stable  (as recorded)
#   * `git archive --format=tar HEAD | gzip`     → ALSO stable, byte-identical
#   * `git archive --format=tar.gz HEAD`         → ALSO stable, byte-identical
#
# All three produced the same SHA256 across runs separated in time. The reason
# the recorded claim ("--format=tar.gz embeds an mtime") does not reproduce is
# mechanical: gzip embeds an MTIME only when it compresses a NAMED FILE; reading
# a pipe it writes MTIME=0. And git's built-in `tar.gz` filter is literally
# `gzip -cn`, so it never embedded one either.
#
# `-n` and the explicit pipe are still the right implementation, for a reason
# that survives the correction: `tar.tar.gz.command` is user-configurable git
# config, so `--format=tar.gz` inherits determinism from the AMBIENT
# ENVIRONMENT, while the explicit pipe does not. A reproducibility claim must
# not depend on the builder's `~/.gitconfig`.
#
# ---------------------------------------------------------------------------
# The pre-flight gate, and why it runs against the EXPORT
# ---------------------------------------------------------------------------
# `./bin/sciagent validate` and `bash tests/run-all.sh` run against a clean
# export of <ref>, not against the working tree. The working tree may sit at a
# different commit than <ref>, and it carries untracked residue the artifact
# will not; testing it would certify bytes nobody ships. The export is exactly
# the shipped bytes.
#
# RECURSION: this script runs the suite, so a test that ran THIS script against
# THIS repo would re-enter the suite. There is deliberately no `--skip-checks`
# escape hatch — an escape hatch on a release gate is a footgun that eventually
# gets used for a real release. Instead `tests/test_build_release_*.sh` build
# from throwaway fixture repos (`tests/_release_lib.sh`) whose `bin/sciagent`
# and `tests/run-all.sh` are instant stubs. The gate code path is therefore
# exercised in both directions (pass and fail) in milliseconds, and recursion is
# impossible by construction rather than by flag.
#
# No network operations. Nothing here resolves a name, opens a socket, or
# consults a registry.
set -euo pipefail

# --- constants -------------------------------------------------------------
ARTIFACT_NAME="scio"          # ADR-D7 (release artifact only; CLI stays `sciagent`)
DEFAULT_VERSION="0.1.0"       # ADR-D7 starting version
SHORT_LEN=12                  # abbreviation length in the human-facing filename
EMBEDDED_METADATA=".scio-release.json"
METADATA_SCHEMA=1

_prog="${0##*/}"

die() { echo "$_prog: $*" >&2; exit 1; }
say() { [[ -n "${QUIET:-}" ]] || echo "$*"; }

usage() {
    cat <<EOF
$_prog — build a deterministic offline release artifact from an explicit Git ref

Usage:
  $_prog <ref> [--out DIR] [--version V] [--quiet]

Arguments:
  <ref>              REQUIRED. Any commit-ish (tag, branch, SHA). There is no
                     default: an implicit "whatever is checked out" is refused.

Options:
  --out DIR          Output directory (default: <repo>/dist)
  --version V        Release version. Default: the tag pointing exactly at <ref>,
                     else $DEFAULT_VERSION (ADR-D7) with a warning.
  --quiet            Only print the produced paths.
  -h, --help         This message.

Produces (in DIR):
  ${ARTIFACT_NAME}-<version>-<short-sha>.tar.gz
  ${ARTIFACT_NAME}-<version>-<short-sha>.tar.gz.sha256
  ${ARTIFACT_NAME}-<version>-<short-sha>.metadata.json   (full 40-char SHA)

Refuses: a dirty working tree, a missing/unresolvable ref, a ref whose exported
tree fails \`sciagent validate\` or \`tests/run-all.sh\`. No network access.
EOF
}

# sha256_of <file> — print the bare hex digest.
sha256_of() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | awk '{print $1}'
    elif command -v shasum >/dev/null 2>&1; then
        shasum -a 256 "$1" | awk '{print $1}'
    else
        die "neither sha256sum nor shasum found; cannot checksum the artifact"
    fi
}

# _json_escape <string> — escape for embedding between JSON double quotes.
# Kept local because this script must run from a bare export with no toolkit sourcing.
_json_escape() {
    local s="$1"
    s="${s//\\/\\\\}"
    s="${s//\"/\\\"}"
    s="${s//$'\n'/\\n}"
    s="${s//$'\r'/\\r}"
    s="${s//$'\t'/\\t}"
    printf '%s' "$s"
}

# --- argument parsing ------------------------------------------------------
ref=""
out_dir=""
version=""
QUIET=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help) usage; exit 0 ;;
        --out)
            [[ $# -ge 2 ]] || die "--out requires a directory"
            out_dir="$2"; shift 2 ;;
        --version)
            [[ $# -ge 2 ]] || die "--version requires a value"
            version="$2"; shift 2 ;;
        --quiet) QUIET=1; shift ;;
        --) shift; break ;;
        -*) die "unknown option: $1 (see --help)" ;;
        *)
            [[ -z "$ref" ]] || die "unexpected extra argument: $1 (exactly one ref)"
            ref="$1"; shift ;;
    esac
done

if [[ -z "$ref" ]]; then
    echo "$_prog: a ref is required — refusing to guess what to release." >&2
    echo "  A release names a commit. Inferring it from the checkout makes the" >&2
    echo "  artifact's identity depend on the builder's working directory." >&2
    echo "  Try:  $_prog HEAD          (explicit, still a real ref)" >&2
    echo "        $_prog v0.1.0" >&2
    exit 2
fi

# --- repo resolution -------------------------------------------------------
# The repo is this script's own repo, not the caller's CWD: build-release.sh
# releases the toolkit it ships with, never whatever directory it is invoked
# from. (Tests copy the script into a fixture repo for exactly this reason.)
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"

command -v git >/dev/null 2>&1 || die "git not found"
git -C "$repo_root" rev-parse --git-dir >/dev/null 2>&1 \
    || die "not a git repository: $repo_root"

# --- refuse a dirty tree ---------------------------------------------------
# `git archive` reads the ref, not the worktree, so a dirty tree cannot corrupt
# the bytes. It is refused anyway, and the reason is the gate rather than the
# archive: uncommitted work means the thing the builder is looking at is not the
# thing being certified, and any "I tested it" claim silently drifts.
dirty="$(git -C "$repo_root" status --porcelain)"
if [[ -n "$dirty" ]]; then
    echo "$_prog: refusing to build from a dirty working tree." >&2
    echo "  Commit, stash, or clean these first:" >&2
    printf '%s\n' "$dirty" | sed 's/^/    /' >&2
    exit 2
fi

# --- resolve the ref -------------------------------------------------------
full_sha="$(git -C "$repo_root" rev-parse --verify --quiet "${ref}^{commit}")" \
    || die "cannot resolve ref to a commit: $ref"
short_sha="$(git -C "$repo_root" rev-parse --short="$SHORT_LEN" "$full_sha")"
commit_date="$(git -C "$repo_root" show -s --format=%cI "$full_sha")"

# --- version ---------------------------------------------------------------
version_source=""
if [[ -n "$version" ]]; then
    version_source="flag"
else
    tag="$(git -C "$repo_root" describe --exact-match --tags "$full_sha" 2>/dev/null || true)"
    if [[ -n "$tag" ]]; then
        version="${tag#v}"
        version_source="tag"
    else
        version="$DEFAULT_VERSION"
        version_source="default"
        echo "$_prog: warning — no tag points at $short_sha; using the ADR-D7 starting" >&2
        echo "  version '$DEFAULT_VERSION'. The artifact records version_source=default so the" >&2
        echo "  label is never mistaken for a tagged release." >&2
    fi
fi
[[ "$version" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] \
    || die "invalid version '$version' (allowed: alphanumerics, dot, dash, underscore)"

stem="${ARTIFACT_NAME}-${version}-${short_sha}"
tarball_name="${stem}.tar.gz"

# --- output directory ------------------------------------------------------
# Absolutised WITHOUT being created: a refused or gate-failed build must leave
# nothing behind, not even an empty directory. It is created just before the
# artifact is written.
[[ -n "$out_dir" ]] || out_dir="$repo_root/dist"
case "$out_dir" in
    /*) ;;
    *)  out_dir="$PWD/$out_dir" ;;
esac
out_dir="${out_dir%/}"

tarball="$out_dir/$tarball_name"
checksum_file="$tarball.sha256"
metadata_file="$out_dir/${stem}.metadata.json"

# --- scratch ---------------------------------------------------------------
work="$(mktemp -d)"
cleanup() { rm -rf "$work"; }
trap cleanup EXIT INT TERM

say "$_prog: ref '$ref' -> $full_sha"
say "  version : $version ($version_source)"
say "  artifact: $tarball_name"

# --- pre-flight gate, against the exported tree ----------------------------
export_dir="$work/export"
mkdir -p "$export_dir"
git -C "$repo_root" archive --format=tar "$full_sha" | tar -x -C "$export_dir"

[[ -x "$export_dir/bin/sciagent" ]] \
    || die "exported tree has no executable bin/sciagent — not a releasable toolkit commit"
[[ -f "$export_dir/tests/run-all.sh" ]] \
    || die "exported tree has no tests/run-all.sh — refusing to release an ungated commit"

say ""
say "$_prog: gate 1/2 — sciagent validate (against the exported tree)"
# `env -u SCIAGENT_TOOLKIT` matters: bin/sciagent honours an inherited
# SCIAGENT_TOOLKIT, so without this the gate could validate the BUILDER's
# toolkit instead of the exported one.
if ! ( cd "$export_dir" && env -u SCIAGENT_TOOLKIT ./bin/sciagent validate ); then
    die "sciagent validate failed for $short_sha — no artifact produced"
fi

say ""
say "$_prog: gate 2/2 — tests/run-all.sh (against the exported tree)"
if ! ( cd "$export_dir" && env -u SCIAGENT_TOOLKIT bash tests/run-all.sh ); then
    die "test suite failed for $short_sha — no artifact produced"
fi

# --- metadata --------------------------------------------------------------
# Deliberately contains NO wall-clock build time, hostname, or builder identity:
# every field is derived from the ref, so the metadata is as reproducible as the
# tarball. `commit_date` (git's own committer date) is the time that matters and
# is a property of the commit, not of the build.
core_meta="$work/meta/$EMBEDDED_METADATA"
mkdir -p "$work/meta"
{
    printf '{\n'
    printf '  "schema": %s,\n' "$METADATA_SCHEMA"
    printf '  "name": "%s",\n'            "$(_json_escape "$ARTIFACT_NAME")"
    printf '  "version": "%s",\n'         "$(_json_escape "$version")"
    printf '  "version_source": "%s",\n'  "$(_json_escape "$version_source")"
    printf '  "ref": "%s",\n'             "$(_json_escape "$ref")"
    printf '  "commit": "%s",\n'          "$full_sha"
    printf '  "commit_short": "%s",\n'    "$short_sha"
    printf '  "commit_date": "%s",\n'     "$(_json_escape "$commit_date")"
    printf '  "artifact": "%s",\n'        "$(_json_escape "$tarball_name")"
    printf '  "archive_prefix": "%s",\n'  "$(_json_escape "$stem/")"
    printf '  "install_subdir": "share/%s/versions/%s",\n' "$ARTIFACT_NAME" "$full_sha"
    printf '  "executable": "bin/sciagent",\n'
    printf '  "build_method": "git archive --format=tar <commit> | gzip -n",\n'
    printf '  "gates": ["sciagent validate", "tests/run-all.sh"]\n'
    printf '}\n'
} > "$core_meta"

# --- the artifact ----------------------------------------------------------
# `--add-file` places the metadata inside the archive with the COMMIT's mtime
# (not the temp file's), verified — so embedding it does not cost determinism.
tmp_tarball="$work/$tarball_name"
git -C "$repo_root" archive --format=tar \
    --prefix="$stem/" \
    --add-file="$core_meta" \
    "$full_sha" \
    | gzip -n > "$tmp_tarball"

digest="$(sha256_of "$tmp_tarball")"
size_bytes="$(wc -c < "$tmp_tarball" | tr -d ' ')"

# sidecar metadata = the embedded core plus the two fields that cannot live
# inside the file they describe.
tmp_meta="$work/${stem}.metadata.json"
sed '$d' "$core_meta" > "$tmp_meta"            # drop the closing brace
{
    printf '  ,"sha256": "%s"\n' "$digest"
    printf '  ,"size_bytes": %s\n' "$size_bytes"
    printf '}\n'
} >> "$tmp_meta"

printf '%s  %s\n' "$digest" "$tarball_name" > "$work/${tarball_name}.sha256"

# Publish all three atomically-ish (each mv is a rename; the artifact lands last
# so a checksum file never exists without its tarball). The output directory is
# created only now, so every refusal above leaves the filesystem untouched.
mkdir -p "$out_dir"
mv "$tmp_meta"                      "$metadata_file"
mv "$work/${tarball_name}.sha256"   "$checksum_file"
mv "$tmp_tarball"                   "$tarball"

say ""
say "$_prog: built"
echo "$tarball"
echo "$checksum_file"
echo "$metadata_file"
say ""
say "  sha256 : $digest"
say "  size   : $size_bytes bytes"
say "  commit : $full_sha"
say ""
say "  verify : (cd $out_dir && sha256sum -c $tarball_name.sha256)"
say "  install: ./install.sh --archive $tarball --checksum $checksum_file --prefix \"\$HOME/.local\""
