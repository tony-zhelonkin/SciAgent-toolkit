#!/usr/bin/env bash
# install.sh — place a `scio` release on THIS machine, from a local file.
#
# Implements `docs/proposals/2026-08-11-offline-distribution/10_packaging_contracts.md` §2
# and ADR-D3. It puts a program on the machine. It does not decide what any
# project gets — that is `scio`'s job, and merging the two is the exact
# mistake `00_INDEX.md` §2 names.
#
# ---------------------------------------------------------------------------
# Usage
# ---------------------------------------------------------------------------
#   ./install.sh --archive  scio-0.1.0-<sha>.tar.gz \
#                --checksum scio-0.1.0-<sha>.tar.gz.sha256 \
#                --prefix   "$HOME/.local"
#
#   ./install.sh --list                     [--prefix DIR]
#   ./install.sh --uninstall <commit-sha>   [--prefix DIR] [--dry-run]
#
# ---------------------------------------------------------------------------
# What this script cannot do, by construction
# ---------------------------------------------------------------------------
#   * No URL is accepted anywhere. `--archive` and `--checksum` must name
#     existing regular files; a URL-shaped argument is rejected with an
#     explanation, not fetched. There is no code path that opens a socket.
#   * No curl, wget, npm, registry, telemetry, or update check. Grep this file
#     for any of those words and the invariant is checkable rather than
#     aspirational (tests/test_install_invariants.sh does exactly that).
#   * No harness detection, no settings.json, no AGENTS.md, no skill mounting,
#     no project mutation of any kind. Grep for a harness name; there is none.
#
# ---------------------------------------------------------------------------
# Layout it writes (everything under --prefix, nothing outside it)
# ---------------------------------------------------------------------------
#   <prefix>/share/scio/versions/<full-git-sha>/     the extracted release
#   <prefix>/share/scio/receipts/<full-git-sha>.json the installation receipt
#   <prefix>/bin/scio -> ../share/scio/versions/<full-git-sha>/bin/scio
#
# Content-addressed on the FULL 40-char commit SHA (ADR-D1, and §1's "why the
# full SHA in metadata"): two releases coexist by construction because their
# directories cannot collide. The SHA is read from `.scio-release.json` inside
# the archive, never from the filename's abbreviation.
#
# Owner ruling 7 aligns the installed command with the `scio` artifact name.
# ADR-D7's original artifact-only decision remains in the decision record with
# a dated rebrand note.
#
# ---------------------------------------------------------------------------
# Fleet precedence is preserved, and that is why this links rather than wraps
# ---------------------------------------------------------------------------
# `bin/scio` resolves its own real path to derive `$SCIO_TOOLKIT`, so a
# symlink here resolves to the version directory and nothing else. This script
# writes no shell profile, exports no variable, and installs no wrapper that
# could pin `$SCIO_TOOLKIT` globally. A project that ships its own
# `01_modules/SciAgent-toolkit` therefore still wins: `scio link` refuses
# to bind such a project from this global copy. Exact-commit
# reproducibility for the fleet depends on that refusal, so nothing installed
# here may weaken it.
#
# ---------------------------------------------------------------------------
# Atomicity and reversibility
# ---------------------------------------------------------------------------
# The archive is verified against its checksum BEFORE anything is extracted.
# Extraction lands in a staging directory under the same filesystem as the final
# location and is moved into place with a single rename, so an interrupted or
# failed install leaves no half-tree — the staging directory is removed by an
# EXIT/INT/TERM trap. Every artifact is recorded in a receipt, and `--uninstall`
# removes only what the receipt proves this script wrote: a symlink is removed
# only if it still points where the receipt says it was pointed. A link that now
# serves another version installed under this prefix is ordinary coexistence —
# that version's receipt owns it, so it stays and the status stays clean. Only an
# unprovable link is left alone and reported (exit 3), the same "cannot verify →
# do not touch" discipline the toolkit's teardown uses.
set -uo pipefail

ARTIFACT_NAME="scio"
EMBEDDED_METADATA=".scio-release.json"
EXECUTABLE="scio"
RECEIPT_SCHEMA=1

_prog="${0##*/}"

die() { echo "$_prog: $*" >&2; exit 1; }
say() { [[ -n "${QUIET:-}" ]] || echo "$*"; }
act() { if [[ -n "${DRY_RUN:-}" ]]; then echo "[dry-run] would $*"; else echo "$*"; fi; }

usage() {
    cat <<EOF
$_prog — install a local ${ARTIFACT_NAME} release tarball on this machine

Install:
  $_prog --archive FILE --checksum FILE [--prefix DIR] [--dry-run] [--force]

Inspect / remove:
  $_prog --list [--prefix DIR]
  $_prog --uninstall <commit-sha> [--prefix DIR] [--dry-run]

Options:
  --archive FILE     Local release tarball. A path, never a URL.
  --checksum FILE    Its .sha256 sidecar (sha256sum format). Verified first.
  --prefix DIR       Install prefix (default: \$HOME/.local).
  --force            Re-extract an already-installed commit; replace a
                     <prefix>/bin/$EXECUTABLE this script does not own.
  --dry-run          Print every intended action; write nothing at all.
  --quiet            Less chatter.
  -h, --help         This message.

Writes only under DIR:
  DIR/share/${ARTIFACT_NAME}/versions/<full-sha>/   the release
  DIR/share/${ARTIFACT_NAME}/receipts/<full-sha>.json  the uninstall receipt
  DIR/bin/$EXECUTABLE -> ../share/${ARTIFACT_NAME}/versions/<full-sha>/bin/$EXECUTABLE

No network access. No project is touched: this installs a program; \`$EXECUTABLE\`
decides what a project gets.
EOF
}

# --- small helpers ---------------------------------------------------------

sha256_of() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | awk '{print $1}'
    elif command -v shasum >/dev/null 2>&1; then
        shasum -a 256 "$1" | awk '{print $1}'
    else
        die "neither sha256sum nor shasum found; cannot verify the archive"
    fi
}

# _json_escape <string>
_json_escape() {
    local s="$1"
    s="${s//\\/\\\\}"; s="${s//\"/\\\"}"
    s="${s//$'\n'/\\n}"; s="${s//$'\r'/\\r}"; s="${s//$'\t'/\\t}"
    printf '%s' "$s"
}

# _json_str_field <file> <key> — read a flat string field. jq when present,
# otherwise a targeted grep/sed against this well-known flat schema. It is not a JSON
# parser and does not pretend to be one).
_json_str_field() {
    local file="$1" key="$2"
    if command -v jq >/dev/null 2>&1; then
        jq -r --arg k "$key" '.[$k] // empty' "$file" 2>/dev/null
    else
        grep -o "\"$key\"[[:space:]]*:[[:space:]]*\"[^\"]*\"" "$file" 2>/dev/null \
            | head -1 | sed 's/.*:[[:space:]]*"//; s/"$//'
    fi
}

# _reject_url <label> <value> — there is no fetch path; say so plainly.
_reject_url() {
    local label="$1" value="$2"
    if [[ "$value" =~ ^[A-Za-z][A-Za-z0-9+.-]*:// || "$value" =~ ^git@ ]]; then
        echo "$_prog: $label looks like a URL: $value" >&2
        echo "  This installer has no network code path — not a disabled one, none." >&2
        echo "  Fetch the release yourself (git bundle, scp, a download, an archival" >&2
        echo "  deposit), then pass the local file." >&2
        exit 1
    fi
}

# --- argument parsing ------------------------------------------------------
archive=""
checksum=""
prefix=""
mode="install"
uninstall_target=""
DRY_RUN=""
FORCE=""
QUIET=""

[[ $# -gt 0 ]] || { usage; exit 1; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help) usage; exit 0 ;;
        --archive)   [[ $# -ge 2 ]] || die "--archive requires a path";   archive="$2"; shift 2 ;;
        --checksum)  [[ $# -ge 2 ]] || die "--checksum requires a path";  checksum="$2"; shift 2 ;;
        --prefix)    [[ $# -ge 2 ]] || die "--prefix requires a path";    prefix="$2"; shift 2 ;;
        --uninstall) [[ $# -ge 2 ]] || die "--uninstall requires a commit sha"
                     mode="uninstall"; uninstall_target="$2"; shift 2 ;;
        --list)      mode="list"; shift ;;
        --dry-run)   DRY_RUN=1; shift ;;
        --force)     FORCE=1; shift ;;
        --quiet)     QUIET=1; shift ;;
        --url|--from-url|--fetch|--download)
            # Present ONLY to fail loudly and name the design rule. Removing
            # these cases would make the same input an "unknown option", which
            # reads like an oversight rather than a decision.
            echo "$_prog: $1 does not exist and will not." >&2
            echo "  Transport is not this script's layer (ADR-D2/D3): it places bytes" >&2
            echo "  that are already on the machine. Pass --archive with a local path." >&2
            exit 1 ;;
        -*) die "unknown option: $1 (see --help)" ;;
        *)  die "unexpected argument: $1 (all inputs are named options; see --help)" ;;
    esac
done

# --- prefix ----------------------------------------------------------------
if [[ -z "$prefix" ]]; then
    [[ -n "${HOME:-}" ]] || die "no --prefix given and \$HOME is unset"
    prefix="$HOME/.local"
fi
_reject_url "--prefix" "$prefix"
# Absolutise without creating anything (dry-run must write nothing, including
# the prefix itself).
case "$prefix" in
    /*) ;;
    *)  prefix="$PWD/$prefix" ;;
esac
prefix="${prefix%/}"

share_dir="$prefix/share/$ARTIFACT_NAME"
versions_dir="$share_dir/versions"
receipts_dir="$share_dir/receipts"
bin_dir="$prefix/bin"
bin_link="$bin_dir/$EXECUTABLE"

# ---------------------------------------------------------------------------
# mode: list
# ---------------------------------------------------------------------------
if [[ "$mode" == "list" ]]; then
    if [[ ! -d "$receipts_dir" ]]; then
        echo "no $ARTIFACT_NAME installations under $prefix"
        exit 0
    fi
    current=""
    [[ -L "$bin_link" ]] && current="$(readlink "$bin_link")"
    found=0
    for r in "$receipts_dir"/*.json; do
        [[ -e "$r" ]] || continue
        found=1
        c="$(_json_str_field "$r" commit)"
        v="$(_json_str_field "$r" version)"
        d="$(_json_str_field "$r" version_dir)"
        marker="  "
        [[ -n "$current" && "$current" == *"$c"* ]] && marker="* "
        state="ok"
        [[ -d "$d" ]] || state="MISSING TREE"
        printf '%s%s  %-10s  %s  (%s)\n' "$marker" "$c" "$v" "$d" "$state"
    done
    (( found )) || echo "no $ARTIFACT_NAME installations under $prefix"
    [[ -n "$current" ]] && echo "" && echo "$bin_link -> $current"
    exit 0
fi

# ---------------------------------------------------------------------------
# mode: uninstall — driven entirely by the receipt
# ---------------------------------------------------------------------------
if [[ "$mode" == "uninstall" ]]; then
    [[ "$uninstall_target" =~ ^[0-9a-fA-F]{4,40}$ ]] \
        || die "--uninstall wants a commit sha (4-40 hex chars), got: $uninstall_target"

    receipt="$receipts_dir/$uninstall_target.json"
    if [[ ! -f "$receipt" ]]; then
        matches=()
        for r in "$receipts_dir/$uninstall_target"*.json; do
            [[ -e "$r" ]] && matches+=("$r")
        done
        case ${#matches[@]} in
            0) die "no receipt for '$uninstall_target' under $receipts_dir — refusing to guess what to delete" ;;
            1) receipt="${matches[0]}" ;;
            *) die "'$uninstall_target' matches ${#matches[@]} receipts; give more of the sha" ;;
        esac
    fi

    commit="$(_json_str_field "$receipt" commit)"
    version_dir="$(_json_str_field "$receipt" version_dir)"
    link_path="$(_json_str_field "$receipt" bin_link)"
    link_target="$(_json_str_field "$receipt" bin_link_target)"
    [[ -n "$commit" && -n "$version_dir" ]] \
        || die "receipt $receipt is unreadable (no commit/version_dir) — refusing to delete on a guess"

    # Containment guard: never delete outside the prefix we were told about,
    # whatever a hand-edited receipt claims.
    case "$version_dir" in
        "$versions_dir"/*) ;;
        *) die "receipt's version_dir ($version_dir) is not under $versions_dir — refusing" ;;
    esac

    rc=0
    say "$_prog: uninstalling $commit from $prefix"

    # 1. the bin symlink — removed only if it still points where we put it.
    if [[ -n "$link_path" ]]; then
        if [[ -L "$link_path" ]]; then
            actual="$(readlink "$link_path")"
            other=""
            if [[ "$actual" != "$link_target" ]]; then
                # A link serving another version installed under this prefix is
                # ordinary housekeeping, not an anomaly: that version's receipt
                # owns it, so it stays and the exit status stays clean. Proven by
                # the receipt rather than by the path's shape, so a half-torn-down
                # tree still reads as unprovable.
                candidate="${actual#../share/$ARTIFACT_NAME/versions/}"
                candidate="${candidate%/bin/$EXECUTABLE}"
                if [[ "$candidate" =~ ^[0-9a-f]{40}$ && -f "$receipts_dir/$candidate.json" ]]; then
                    other="$candidate"
                fi
            fi
            if [[ "$actual" == "$link_target" ]]; then
                act "remove symlink $link_path"
                [[ -n "$DRY_RUN" ]] || rm -f "$link_path"
            elif [[ -n "$other" ]]; then
                say "  ($link_path serves ${other:0:12}, still installed here — left in place)"
            else
                echo "$_prog: $link_path points at '$actual', not the recorded" >&2
                echo "  '$link_target' — this installer cannot prove it owns it; left alone." >&2
                rc=3
            fi
        elif [[ -e "$link_path" ]]; then
            echo "$_prog: $link_path exists but is not a symlink — left alone." >&2
            rc=3
        fi
    fi

    # 2. the version tree
    if [[ -d "$version_dir" ]]; then
        act "remove tree $version_dir"
        [[ -n "$DRY_RUN" ]] || rm -rf "$version_dir"
    else
        say "  (tree already absent: $version_dir)"
    fi

    # 3. the receipt last, so an interrupted uninstall stays resumable
    act "remove receipt $receipt"
    [[ -n "$DRY_RUN" ]] || rm -f "$receipt"

    if [[ -z "$DRY_RUN" ]]; then
        remaining=()
        for r in "$receipts_dir"/*.json; do [[ -e "$r" ]] && remaining+=("$r"); done
        if (( ${#remaining[@]} > 0 )); then
            say ""
            say "  still installed: ${#remaining[@]} version(s) — see $_prog --list"
            if [[ ! -e "$bin_link" ]]; then
                say "  $bin_link was removed; re-run the installer for another version to restore it."
            fi
        fi
    fi
    exit $rc
fi

# ---------------------------------------------------------------------------
# mode: install
# ---------------------------------------------------------------------------
[[ -n "$archive" ]]  || die "--archive is required (a local file; see --help)"
[[ -n "$checksum" ]] || die "--checksum is required — an unverified archive is not installable"

_reject_url "--archive"  "$archive"
_reject_url "--checksum" "$checksum"

[[ -f "$archive" ]]  || die "no such file: $archive"
[[ -f "$checksum" ]] || die "no such file: $checksum"

archive="$(cd "$(dirname "$archive")" && pwd)/$(basename "$archive")"
checksum="$(cd "$(dirname "$checksum")" && pwd)/$(basename "$checksum")"

# --- 1. verify the checksum BEFORE touching anything -----------------------
expected="$(awk 'NR==1{print $1}' "$checksum")"
[[ "$expected" =~ ^[0-9a-fA-F]{64}$ ]] \
    || die "checksum file $checksum does not start with a sha256 digest"
actual="$(sha256_of "$archive")"
if [[ "${expected,,}" != "${actual,,}" ]]; then
    echo "$_prog: CHECKSUM MISMATCH — nothing was extracted." >&2
    echo "  archive : $archive" >&2
    echo "  expected: $expected  (from $checksum)" >&2
    echo "  actual  : $actual" >&2
    echo "  The archive is corrupt or is not the one this checksum describes." >&2
    exit 1
fi
say "$_prog: checksum ok ($actual)"

# --- 2. read the release identity out of the archive -----------------------
# Read, do not extract: `tar -O` streams one member to stdout and writes nothing
# to disk, so --dry-run stays a pure read even here.
top="$(tar tzf "$archive" 2>/dev/null | head -1 | cut -d/ -f1)"
[[ -n "$top" ]] || die "cannot list $archive — is it a gzip tarball?"
meta_member="$top/$EMBEDDED_METADATA"
meta_tmp=""
if ! tar xzOf "$archive" "$meta_member" > /dev/null 2>&1; then
    die "archive has no $meta_member — not a $ARTIFACT_NAME release (build it with scripts/build-release.sh)"
fi
meta_tmp="$(mktemp)"
tar xzOf "$archive" "$meta_member" > "$meta_tmp" 2>/dev/null

commit="$(_json_str_field "$meta_tmp" commit)"
version="$(_json_str_field "$meta_tmp" version)"
rel_name="$(_json_str_field "$meta_tmp" name)"
rm -f "$meta_tmp"

[[ "$commit" =~ ^[0-9a-f]{40}$ ]] \
    || die "release metadata carries no full 40-char commit sha (got '$commit')"
[[ -n "$version" ]] || version="unknown"
[[ -n "$rel_name" ]] || rel_name="$ARTIFACT_NAME"

version_dir="$versions_dir/$commit"
receipt="$receipts_dir/$commit.json"
link_target="../share/$ARTIFACT_NAME/versions/$commit/bin/$EXECUTABLE"

say "  release : $rel_name $version"
say "  commit  : $commit"
say "  prefix  : $prefix"

# --- 3. refuse to clobber a bin/scio we do not own -------------------------
# Checked before any write so the refusal cannot leave a half-install behind.
if [[ -e "$bin_link" || -L "$bin_link" ]]; then
    if [[ -L "$bin_link" ]]; then
        cur="$(readlink "$bin_link")"
        if [[ "$cur" != ../share/"$ARTIFACT_NAME"/versions/* && -z "$FORCE" ]]; then
            echo "$_prog: $bin_link is a symlink to '$cur', which this installer did not create." >&2
            echo "  Refusing to replace another tool's link. Re-run with --force to take it over" >&2
            echo "  (the previous target is recorded in the receipt so it can be restored)." >&2
            exit 1
        fi
    elif [[ -d "$bin_link" ]]; then
        die "$bin_link is a directory — refusing to touch it"
    else
        [[ -n "$FORCE" ]] || die "$bin_link exists and is a regular file (another program named $EXECUTABLE?). Re-run with --force to replace it."
    fi
fi
previous_target=""
[[ -L "$bin_link" ]] && previous_target="$(readlink "$bin_link")"

# --- 4. dry-run stops here, having written nothing -------------------------
if [[ -n "$DRY_RUN" ]]; then
    echo "[dry-run] would create directory $versions_dir"
    echo "[dry-run] would extract $(basename "$archive") to $version_dir"
    echo "[dry-run] would write receipt $receipt"
    echo "[dry-run] would link $bin_link -> $link_target"
    [[ -d "$version_dir" ]] && echo "[dry-run] note: $version_dir already exists (would be kept; --force re-extracts)"
    echo "[dry-run] nothing was written"
    exit 0
fi

# --- 5. extract into staging, then move atomically -------------------------
mkdir -p "$versions_dir" "$receipts_dir" "$share_dir/.staging" \
    || die "cannot create $share_dir (check --prefix permissions)"

staging="$(mktemp -d "$share_dir/.staging/install.XXXXXX")" \
    || die "cannot create a staging directory under $share_dir/.staging"
cleanup() {
    # An interrupted or failed install leaves no half-tree: the only thing on
    # disk before the final rename is this staging directory.
    [[ -n "${staging:-}" && -d "$staging" ]] && rm -rf "$staging"
    rmdir "$share_dir/.staging" 2>/dev/null || true
}
trap cleanup EXIT INT TERM

if ! tar xzf "$archive" -C "$staging"; then
    # Covers every way extraction can die part-way — truncated stream, full
    # disk, a resource limit, a signal. The trap has already removed the
    # partial tree by the time this message is read.
    die "extraction failed part-way (truncated archive, or no space / limit reached) — nothing was installed"
fi

src="$staging/$top"
[[ -d "$src" ]] || die "archive does not contain the expected top-level directory '$top'"
[[ -x "$src/bin/$EXECUTABLE" ]] \
    || die "archive has no executable bin/$EXECUTABLE — refusing to install it"

if [[ -d "$version_dir" ]]; then
    if [[ -n "$FORCE" ]]; then
        old="$staging/.replaced"
        mv "$version_dir" "$old" || die "cannot move aside the existing $version_dir"
        mv "$src" "$version_dir" || { mv "$old" "$version_dir"; die "install failed; previous tree restored"; }
        rm -rf "$old"
        say "  replaced the existing tree for $commit (--force)"
    else
        # Content-addressed: the same commit is the same bytes. Converge on the
        # link and the receipt instead of failing.
        say "  $version_dir already present — keeping it (use --force to re-extract)"
    fi
else
    mv "$src" "$version_dir" || die "cannot move the extracted tree into $version_dir"
fi

# --- 6. the receipt (atomic write) -----------------------------------------
installed_at="$(date -u +%Y-%m-%dT%H:%M:%SZ 2>/dev/null || echo unknown)"
archive_sha="$actual"
archive_size="$(wc -c < "$archive" | tr -d ' ')"

receipt_tmp="$staging/receipt.json"
{
    printf '{\n'
    printf '  "schema": %s,\n' "$RECEIPT_SCHEMA"
    printf '  "name": "%s",\n'            "$(_json_escape "$rel_name")"
    printf '  "version": "%s",\n'         "$(_json_escape "$version")"
    printf '  "commit": "%s",\n'          "$commit"
    printf '  "prefix": "%s",\n'          "$(_json_escape "$prefix")"
    printf '  "version_dir": "%s",\n'     "$(_json_escape "$version_dir")"
    printf '  "receipt": "%s",\n'         "$(_json_escape "$receipt")"
    printf '  "bin_link": "%s",\n'        "$(_json_escape "$bin_link")"
    printf '  "bin_link_target": "%s",\n' "$(_json_escape "$link_target")"
    printf '  "replaced_bin_link_target": "%s",\n' "$(_json_escape "$previous_target")"
    printf '  "archive": "%s",\n'         "$(_json_escape "$archive")"
    printf '  "archive_sha256": "%s",\n'  "$archive_sha"
    printf '  "archive_size_bytes": %s,\n' "$archive_size"
    printf '  "installed_at": "%s",\n'    "$installed_at"
    printf '  "installer": "%s"\n'        "$_prog"
    printf '}\n'
} > "$receipt_tmp"
mv "$receipt_tmp" "$receipt"

# --- 7. link only the executable -------------------------------------------
mkdir -p "$bin_dir"
link_tmp="$bin_dir/.$EXECUTABLE.install.$$"
ln -s "$link_target" "$link_tmp"
mv -f "$link_tmp" "$bin_link"    # rename over the old link: atomic swap

say ""
say "$_prog: installed $rel_name $version"
say "  tree    : $version_dir"
say "  receipt : $receipt"
say "  command : $bin_link -> $link_target"
say ""
say "  If $bin_dir is on your PATH, \`$EXECUTABLE --help\` works now."
say "  A project that ships its own toolkit still overrides this copy — that is"
say "  deliberate, and it is what keeps a published analysis pinned to its commit."
say "  Remove with: $_prog --uninstall $commit --prefix $prefix"
