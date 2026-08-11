# lib/sciagent/block.sh — AGENTS.md managed-block primitives.
#
# Markers (parametric; id is MANDATORY on every public call — there is no
# default):
#   <!-- BEGIN SCIAGENT:<ID> v1 hash=<sha1> -->
#   <!-- END SCIAGENT:<ID> -->
#
# For id=ROLES this is byte-identical to the original hard-coded markers.
# Three ids are in production use: ROLES (the effective-stack block written
# by activate.sh, checked by status.sh), CRAFT (craft.sh, the owner's
# standing craft conventions), and CONTEXT (provision.sh). All three coexist
# in the same AGENTS.md, each independently drift-checked.
#
# id used to default to ROLES on every function below. That default is
# GONE (2026-08, Phase 5c follow-up): a bare `block_write AGENTS.md "$body"`
# silently meant "write the ROLES block", which is exactly the kind of trap
# that outlives the person who understood it — it briefly looked, from the
# call sites alone, like ROLES might be a dead id nothing writes any more.
# It isn't; the default was just hiding it. Every call site now names its
# id explicitly. Public functions below reject a missing/empty id outright.
#
# Exit conventions:
#   block_read         0=ok, 1=no markers, 2=one marker only
#   block_hash_check   0=match, 1=no block, 2=corrupted, 3=drift
#   block_write/remove always 0 unless I/O fails
#   any public function with a missing/empty id: returns 4, prints to stderr
#
# rc=4 is its own code, distinct from every function's own rc=1 ("no block" /
# "no markers"). A caller bug (forgot the id) must never be mistaken for a
# legitimate "there is no block here" — that collision would let e.g.
# `block_hash_check "$f" "$maybe_empty_var"` take the "safe to write, no
# block yet" branch on what was actually a bug, which is the most
# destructive branch available. See test_block_id_required.sh.

# shellcheck shell=bash

# Marker-string helpers. Both functions live ONLY in block.sh so the
# "managed-block framing must not leak" invariant (test_block_marker_boundary)
# is maintained. Callers outside block.sh always go through the public API.
_block_begin_prefix() {
    # Emit the fixed portion of the BEGIN marker for a given block id.
    # Format: <!-- BEGIN SCIAGENT:<id> v1 hash=
    local id="${1:-ROLES}"
    printf '<!-- BEGIN SCIAGENT:%s v1 hash=' "$id"
}

_block_end_marker() {
    # Emit the full END marker for a given block id.
    # Format: <!-- END SCIAGENT:<id> -->
    local id="${1:-ROLES}"
    printf '<!-- END SCIAGENT:%s -->' "$id"
}

_sha1() {
    sha1sum | awk '{print $1}'
}

# _block_require_id <id> <caller-name>
# Every public function below calls this before touching the filesystem.
# There is no default id (see file header) — a missing/empty id is a caller
# bug, not "assume ROLES". Callers propagate this as rc=4, NOT rc=1 — rc=1
# already means "no block"/"no markers" for block_read/block_hash_check, and
# a caller-bug code must never collide with a legitimate result code.
_block_require_id() {
    local id="$1" caller="$2"
    if [[ -z "$id" ]]; then
        echo "$caller: missing required <id> argument (no default — pass ROLES/CRAFT/CONTEXT explicitly)" >&2
        return 1
    fi
    return 0
}

# _block_replace_preserving_mode <tmp> <target>
# Move <tmp> over <target>, keeping <target>'s permission bits.
#
# mktemp creates 0600, and a bare `mv` carries that mode onto the target — so
# every in-place block rewrite silently stripped group/other read. That already
# happened across the fleet: AGENTS.md sits at 0600 in eleven analysis repos,
# including a shared lab tree where its own siblings are 0644 and colleagues
# consequently cannot read it. AGENTS.md is a file other people are meant to
# read; a managed-block rewrite has no business changing who can.
#
# The mode is captured BEFORE the mv, since the target is replaced. If stat is
# unavailable or the target vanished, fall back to the plain mv rather than
# failing the write — losing the mode is bad, losing the block is worse.
_block_replace_preserving_mode() {
    local tmp="$1" target="$2"
    local mode=""
    [[ -f "$target" ]] && mode=$(stat -c '%a' "$target" 2>/dev/null || true)
    mv "$tmp" "$target" || return 1
    [[ -n "$mode" ]] && chmod "$mode" "$target" 2>/dev/null
    return 0
}

_count_lines() {
    # echo a numeric count of matching lines for a fixed-string pattern.
    local pat="$1" file="$2"
    if [[ ! -f "$file" ]]; then
        echo 0
        return
    fi
    local n
    n=$(grep -cF -- "$pat" "$file" 2>/dev/null) || n=0
    echo "$n"
}

block_read() {
    local file="$1"
    local id="${2:-}"
    _block_require_id "$id" "block_read" || return 4
    [[ -f "$file" ]] || return 1
    local beg_prefix end_marker has_begin has_end
    beg_prefix=$(_block_begin_prefix "$id")
    end_marker=$(_block_end_marker "$id")
    has_begin=$(_count_lines "$beg_prefix" "$file")
    has_end=$(_count_lines "$end_marker" "$file")
    if (( has_begin == 0 && has_end == 0 )); then
        return 1
    fi
    if (( has_begin == 0 || has_end == 0 )); then
        return 2
    fi
    awk -v BEG="$beg_prefix" -v END_MARK="$end_marker" '
        index($0, BEG) == 1 { inblock=1; next }
        $0 == END_MARK     { inblock=0; next }
        inblock { print }
    ' "$file"
    return 0
}

# Public: emit the hash stored in the BEGIN marker (empty if no block).
# Activate/inject reuse this to populate the manifest's BLOCK_HASH after
# block_write, so canonicalisation rules live in exactly one place.
block_stored_hash() {
    local file="$1"
    local id="${2:-}"
    _block_require_id "$id" "block_stored_hash" || return 4
    local beg_prefix end_marker
    beg_prefix=$(_block_begin_prefix "$id")
    end_marker=$(_block_end_marker "$id")
    # Strip prefix and trailing " -->" from the BEGIN line.
    grep -F -- "$beg_prefix" "$file" 2>/dev/null \
        | head -n1 \
        | sed -e "s|^$beg_prefix||" -e 's| -->$||'
}

# block_line_range <file> <id>
# Print "<begin-line> <end-line>" (1-based) for the managed block, or nothing
# (return 1) when there is no complete block. Keeps marker-string knowledge in
# block.sh so consumers (e.g. status.sh) never grep the literals themselves.
block_line_range() {
    local file="$1"
    local id="${2:-}"
    _block_require_id "$id" "block_line_range" || return 4
    [[ -f "$file" ]] || return 1
    local beg_prefix end_marker
    beg_prefix=$(_block_begin_prefix "$id")
    end_marker=$(_block_end_marker "$id")
    grep -qF -- "$beg_prefix" "$file" 2>/dev/null || return 1
    local lb le
    lb=$(grep -nF -- "$beg_prefix" "$file" | head -n1 | cut -d: -f1)
    le=$(grep -nF -- "$end_marker" "$file" | head -n1 | cut -d: -f1)
    [[ -n "$lb" && -n "$le" ]] || return 1
    printf '%s %s\n' "$lb" "$le"
}

block_hash_check() {
    local file="$1"
    local id="${2:-}"
    _block_require_id "$id" "block_hash_check" || return 4
    [[ -f "$file" ]] || return 1
    local beg_prefix end_marker has_begin has_end
    beg_prefix=$(_block_begin_prefix "$id")
    end_marker=$(_block_end_marker "$id")
    has_begin=$(_count_lines "$beg_prefix" "$file")
    has_end=$(_count_lines "$end_marker" "$file")
    if (( has_begin == 0 && has_end == 0 )); then
        return 1
    fi
    if (( has_begin == 0 || has_end == 0 )); then
        return 2
    fi
    local stored actual
    stored=$(block_stored_hash "$file" "$id")
    actual=$(block_read "$file" "$id" | _sha1)
    if [[ "$stored" == "$actual" ]]; then
        return 0
    else
        return 3
    fi
}

# block_write <file> <body> <id>
# Idempotent. Preserves bytes outside markers.
block_write() {
    local file="$1"
    local body="$2"
    local id="${3:-}"
    _block_require_id "$id" "block_write" || return 4
    # INVARIANT: body must be hashed AFTER trailing-newline canonicalisation,
    # so the stored hash matches block_read's awk-based reconstruction (awk
    # `print` always emits a trailing \n).
    [[ "${body: -1}" == $'\n' ]] || body="${body}"$'\n'
    local hash
    hash=$(printf '%s' "$body" | _sha1)
    local beg_prefix end_marker
    beg_prefix=$(_block_begin_prefix "$id")
    end_marker=$(_block_end_marker "$id")
    local begin="${beg_prefix}${hash} -->"

    if [[ ! -f "$file" ]]; then
        printf '%s\n%s%s\n' "$begin" "$body" "$end_marker" > "$file"
        return 0
    fi

    local has_begin has_end
    has_begin=$(_count_lines "$beg_prefix" "$file")
    has_end=$(_count_lines "$end_marker" "$file")
    if (( has_begin > 0 && has_end > 0 )); then
        local tmp
        tmp=$(mktemp)
        awk -v BEG="$beg_prefix" -v END_MARK="$end_marker" \
            -v NEWBEG="$begin" -v BODY="$body" '
            BEGIN { state=0 }
            state==0 && index($0, BEG) == 1 {
                print NEWBEG
                printf "%s", BODY    # BODY already ends with \n (canonicalised)
                print END_MARK
                state=1
                next
            }
            state==1 && $0 == END_MARK { state=2; next }
            state==1 { next }
            { print }
        ' "$file" > "$tmp"
        _block_replace_preserving_mode "$tmp" "$file" || return 1
        return 0
    fi

    # No markers — append fresh block, preceded by one blank line.
    # File must end `...\n\n` before the BEGIN marker.
    if [[ -s "$file" ]]; then
        local last_byte
        last_byte=$(tail -c 1 "$file" | od -An -tx1 | tr -d ' \n')
        [[ "$last_byte" == "0a" ]] || printf '\n' >> "$file"
        printf '\n' >> "$file"
    fi
    printf '%s\n%s%s\n' "$begin" "$body" "$end_marker" >> "$file"
}

# block_remove <file> <id>
# Removes block plus the single blank line immediately preceding BEGIN
# (the separator block_write inserts). Bytes elsewhere unchanged.
block_remove() {
    local file="$1"
    local id="${2:-}"
    _block_require_id "$id" "block_remove" || return 4
    [[ -f "$file" ]] || return 0
    local beg_prefix end_marker
    beg_prefix=$(_block_begin_prefix "$id")
    end_marker=$(_block_end_marker "$id")
    grep -qF -- "$beg_prefix" "$file" || return 0
    local tmp
    tmp=$(mktemp)
    awk -v BEG="$beg_prefix" -v END_MARK="$end_marker" '
        {
            lines[NR] = $0
        }
        END {
            for (i=1; i<=NR; i++) {
                if (index(lines[i], BEG) == 1) begin_n=i
                if (lines[i] == END_MARK)      { end_n=i; break }
            }
            if (!begin_n || !end_n) {
                for (i=1; i<=NR; i++) print lines[i]
                exit
            }
            start = begin_n
            if (start > 1 && lines[start-1] == "") start = begin_n - 1
            for (i=1; i<start; i++) print lines[i]
            for (i=end_n+1; i<=NR; i++) print lines[i]
        }
    ' "$file" > "$tmp"
    _block_replace_preserving_mode "$tmp" "$file"
}
