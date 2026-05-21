# lib/sciagent/block.sh — AGENTS.md managed-block primitives.
#
# Markers (literal):
#   <!-- BEGIN SCIAGENT:ROLES v1 hash=<sha1> -->
#   <!-- END SCIAGENT:ROLES -->
#
# Exit conventions:
#   block_read         0=ok, 1=no markers, 2=one marker only
#   block_hash_check   0=match, 1=no block, 2=corrupted, 3=drift
#   block_write/remove always 0 unless I/O fails

# shellcheck shell=bash

# Plain substrings we search for. grep -F avoids regex-escape pain.
_BLOCK_BEGIN_PREFIX='<!-- BEGIN SCIAGENT:ROLES v1 hash='
_BLOCK_END='<!-- END SCIAGENT:ROLES -->'

_sha1() {
    sha1sum | awk '{print $1}'
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
    [[ -f "$file" ]] || return 1
    local has_begin has_end
    has_begin=$(_count_lines "$_BLOCK_BEGIN_PREFIX" "$file")
    has_end=$(_count_lines "$_BLOCK_END" "$file")
    if (( has_begin == 0 && has_end == 0 )); then
        return 1
    fi
    if (( has_begin == 0 || has_end == 0 )); then
        return 2
    fi
    awk -v BEG="$_BLOCK_BEGIN_PREFIX" -v END_MARK="$_BLOCK_END" '
        index($0, BEG) == 1 { inblock=1; next }
        $0 == END_MARK     { inblock=0; next }
        inblock { print }
    ' "$file"
    return 0
}

_block_stored_hash() {
    local file="$1"
    # Strip prefix and trailing " -->" from the BEGIN line.
    grep -F -- "$_BLOCK_BEGIN_PREFIX" "$file" 2>/dev/null \
        | head -n1 \
        | sed -e "s|^$_BLOCK_BEGIN_PREFIX||" -e 's| -->$||'
}

block_hash_check() {
    local file="$1"
    [[ -f "$file" ]] || return 1
    local has_begin has_end
    has_begin=$(_count_lines "$_BLOCK_BEGIN_PREFIX" "$file")
    has_end=$(_count_lines "$_BLOCK_END" "$file")
    if (( has_begin == 0 && has_end == 0 )); then
        return 1
    fi
    if (( has_begin == 0 || has_end == 0 )); then
        return 2
    fi
    local stored actual
    stored=$(_block_stored_hash "$file")
    actual=$(block_read "$file" | _sha1)
    if [[ "$stored" == "$actual" ]]; then
        return 0
    else
        return 3
    fi
}

# block_write <file> <body>
# Idempotent. Preserves bytes outside markers.
block_write() {
    local file="$1"
    local body="$2"
    # Canonicalise: body always ends with exactly one newline before hashing,
    # so the stored hash matches block_read's awk-based reconstruction.
    [[ "${body: -1}" == $'\n' ]] || body="${body}"$'\n'
    local hash
    hash=$(printf '%s' "$body" | _sha1)
    local begin="${_BLOCK_BEGIN_PREFIX}${hash} -->"
    local end="$_BLOCK_END"

    if [[ ! -f "$file" ]]; then
        { printf '%s\n' "$begin"
          printf '%s' "$body"
          [[ "${body: -1}" == $'\n' ]] || printf '\n'
          printf '%s\n' "$end"
        } > "$file"
        return 0
    fi

    local has_begin has_end
    has_begin=$(_count_lines "$_BLOCK_BEGIN_PREFIX" "$file")
    has_end=$(_count_lines "$_BLOCK_END" "$file")
    if (( has_begin > 0 && has_end > 0 )); then
        local tmp
        tmp=$(mktemp)
        awk -v BEG="$_BLOCK_BEGIN_PREFIX" -v END_MARK="$_BLOCK_END" \
            -v NEWBEG="$begin" -v BODY="$body" '
            BEGIN { state=0 }
            state==0 && index($0, BEG) == 1 {
                print NEWBEG
                printf "%s", BODY
                if (substr(BODY, length(BODY), 1) != "\n") print ""
                print END_MARK
                state=1
                next
            }
            state==1 && $0 == END_MARK { state=2; next }
            state==1 { next }
            { print }
        ' "$file" > "$tmp"
        mv "$tmp" "$file"
        return 0
    fi

    # No markers — append fresh block, preceded by one blank line.
    # We want file to end with `...\n\n` before the BEGIN marker.
    if [[ -s "$file" ]]; then
        local last_byte
        # od strips trailing-newline pitfalls of $(...)
        last_byte=$(tail -c 1 "$file" | od -An -tx1 | tr -d ' \n')
        if [[ "$last_byte" != "0a" ]]; then
            printf '\n' >> "$file"
        fi
        printf '\n' >> "$file"
    fi
    {
        printf '%s\n' "$begin"
        printf '%s' "$body"
        [[ "${body: -1}" == $'\n' ]] || printf '\n'
        printf '%s\n' "$end"
    } >> "$file"
}

# block_remove <file>
# Removes block plus the single blank line immediately preceding BEGIN
# (the separator block_write inserts). Bytes elsewhere unchanged.
block_remove() {
    local file="$1"
    [[ -f "$file" ]] || return 0
    grep -qF -- "$_BLOCK_BEGIN_PREFIX" "$file" || return 0
    local tmp
    tmp=$(mktemp)
    awk -v BEG="$_BLOCK_BEGIN_PREFIX" -v END_MARK="$_BLOCK_END" '
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
    mv "$tmp" "$file"
}
