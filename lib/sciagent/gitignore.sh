# lib/sciagent/gitignore.sh — sciagent gitignore [<path>]
# Manages a SCIAGENT:GITIGNORE block inside a .gitignore file.

# shellcheck shell=bash

_SCIAGENT_GITIGNORE_BEGIN='# BEGIN SCIAGENT:GITIGNORE'
_SCIAGENT_GITIGNORE_END='# END SCIAGENT:GITIGNORE'

_sciagent_gitignore_block() {
    cat <<'EOF'
# BEGIN SCIAGENT:GITIGNORE
docs/_internal/
.claude/
.agents/
.gemini/
.sciagent/
.mcp.json
.env
.env.*
EOF
    printf '%s\n' '# END SCIAGENT:GITIGNORE'
}

cmd_gitignore() {
    local target="${1:-$(pwd)/.gitignore}"

    # Create the file if it doesn't exist.
    if [[ ! -e "$target" ]]; then
        touch "$target"
    fi

    if grep -qF "$_SCIAGENT_GITIGNORE_BEGIN" "$target"; then
        # Replace existing block in-place using awk.
        local tmpfile
        tmpfile=$(mktemp)
        awk -v begin="$_SCIAGENT_GITIGNORE_BEGIN" \
            -v end="$_SCIAGENT_GITIGNORE_END" \
            '
            $0 == begin { in_block=1; next }
            $0 == end   { in_block=0; next }
            !in_block   { print }
            ' "$target" > "$tmpfile"

        # Build replacement file: original (minus old block) + new block.
        {
            cat "$tmpfile"
            _sciagent_gitignore_block
        } > "$target.tmp"
        mv "$target.tmp" "$target"
        rm -f "$tmpfile"
        echo "updated SCIAGENT:GITIGNORE block in: $target"
    else
        # Append the block.
        {
            # Ensure there is a newline before the block if file is non-empty.
            if [[ -s "$target" ]]; then
                printf '\n'
            fi
            _sciagent_gitignore_block
        } >> "$target"
        echo "appended SCIAGENT:GITIGNORE block to: $target"
    fi
}
