# lib/sciagent/symlinks.sh — dual-track symlink + manifest helpers.
#
# Manifest format: line-delimited records (no jq dependency). The file name
# keeps the .json extension for forward-compatibility, but its contents are
# explicitly bash-parseable to avoid pulling in jq. State-of-the-art schema:
#
#   VERSION 1
#   STACK <role1> [<role2>]
#   BLOCK_HASH <sha1>
#   SYMLINK <relative-path>
#   SYMLINK <relative-path>
#   ...
#
# Rationale: simple grep+awk parsing; jq is optional and not assumed.

# shellcheck shell=bash

_MANIFEST_PATH=".sciagent/manifest.json"

# symlink_create_dual <category> <name> <canonical_path>
# Categories: skills (dir symlink), agents/commands/output-styles (file symlink).
# Records both .claude/<cat>/<name> and .agents/<cat>/<name>.
symlink_create_dual() {
    local category="$1"
    local name="$2"
    local canonical="$3"

    local claude_path agents_path
    case "$category" in
        skills)
            claude_path=".claude/skills/$name"
            agents_path=".agents/skills/$name"
            ;;
        agents)
            claude_path=".claude/agents/${name}.md"
            agents_path=".agents/agents/${name}.md"
            ;;
        commands)
            claude_path=".claude/commands/${name}.md"
            agents_path=".agents/commands/${name}.md"
            ;;
        output-styles)
            # Claude-only. No .agents/ mirror.
            claude_path=".claude/output-styles/${name}.md"
            agents_path=""
            ;;
        *)
            echo "symlink_create_dual: unknown category '$category'" >&2
            return 1
            ;;
    esac

    mkdir -p "$(dirname "$claude_path")"
    ln -sfn "$canonical" "$claude_path"
    _manifest_record_symlink "$claude_path"

    if [[ -n "$agents_path" ]]; then
        mkdir -p "$(dirname "$agents_path")"
        ln -sfn "$canonical" "$agents_path"
        _manifest_record_symlink "$agents_path"
    fi
}

# Internal: append a SYMLINK line to a staging buffer file. Activate code
# accumulates symlinks and finalises the manifest in one write.
_manifest_staging=""

manifest_begin() {
    local stack="$1"   # space-separated role names
    _manifest_staging=$(mktemp)
    {
        printf 'VERSION 1\n'
        printf 'STACK %s\n' "$stack"
    } > "$_manifest_staging"
}

_manifest_record_symlink() {
    [[ -n "$_manifest_staging" ]] || return 0
    printf 'SYMLINK %s\n' "$1" >> "$_manifest_staging"
}

manifest_finalize() {
    local block_hash="$1"
    [[ -n "$_manifest_staging" ]] || return 1
    printf 'BLOCK_HASH %s\n' "$block_hash" >> "$_manifest_staging"
    mkdir -p "$(dirname "$_MANIFEST_PATH")"
    mv "$_manifest_staging" "$_MANIFEST_PATH"
    _manifest_staging=""
}

# symlink_teardown_all — read manifest, remove only the symlinks we own.
symlink_teardown_all() {
    [[ -f "$_MANIFEST_PATH" ]] || return 0
    local path
    while read -r kind path; do
        [[ "$kind" == "SYMLINK" ]] || continue
        if [[ ! -e "$path" && ! -L "$path" ]]; then
            echo "warning: $path already gone, skipping" >&2
            continue
        fi
        if [[ ! -L "$path" ]]; then
            echo "warning: $path is not a symlink, skipping" >&2
            continue
        fi
        rm "$path"
    done < "$_MANIFEST_PATH"
    rm -f "$_MANIFEST_PATH"
    rmdir .sciagent 2>/dev/null || true
    # Best-effort: clean empty .claude/* and .agents/* dirs we created.
    for d in .claude/skills .claude/agents .claude/commands .claude/output-styles \
             .agents/skills .agents/agents .agents/commands; do
        [[ -d "$d" ]] && rmdir "$d" 2>/dev/null || true
    done
    rmdir .agents 2>/dev/null || true
}

manifest_stack() {
    [[ -f "$_MANIFEST_PATH" ]] || return 1
    grep '^STACK ' "$_MANIFEST_PATH" | head -n1 | cut -d' ' -f2-
}

manifest_exists() {
    [[ -f "$_MANIFEST_PATH" ]]
}
