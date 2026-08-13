# lib/scio/ownership.sh — hash-and-cede discipline for materialized bodies.

# shellcheck shell=bash

_SCIO_TEMPLATE_PROVENANCE="templates/PROVENANCE.sha1"

# True when the body matches a version the toolkit has shipped.
ownership_template_hash_known() {
    local hash="$1" rel="$2"
    local manifest="$SCIO_TOOLKIT/$_SCIO_TEMPLATE_PROVENANCE"
    [[ -n "$hash" && -f "$manifest" ]] || return 1
    local known_hash known_path
    while read -r known_hash known_path; do
        [[ "$known_hash" == \#* || -z "$known_hash" ]] && continue
        [[ "$known_hash" == "$hash" && "$known_path" == "$rel" ]] && return 0
    done < "$manifest"
    return 1
}

# Materialize a body, refresh toolkit-owned bytes, and cede user-edited bytes.
ownership_ensure_body() {
    local src="$1" dst="$2" rel="$3" state="$4" ceded="$5" mode="${6:-exec}"
    [[ -f "$src" ]] || return 0

    mkdir -p "$(dirname "$dst")"

    if [[ ! -f "$dst" ]]; then
        cp "$src" "$dst" || { echo "scio: failed to write $dst" >&2; return 1; }
        [[ "$mode" == exec ]] && chmod +x "$dst"
        echo "wrote: $dst"
        mkdir -p "$(dirname "$state")"
        scio_sha1_file "$dst" > "$state"
        rm -f "$ceded"
        return 0
    fi

    [[ "$mode" == exec ]] && { chmod +x "$dst" 2>/dev/null || true; }

    local current template
    current=$(scio_sha1_file "$dst")
    template=$(scio_sha1_file "$src")

    if [[ "$current" == "$template" ]]; then
        if [[ ! -f "$state" && ! -f "$ceded" ]]; then
            mkdir -p "$(dirname "$state")"
            printf '%s\n' "$current" > "$state"
        fi
        return 0
    fi

    if [[ -f "$ceded" && "$(<"$ceded")" == "$current" ]]; then
        return 0
    fi

    local ours=false
    if [[ -f "$state" && "$(<"$state")" == "$current" ]]; then
        ours=true
    elif ownership_template_hash_known "$current" "$rel"; then
        ours=true
    fi

    if [[ "$ours" == true ]]; then
        cp "$src" "$dst" || { echo "scio: failed to refresh $dst" >&2; return 1; }
        [[ "$mode" == exec ]] && chmod +x "$dst"
        echo "refreshed: $dst (was an older toolkit version)"
        mkdir -p "$(dirname "$state")"
        printf '%s\n' "$template" > "$state"
        rm -f "$ceded"
        return 0
    fi

    echo "scio: $dst differs from the toolkit's version and was not written by scio" >&2
    echo "  leaving it as yours; scio will not manage it from now on" >&2
    rm -f "$state"
    mkdir -p "$(dirname "$ceded")"
    printf '%s\n' "$current" > "$ceded"
}
