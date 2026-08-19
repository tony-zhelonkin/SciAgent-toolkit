# lib/scio/catalog.sh — toolkit catalog checks used by `scio lint`.

# shellcheck shell=bash

: "${SCIO_DESC_MAX:=350}"

_catalog_root() {
    if [[ -n "${SCIO_TOOLKIT:-}" ]]; then
        printf '%s\n' "$SCIO_TOOLKIT"
        return
    fi
    local self_dir
    self_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    (cd "$self_dir/../.." && pwd)
}

# Print a description value with YAML block-scalar folding and chomping.
_catalog_description() {
    awk '
        function leadspaces(s) {
            n = 0
            while (n < length(s) && substr(s, n+1, 1) == " ") n++
            return n
        }
        function finalize() {
            out = ""
            last = nlines
            while (last > 0 && blank[last]) last--
            if (style == "|") {
                for (i=1; i<=last; i++) out = out (i>1 ? "\n" : "") content[i]
                if (chomp == "clip" && last > 0) out = out "\n"
                else if (chomp == "keep") {
                    for (i=last+1; i<=nlines; i++) out = out "\n"
                    if (last > 0 || nlines > 0) out = out "\n"
                }
            } else {
                blankrun = 0
                first = 1
                for (i=1; i<=last; i++) {
                    if (blank[i]) { blankrun++; continue }
                    if (!first) {
                        if (blankrun > 0) for (b=0; b<blankrun; b++) out = out "\n"
                        else out = out " "
                    }
                    out = out content[i]
                    first = 0
                    blankrun = 0
                }
                if (chomp == "clip" && last > 0) out = out "\n"
                else if (chomp == "keep") {
                    for (i=last+1; i<=nlines; i++) out = out "\n"
                    if (last > 0 || nlines > 0) out = out "\n"
                }
            }
            printf "%s", out
            state = 2
        }
        NR==1 && /^---[ \t]*$/ { infm=1; next }
        infm && /^---[ \t]*$/ { if (state == 1) finalize(); exit }
        !infm { next }
        state==0 && /^description:[ \t]*/ {
            found = 1
            key_indent = leadspaces($0)
            line = $0
            sub(/^description:[ \t]*/, "", line)
            sub(/[ \t]+$/, "", line)
            if (line ~ /^[|>][+-]?[0-9]*[ \t]*(#.*)?$/) {
                style = substr(line, 1, 1)
                rest = substr(line, 2)
                chomp = "clip"
                if (rest ~ /^-/) chomp = "strip"
                else if (rest ~ /^\+/) chomp = "keep"
                state = 1
                next
            }
            val = line
            if (val ~ /^".*"$/ && length(val) >= 2) {
                val = substr(val, 2, length(val)-2)
                gsub(/\\"/, "\"", val)
            } else if (val ~ /^\047.*\047$/ && length(val) >= 2) {
                val = substr(val, 2, length(val)-2)
                gsub(/\047\047/, "\047", val)
            }
            printf "%s", val
            state = 2
            exit
        }
        state==1 {
            if ($0 ~ /^[ \t]*$/) {
                nlines++
                blank[nlines] = 1
                content[nlines] = ""
                next
            }
            ls = leadspaces($0)
            if (!have_base) {
                if (ls <= key_indent) { finalize(); exit }
                base = ls
                have_base = 1
            }
            if (ls < base) { finalize(); exit }
            nlines++
            blank[nlines] = 0
            content[nlines] = substr($0, base+1)
        }
        END { if (state == 1) finalize() }
    ' "$1"
}

# Emit name<TAB>kind-csv for collisions across the three mounted namespaces.
_catalog_collisions() {
    local root d f name
    root="$(_catalog_root)"
    {
        for d in "$root"/skills/*/; do
            [[ -f "$d/SKILL.md" ]] || continue
            name=$(basename "$d")
            [[ "$name" == _* ]] || printf '%s skill\n' "$name"
        done
        for f in "$root"/agents/*.md; do
            [[ -f "$f" ]] || continue
            name=$(basename "$f" .md)
            [[ "$name" == README ]] || printf '%s agent\n' "$name"
        done
        for f in "$root"/commands/*.md; do
            [[ -f "$f" ]] || continue
            name=$(basename "$f" .md)
            [[ "$name" == README ]] || printf '%s command\n' "$name"
        done
    } | LC_ALL=C sort | awk '
        function emit(   out) {
            if (count < 2) return
            if (skill) out = "skill"
            if (agent) out = out (out ? "," : "") "agent"
            if (command) out = out (out ? "," : "") "command"
            printf "%s\t%s\n", current, out
        }
        $1 != current {
            if (current != "") emit()
            current=$1; count=0; skill=0; agent=0; command=0
        }
        $2 == "skill" && !skill { skill=1; count++ }
        $2 == "agent" && !agent { agent=1; count++ }
        $2 == "command" && !command { command=1; count++ }
        END { if (current != "") emit() }
    '
}

_catalog_join_kinds() {
    local csv="$1" first="${1%%,*}" rest="${1#*,}"
    if [[ "$rest" != *,* ]]; then
        printf 'both %s and %s' "$first" "$rest"
    else
        printf '%s, %s, and %s' "$first" "${rest%%,*}" "${rest##*,}"
    fi
}

_catalog_check() {
    local quiet="$1" root skill_dir skill_name skill_file fm_name description desc_len
    local fail=0
    root="$(_catalog_root)"

    for skill_dir in "$root"/skills/*/; do
        [[ -d "$skill_dir" ]] || continue
        skill_name=$(basename "$skill_dir")
        [[ "$skill_name" == _* ]] && continue
        skill_file="$skill_dir/SKILL.md"
        [[ -f "$skill_file" ]] || continue

        fm_name=$(awk '
            NR==1 && /^---[ \t]*$/ { infm=1; next }
            infm && /^---[ \t]*$/ { exit }
            infm && /^name:[ \t]*/ {
                sub(/^name:[ \t]*/, ""); sub(/[ \t]+$/, "")
                sub(/^"/, ""); sub(/"$/, "")
                sub(/^\047/, ""); sub(/\047$/, ""); print; exit
            }
        ' "$skill_file")
        if [[ -z "$fm_name" ]]; then
            echo "ERROR toolkit: $skill_name: SKILL.md frontmatter has no name:" >&2
            fail=1
        elif [[ "$fm_name" != "$skill_name" ]]; then
            echo "ERROR toolkit: $skill_name: frontmatter name: '$fm_name' does not match its directory" >&2
            fail=1
        fi

        if ! awk '
            NR==1 && /^---[ \t]*$/ { infm=1; next }
            infm && /^---[ \t]*$/ { exit }
            infm && /^description:/ { found=1 }
            END { exit !found }
        ' "$skill_file"; then
            echo "ERROR toolkit: $skill_name: SKILL.md frontmatter has no description:" >&2
            fail=1
            continue
        fi

        description="$(_catalog_description "$skill_file"; printf X)"
        description="${description%X}"
        desc_len=${#description}
        if (( desc_len > SCIO_DESC_MAX )); then
            echo "ERROR toolkit: $skill_name: description is $desc_len chars (max $SCIO_DESC_MAX)" >&2
            fail=1
        fi
    done

    # The CRAFT body is always-on text in every consumer's AGENTS.md, so its
    # declared budget binds here. Silent without a craft.yaml, matching
    # craft_render_and_write, so a toolkit with no craft SSOT is unaffected.
    local craft_body craft_lines craft_max
    if craft_body=$(_craft_render_body); then
        craft_max=$(craft_max_lines)
        craft_lines=$(printf '%s\n' "$craft_body" | awk 'END { print NR }')
        if (( craft_lines > craft_max )); then
            echo "ERROR toolkit: craft.yaml: rendered CRAFT body is $craft_lines lines (max $craft_max) — move depth into a skill" >&2
            fail=1
        fi
    fi

    if [[ "$quiet" -eq 0 ]]; then
        local name kinds phrase
        while IFS=$'\t' read -r name kinds; do
            [[ -n "$name" ]] || continue
            phrase=$(_catalog_join_kinds "$kinds")
            echo "scio lint: warning — name '$name' appears as $phrase" >&2
        done < <(_catalog_collisions)
    fi

    return "$fail"
}
