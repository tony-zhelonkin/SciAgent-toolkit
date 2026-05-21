# lib/sciagent/status.sh — sciagent status [--json|--effective|--source <name>]
# Reports active stack, effective tables, shadow info, block-hash drift,
# symlink integrity, and harness presence (per architecture §5).

# shellcheck shell=bash

cmd_status() {
    local mode="text"
    local source_query=""
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --json)       mode="json"; shift ;;
            --effective)  mode="effective"; shift ;;
            --source)
                mode="source"
                source_query="${2:-}"
                if [[ -z "$source_query" ]]; then
                    echo "usage: sciagent status --source <name>" >&2
                    return 1
                fi
                shift 2 ;;
            *)
                echo "sciagent status: unknown flag '$1'" >&2
                return 1 ;;
        esac
    done

    if ! manifest_exists; then
        case "$mode" in
            json)      printf '{"active":false}\n' ;;
            effective) ;;  # silent
            source)    echo "no role active" ;;
            *)         echo "No role active." ;;
        esac
        return 0
    fi

    # Load state into globals consumed by renderers below.
    _status_load_state

    case "$mode" in
        json)      _status_render_json ;;
        effective) _status_render_effective ;;
        source)    _status_render_source "$source_query" ;;
        *)         _status_render_text ;;
    esac
}

# Globals populated by _status_load_state:
#   _STK_BASE, _STK_OVERLAY     stack role names
#   _BLOCK_HASH_OK              "ok" | "drift" | "missing" | "corrupt"
#   _SYMLINKS_OK                "ok" | "broken"
#   SKILL_ORDER[], AGENTS_ORDER[], COMMANDS_ORDER[], INJECTED_LIST[]
#   SKILLS[name]=role, AGENTS_M[name]=role, COMMANDS_M[name]=role
#   SKILL_SHADOWS[name]=role, AGENT_SHADOWS[name]=role, COMMAND_SHADOWS[name]=role
#   OUTPUT_STYLE, OUTPUT_STYLE_ROLE
_status_load_state() {
    local stack
    stack=$(manifest_stack)
    _STK_BASE=$(printf '%s\n' "$stack" | awk '{print $1}')
    _STK_OVERLAY=$(printf '%s\n' "$stack" | awk '{print $2}')

    # Drift check.
    if [[ -f AGENTS.md ]]; then
        block_hash_check AGENTS.md
        case $? in
            0) _BLOCK_HASH_OK="ok" ;;
            1) _BLOCK_HASH_OK="missing" ;;
            2) _BLOCK_HASH_OK="corrupt" ;;
            3) _BLOCK_HASH_OK="drift" ;;
        esac
    else
        _BLOCK_HASH_OK="missing"
    fi

    # Symlink integrity: every manifest SYMLINK should exist and be a symlink.
    _SYMLINKS_OK="ok"
    local kind path
    while read -r kind path _rest; do
        [[ "$kind" == "SYMLINK" ]] || continue
        if [[ ! -L "$path" ]]; then
            _SYMLINKS_OK="broken"
            break
        fi
    done < "$_MANIFEST_PATH"

    # Collect injected skill names.
    INJECTED_LIST=()
    local _ov nm
    while read -r kind _ov nm; do
        [[ "$kind" == "INJECTED" ]] || continue
        INJECTED_LIST+=("$nm")
    done < "$_MANIFEST_PATH"

    # Walk role YAMLs to fill effective tables.
    declare -gA SKILLS=() AGENTS_M=() COMMANDS_M=()
    declare -gA SKILL_SHADOWS=() AGENT_SHADOWS=() COMMAND_SHADOWS=()
    SKILL_ORDER=(); AGENTS_ORDER=(); COMMANDS_ORDER=()
    OUTPUT_STYLE=""; OUTPUT_STYLE_ROLE=""

    local role line k2 n2
    for role in $_STK_BASE ${_STK_OVERLAY:+$_STK_OVERLAY}; do
        [[ "$role" == "_injected" ]] && continue
        while IFS= read -r line; do
            k2="${line%% *}"; n2="${line#* }"
            case "$k2" in
                SKILL)
                    if [[ -n "${SKILLS[$n2]:-}" ]]; then
                        SKILL_SHADOWS[$n2]="${SKILLS[$n2]}"
                    else
                        SKILL_ORDER+=("$n2")
                    fi
                    SKILLS[$n2]="$role" ;;
                AGENT)
                    if [[ -n "${AGENTS_M[$n2]:-}" ]]; then
                        AGENT_SHADOWS[$n2]="${AGENTS_M[$n2]}"
                    else
                        AGENTS_ORDER+=("$n2")
                    fi
                    AGENTS_M[$n2]="$role" ;;
                COMMAND)
                    if [[ -n "${COMMANDS_M[$n2]:-}" ]]; then
                        COMMAND_SHADOWS[$n2]="${COMMANDS_M[$n2]}"
                    else
                        COMMANDS_ORDER+=("$n2")
                    fi
                    COMMANDS_M[$n2]="$role" ;;
                OUTPUT_STYLE)
                    OUTPUT_STYLE="$n2"; OUTPUT_STYLE_ROLE="$role" ;;
            esac
        done < <(role_load "$role")
    done
}

_status_render_text() {
    printf 'Stack:\n'
    local desc
    desc=$(role_description "$_STK_BASE")
    printf '  1. %-10s roles/%s.yaml  %s\n' "$_STK_BASE" "$_STK_BASE" "$desc"
    if [[ -n "$_STK_OVERLAY" ]]; then
        if [[ "$_STK_OVERLAY" == "_injected" ]]; then
            printf '  2. %-10s (synthetic)            [overlay]\n' "_injected"
        else
            desc=$(role_description "$_STK_OVERLAY")
            printf '  2. %-10s roles/%s.yaml  %s  [overlay]\n' "$_STK_OVERLAY" "$_STK_OVERLAY" "$desc"
        fi
    fi
    printf '\n'

    local n note total
    total=$(( ${#SKILL_ORDER[@]} + ${#INJECTED_LIST[@]} ))
    printf 'Skills (%d effective):\n' "$total"
    for n in "${SKILL_ORDER[@]:-}"; do
        [[ -z "$n" ]] && continue
        note=""
        [[ -n "${SKILL_SHADOWS[$n]:-}" ]] && note="    (shadows ${SKILL_SHADOWS[$n]})"
        printf '  %-24s %s%s\n' "$n" "${SKILLS[$n]}" "$note"
    done
    for n in "${INJECTED_LIST[@]:-}"; do
        [[ -z "$n" ]] && continue
        printf '  %-24s %s\n' "$n" "injected"
    done
    printf '\n'

    printf 'Sub-agents (%d effective, Claude-only):\n' "${#AGENTS_ORDER[@]}"
    for n in "${AGENTS_ORDER[@]:-}"; do
        [[ -z "$n" ]] && continue
        note=""
        [[ -n "${AGENT_SHADOWS[$n]:-}" ]] && note="    (shadows ${AGENT_SHADOWS[$n]})"
        printf '  %-24s %s%s\n' "$n" "${AGENTS_M[$n]}" "$note"
    done
    printf '\n'

    printf 'Slash commands (%d effective, Claude-only):\n' "${#COMMANDS_ORDER[@]}"
    for n in "${COMMANDS_ORDER[@]:-}"; do
        [[ -z "$n" ]] && continue
        note=""
        [[ -n "${COMMAND_SHADOWS[$n]:-}" ]] && note="    (shadows ${COMMAND_SHADOWS[$n]})"
        printf '  /%-23s %s%s\n' "$n" "${COMMANDS_M[$n]}" "$note"
    done
    printf '\n'

    # Block + symlinks + harness summary.
    local hash_label
    case "$_BLOCK_HASH_OK" in
        ok)      hash_label="hash OK" ;;
        drift)   hash_label="hash DRIFT" ;;
        corrupt) hash_label="markers corrupted" ;;
        missing) hash_label="block missing" ;;
    esac
    local lines=""
    if [[ -f AGENTS.md ]] && grep -qF '<!-- BEGIN SCIAGENT:ROLES' AGENTS.md; then
        local lb le
        lb=$(grep -nF '<!-- BEGIN SCIAGENT:ROLES' AGENTS.md | head -n1 | cut -d: -f1)
        le=$(grep -nF '<!-- END SCIAGENT:ROLES'   AGENTS.md | head -n1 | cut -d: -f1)
        lines="lines ${lb}-${le}  "
    fi
    printf 'Managed block: AGENTS.md  %s%s\n' "$lines" "$hash_label"

    local sym_label
    case "$_SYMLINKS_OK" in
        ok)     sym_label=".claude/* OK   .agents/* OK" ;;
        broken) sym_label="BROKEN — run \`sciagent activate $_STK_BASE${_STK_OVERLAY:+ $_STK_OVERLAY}\` to repair" ;;
    esac
    printf 'Symlinks:      %s\n' "$sym_label"

    local claude_state pi_state
    if [[ -d .claude ]]; then claude_state="Claude Code detected (.claude/ present)"
    else claude_state="Claude Code not detected (.claude/ absent)"; fi
    if [[ -d .pi ]]; then pi_state="detected (.pi/ present)"
    else pi_state="not detected (.pi/ absent)"; fi
    printf 'Harness:       %s   Pi: %s\n' "$claude_state" "$pi_state"
}

_status_render_effective() {
    local n
    for n in "${SKILL_ORDER[@]:-}";    do [[ -n "$n" ]] && echo "$n"; done
    for n in "${INJECTED_LIST[@]:-}";  do [[ -n "$n" ]] && echo "$n"; done
    for n in "${AGENTS_ORDER[@]:-}";   do [[ -n "$n" ]] && echo "$n"; done
    for n in "${COMMANDS_ORDER[@]:-}"; do [[ -n "$n" ]] && echo "$n"; done
}

_status_render_source() {
    local q="$1" n
    for n in "${SKILL_ORDER[@]:-}"; do
        [[ "$n" == "$q" ]] || continue
        if [[ -n "${SKILL_SHADOWS[$n]:-}" ]]; then
            echo "${SKILLS[$n]} (shadowed by ${SKILLS[$n]}, also in ${SKILL_SHADOWS[$n]})"
        else
            echo "${SKILLS[$n]}"
        fi
        return 0
    done
    for n in "${INJECTED_LIST[@]:-}"; do
        [[ "$n" == "$q" ]] || continue
        echo "injected"; return 0
    done
    for n in "${AGENTS_ORDER[@]:-}"; do
        [[ "$n" == "$q" ]] || continue
        echo "${AGENTS_M[$n]}"; return 0
    done
    for n in "${COMMANDS_ORDER[@]:-}"; do
        [[ "$n" == "$q" ]] || continue
        echo "${COMMANDS_M[$n]}"; return 0
    done
    echo "not in active stack: $q" >&2
    return 1
}

# JSON escape: backslash and double-quote only (control chars not expected
# in role/skill/agent/command names).
_json_esc() {
    local s="$1"
    s="${s//\\/\\\\}"
    s="${s//\"/\\\"}"
    printf '%s' "$s"
}

_json_array_names() {
    local first=1 n
    printf '['
    for n in "$@"; do
        [[ -z "$n" ]] && continue
        if (( first )); then first=0; else printf ','; fi
        printf '"%s"' "$(_json_esc "$n")"
    done
    printf ']'
}

_status_render_json() {
    local first
    printf '{'
    printf '"active":true,'

    # stack
    printf '"stack":['
    printf '{"role":"%s","kind":"base"}' "$(_json_esc "$_STK_BASE")"
    if [[ -n "$_STK_OVERLAY" ]]; then
        printf ',{"role":"%s","kind":"overlay"}' "$(_json_esc "$_STK_OVERLAY")"
    fi
    printf '],'

    # skills (with role + injected list embedded)
    printf '"skills":['
    first=1
    local n
    for n in "${SKILL_ORDER[@]:-}"; do
        [[ -z "$n" ]] && continue
        if (( first )); then first=0; else printf ','; fi
        printf '{"name":"%s","source":"%s"}' "$(_json_esc "$n")" "$(_json_esc "${SKILLS[$n]}")"
    done
    for n in "${INJECTED_LIST[@]:-}"; do
        [[ -z "$n" ]] && continue
        if (( first )); then first=0; else printf ','; fi
        printf '{"name":"%s","source":"injected"}' "$(_json_esc "$n")"
    done
    printf '],'

    # agents
    printf '"agents":['
    first=1
    for n in "${AGENTS_ORDER[@]:-}"; do
        [[ -z "$n" ]] && continue
        if (( first )); then first=0; else printf ','; fi
        printf '{"name":"%s","source":"%s"}' "$(_json_esc "$n")" "$(_json_esc "${AGENTS_M[$n]}")"
    done
    printf '],'

    # commands
    printf '"commands":['
    first=1
    for n in "${COMMANDS_ORDER[@]:-}"; do
        [[ -z "$n" ]] && continue
        if (( first )); then first=0; else printf ','; fi
        printf '{"name":"%s","source":"%s"}' "$(_json_esc "$n")" "$(_json_esc "${COMMANDS_M[$n]}")"
    done
    printf '],'

    # shadowed
    printf '"shadowed":['
    first=1
    for n in "${!SKILL_SHADOWS[@]}"; do
        if (( first )); then first=0; else printf ','; fi
        printf '{"name":"%s","kind":"skill","hidden":"%s","winner":"%s"}' \
            "$(_json_esc "$n")" "$(_json_esc "${SKILL_SHADOWS[$n]}")" "$(_json_esc "${SKILLS[$n]}")"
    done
    for n in "${!AGENT_SHADOWS[@]}"; do
        if (( first )); then first=0; else printf ','; fi
        printf '{"name":"%s","kind":"agent","hidden":"%s","winner":"%s"}' \
            "$(_json_esc "$n")" "$(_json_esc "${AGENT_SHADOWS[$n]}")" "$(_json_esc "${AGENTS_M[$n]}")"
    done
    for n in "${!COMMAND_SHADOWS[@]}"; do
        if (( first )); then first=0; else printf ','; fi
        printf '{"name":"%s","kind":"command","hidden":"%s","winner":"%s"}' \
            "$(_json_esc "$n")" "$(_json_esc "${COMMAND_SHADOWS[$n]}")" "$(_json_esc "${COMMANDS_M[$n]}")"
    done
    printf '],'

    local hash_ok="false"; [[ "$_BLOCK_HASH_OK" == "ok" ]] && hash_ok="true"
    local sym_ok="false";  [[ "$_SYMLINKS_OK"  == "ok" ]] && sym_ok="true"
    local claude="false";  [[ -d .claude ]] && claude="true"
    local pi="false";      [[ -d .pi     ]] && pi="true"
    printf '"block":{"hash_ok":%s,"state":"%s"},' "$hash_ok" "$_BLOCK_HASH_OK"
    printf '"symlinks_ok":%s,' "$sym_ok"
    printf '"harness":{"claude":%s,"pi":%s}' "$claude" "$pi"
    printf '}\n'
}

# cmd_list [roles|skills|agents|commands]
cmd_list() {
    local what="${1:-all}"
    case "$what" in
        roles)    _list_roles ;;
        skills)   _list_skills ;;
        agents)   _list_agents ;;
        commands) _list_commands ;;
        all)
            printf 'Roles:\n';    _list_roles
            printf '\nSkills:\n';   _list_skills
            printf '\nAgents:\n';   _list_agents
            printf '\nCommands:\n'; _list_commands
            ;;
        *)
            echo "sciagent list: unknown category '$what' (use: roles|skills|agents|commands)" >&2
            return 1 ;;
    esac
}

_list_roles() {
    local f name desc
    for f in "$SCIAGENT_TOOLKIT"/roles/*.yaml; do
        [[ -f "$f" ]] || continue
        name=$(basename "$f" .yaml)
        desc=$(role_scalar "$f" description)
        printf '  %-26s %s\n' "$name" "$desc"
    done
}

_list_skills() {
    local d name
    for d in "$SCIAGENT_TOOLKIT"/skills/*/; do
        [[ -f "$d/SKILL.md" ]] || continue
        name=$(basename "$d")
        [[ "$name" == "_TEMPLATE" ]] && continue
        printf '  %s\n' "$name"
    done
}

_list_agents() {
    local f name
    for f in "$SCIAGENT_TOOLKIT"/agents/*.md; do
        [[ -f "$f" ]] || continue
        name=$(basename "$f" .md)
        [[ "$name" == "README" ]] && continue
        printf '  %s\n' "$name"
    done
}

_list_commands() {
    local f name
    for f in "$SCIAGENT_TOOLKIT"/commands/*.md; do
        [[ -f "$f" ]] || continue
        name=$(basename "$f" .md)
        [[ "$name" == "README" ]] && continue
        printf '  /%s\n' "$name"
    done
}
