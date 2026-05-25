# lib/sciagent/status.sh — sciagent status [--json|--effective|--source <name>]
# Reports active stack, effective tables, shadow info, block-hash drift,
# symlink integrity, and harness presence (per architecture §5).
# Stack-walking is delegated to stack.sh:stack_walk.

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
#   INJECTED_OVERLAY[]          parallel array: which overlay each injected skill is in
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

    # Symlink integrity: every manifest symlink should exist and be a symlink.
    _SYMLINKS_OK="ok"
    local path
    while IFS= read -r path; do
        [[ -z "$path" ]] && continue
        if [[ ! -L "$path" ]]; then
            _SYMLINKS_OK="broken"
            break
        fi
    done < <(manifest_symlinks)

    # Collect injected skill names and their overlay from the manifest.
    INJECTED_LIST=()
    INJECTED_OVERLAY=()
    local _ov nm _via _kind
    while IFS='|' read -r _ov nm _via _kind; do
        [[ -n "$_ov" ]] || continue
        INJECTED_LIST+=("$nm")
        INJECTED_OVERLAY+=("$_ov")
    done < <(manifest_injected)

    # Collect skill names mounted on disk per manifest (one per `.claude/skills/*`
    # symlink). Skills can land here via the role yaml, via `inject`, OR via
    # the `metadata.requires:` transitive closure (activate.sh mounts every
    # dep). The walker above only sees role-yaml entries — anything mounted
    # purely via requires-inheritance shows up here and not in SKILL_ORDER.
    MANIFEST_SKILLS=()
    local _ms_path _ms_name
    while IFS= read -r _ms_path; do
        [[ -z "$_ms_path" ]] && continue
        [[ "$_ms_path" == .claude/skills/* ]] || continue
        _ms_name="${_ms_path#.claude/skills/}"
        MANIFEST_SKILLS+=("$_ms_name")
    done < <(manifest_symlinks)

    # Walk role YAMLs via stack_walk to fill effective tables.
    declare -gA SKILLS=() AGENTS_M=() COMMANDS_M=()
    declare -gA SKILL_SHADOWS=() AGENT_SHADOWS=() COMMAND_SHADOWS=()
    SKILL_ORDER=(); AGENTS_ORDER=(); COMMANDS_ORDER=()
    OUTPUT_STYLE=""; OUTPUT_STYLE_ROLE=""

    local kind n2 provider shadow
    while IFS=$'\t' read -r kind n2 provider shadow; do
        case "$kind" in
            SKILL)
                SKILLS[$n2]="$provider"
                SKILL_SHADOWS[$n2]="$shadow"
                SKILL_ORDER+=("$n2")
                ;;
            AGENT)
                AGENTS_M[$n2]="$provider"
                AGENT_SHADOWS[$n2]="$shadow"
                AGENTS_ORDER+=("$n2")
                ;;
            COMMAND)
                COMMANDS_M[$n2]="$provider"
                COMMAND_SHADOWS[$n2]="$shadow"
                COMMANDS_ORDER+=("$n2")
                ;;
            STYLE)
                OUTPUT_STYLE="$n2"
                OUTPUT_STYLE_ROLE="$provider"
                ;;
        esac
    done < <(stack_walk "$_STK_BASE" "$_STK_OVERLAY")
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

    # Compute skills mounted purely via `metadata.requires:` inheritance —
    # they sit on disk per the manifest but are not in any role yaml or the
    # injected list. Use a temp associative array as a set.
    declare -A _declared_set=()
    local n
    for n in "${SKILL_ORDER[@]:-}";   do [[ -n "$n" ]] && _declared_set[$n]=1; done
    for n in "${INJECTED_LIST[@]:-}"; do [[ -n "$n" ]] && _declared_set[$n]=1; done
    local -a INHERITED_SKILLS=()
    for n in "${MANIFEST_SKILLS[@]:-}"; do
        [[ -z "$n" ]] && continue
        [[ -n "${_declared_set[$n]:-}" ]] && continue
        INHERITED_SKILLS+=("$n")
    done

    local note total
    total=$(( ${#SKILL_ORDER[@]} + ${#INJECTED_LIST[@]} + ${#INHERITED_SKILLS[@]} ))
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
    for n in "${INHERITED_SKILLS[@]:-}"; do
        [[ -z "$n" ]] && continue
        printf '  %-24s %s\n' "$n" "inherited via requires:"
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
            echo "${SKILLS[$n]} (shadows ${SKILL_SHADOWS[$n]})"
        else
            echo "${SKILLS[$n]}"
        fi
        return 0
    done
    # Injected skills: report overlay name for clarity.
    local i _ninj_list=${#INJECTED_LIST[@]}
    for (( i=0; i<_ninj_list; i++ )); do
        [[ "${INJECTED_LIST[$i]}" == "$q" ]] || continue
        local ov="${INJECTED_OVERLAY[$i]:-_injected}"
        echo "injected (into $ov)"
        return 0
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
    local i _ninj_json=${#INJECTED_LIST[@]}
    for (( i=0; i<_ninj_json; i++ )); do
        n="${INJECTED_LIST[$i]}"
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
        [[ -z "${SKILL_SHADOWS[$n]:-}" ]] && continue
        if (( first )); then first=0; else printf ','; fi
        printf '{"name":"%s","kind":"skill","hidden":"%s","winner":"%s"}' \
            "$(_json_esc "$n")" "$(_json_esc "${SKILL_SHADOWS[$n]}")" "$(_json_esc "${SKILLS[$n]}")"
    done
    for n in "${!AGENT_SHADOWS[@]}"; do
        [[ -z "${AGENT_SHADOWS[$n]:-}" ]] && continue
        if (( first )); then first=0; else printf ','; fi
        printf '{"name":"%s","kind":"agent","hidden":"%s","winner":"%s"}' \
            "$(_json_esc "$n")" "$(_json_esc "${AGENT_SHADOWS[$n]}")" "$(_json_esc "${AGENTS_M[$n]}")"
    done
    for n in "${!COMMAND_SHADOWS[@]}"; do
        [[ -z "${COMMAND_SHADOWS[$n]:-}" ]] && continue
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

# cmd_list [roles|skills|agents|commands|deps <skill>|dependents <skill>]
cmd_list() {
    local what="${1:-all}"
    case "$what" in
        roles)    _list_roles ;;
        skills)   _list_skills ;;
        agents)   _list_agents ;;
        commands) _list_commands ;;
        deps)
            shift || true
            if [[ -z "${1:-}" ]]; then
                echo "sciagent list deps: missing <skill> argument" >&2
                return 1
            fi
            _list_deps "$1"
            ;;
        dependents)
            shift || true
            if [[ -z "${1:-}" ]]; then
                echo "sciagent list dependents: missing <skill> argument" >&2
                return 1
            fi
            _list_dependents "$1"
            ;;
        all)
            printf 'Roles:\n';    _list_roles
            printf '\nSkills:\n';   _list_skills
            printf '\nAgents:\n';   _list_agents
            printf '\nCommands:\n'; _list_commands
            ;;
        *)
            echo "sciagent list: unknown category '$what' (use: roles|skills|agents|commands|deps <skill>|dependents <skill>)" >&2
            return 1 ;;
    esac
}

# _list_deps <skill> — print transitive requires-closure, topo-sorted (leaves first).
# Excludes the input skill itself from output. Returns 1 with a stderr error
# if the target skill (or any of its requires) cannot be resolved.
_list_deps() {
    local target="$1"
    # Validate up-front: the resolver writes its own stderr; we just need to
    # honour its exit code so an unknown skill propagates as exit 1.
    local resolved
    if ! resolved=$(skill_resolve_transitive "$target" 2>&1 >/dev/null) \
        && [[ -n "$resolved" ]]; then
        # Print the captured stderr stanza and propagate failure.
        printf '%s\n' "$resolved" >&2
        return 1
    fi
    # Re-run for stdout collection (cheap; the toolkit's skill set is small).
    local n
    while IFS= read -r n; do
        [[ "$n" == "$target" ]] && continue
        printf '%s\n' "$n"
    done < <(skill_resolve_transitive "$target")
}

# _list_dependents <skill> — print direct (one-level) dependents: skills whose
# `metadata.requires:` includes the input skill. One name per line, sorted.
# Returns 1 with a stderr error if the target skill itself doesn't exist;
# an empty dependents set on a known skill still exits 0.
_list_dependents() {
    local target="$1"
    if ! skill_frontmatter_path "$target" >/dev/null 2>&1; then
        echo "sciagent list dependents: skill '$target' not found" >&2
        return 1
    fi
    local d name
    local -a results=()
    for d in "$SCIAGENT_TOOLKIT"/skills/*/; do
        [[ -f "$d/SKILL.md" ]] || continue
        name=$(basename "$d")
        [[ "$name" == "_TEMPLATE" ]] && continue
        [[ "$name" == "$target" ]] && continue
        if skill_read_requires "$name" 2>/dev/null | grep -Fxq "$target"; then
            results+=("$name")
        fi
    done
    printf '%s\n' "${results[@]+"${results[@]}"}" | sort -u | grep -v '^$' || true
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
    # Recursive walk so subfolders under agents/ are transparent to `list`.
    local f name
    while IFS= read -r f; do
        name=$(basename "$f" .md)
        [[ "$name" == "README" ]] && continue
        printf '  %s\n' "$name"
    done < <(find "$SCIAGENT_TOOLKIT/agents" -type f -name '*.md' -not -path '*/.*' 2>/dev/null | sort)
}

_list_commands() {
    # Recursive walk so subfolders under commands/ are transparent to `list`.
    local f name
    while IFS= read -r f; do
        name=$(basename "$f" .md)
        [[ "$name" == "README" ]] && continue
        printf '  /%s\n' "$name"
    done < <(find "$SCIAGENT_TOOLKIT/commands" -type f -name '*.md' -not -path '*/.*' 2>/dev/null | sort)
}
