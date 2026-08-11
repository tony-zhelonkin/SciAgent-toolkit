# lib/sciagent/status.sh — sciagent status [--json|--effective|--source <name>]
# Reports active stack, effective tables, shadow info, block-hash drift,
# symlink integrity, and harness presence.
# Stack-walking is delegated to stack.sh:stack_walk.
#
# Also hosts cmd_list (sciagent list ...).

# shellcheck shell=bash

# Source the frontmatter parser if available (graceful degradation if absent).
_status_load_frontmatter() {
    local fm_path="${SCIAGENT_TOOLKIT}/lib/sciagent/frontmatter.sh"
    if [[ -f "$fm_path" ]] && ! declare -F _fm_extract >/dev/null 2>&1; then
        # shellcheck source=/dev/null
        . "$fm_path"
    fi
}

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
#   SKILL_ORDER[], AGENTS_ORDER[], COMMANDS_ORDER[]
#   SKILLS[name]=role, AGENTS_M[name]=role, COMMANDS_M[name]=role
#   SKILL_SHADOWS[name]=role, AGENT_SHADOWS[name]=role, COMMAND_SHADOWS[name]=role
# output_style is no longer role-scoped (Phase 5d) — status does not track it;
# see `.claude/settings.local.json`'s outputStyle key for the applied value.
_status_load_state() {
    local stack
    stack=$(manifest_stack)
    _STK_BASE=$(printf '%s\n' "$stack" | awk '{print $1}')
    _STK_OVERLAY=$(printf '%s\n' "$stack" | awk '{print $2}')

    # Drift check. ROLES explicitly — this is the effective-stack block
    # activate.sh writes (see block.sh header: three ids exist — ROLES,
    # CRAFT, CONTEXT — status only reports ROLES drift here; CRAFT's own
    # drift is surfaced by craft_verb.sh/lint.sh, not duplicated here).
    if [[ -f AGENTS.md ]]; then
        block_hash_check AGENTS.md ROLES
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

    # Collect skill names mounted on disk per manifest (one per `.claude/skills/*`
    # symlink). The walker above only sees role-yaml entries, so anything mounted
    # by some earlier verb shows up here and not in SKILL_ORDER.
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
        esac
    done < <(stack_walk "$_STK_BASE" "$_STK_OVERLAY")

    # Skills present on disk per the manifest but absent from every role yaml.
    # Normally empty now that the whole catalog is mounted; a legacy manifest
    # from a curated mount can still carry entries. Depends on SKILL_ORDER and
    # MANIFEST_SKILLS being populated above. Use a temp associative array as a
    # set.
    declare -A _inh_declared_set=()
    local _inh_n
    for _inh_n in "${SKILL_ORDER[@]:-}";     do [[ -n "$_inh_n" ]] && _inh_declared_set[$_inh_n]=1; done
    INHERITED_SKILLS=()
    for _inh_n in "${MANIFEST_SKILLS[@]:-}"; do
        [[ -z "$_inh_n" ]] && continue
        [[ -n "${_inh_declared_set[$_inh_n]:-}" ]] && continue
        INHERITED_SKILLS+=("$_inh_n")
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
    total=$(( ${#SKILL_ORDER[@]} + ${#INHERITED_SKILLS[@]} ))
    printf 'Skills (%d effective):\n' "$total"
    for n in "${SKILL_ORDER[@]:-}"; do
        [[ -z "$n" ]] && continue
        note=""
        [[ -n "${SKILL_SHADOWS[$n]:-}" ]] && note="    (shadows ${SKILL_SHADOWS[$n]})"
        printf '  %-24s %s%s\n' "$n" "${SKILLS[$n]}" "$note"
    done
    for n in "${INHERITED_SKILLS[@]:-}"; do
        [[ -z "$n" ]] && continue
        printf '  %-24s %s\n' "$n" "on disk, not declared"
    done
    printf '\n'

    local agent_total=${#AGENTS_ORDER[@]}
    printf 'Sub-agents (%d effective, Claude-only):\n' "$agent_total"

    # Try to load frontmatter parser for enriched agent info.
    _status_load_frontmatter

    for n in "${AGENTS_ORDER[@]:-}"; do
        [[ -z "$n" ]] && continue
        note=""
        [[ -n "${AGENT_SHADOWS[$n]:-}" ]] && note="    (shadows ${AGENT_SHADOWS[$n]})"
        # Enrich with domain[0] and description brief if frontmatter.sh is loaded
        # and the activated agent file exists.
        local _agent_extra=""
        if declare -F _fm_extract >/dev/null 2>&1; then
            local _agent_md=".claude/agents/${n}.md"
            if [[ -f "$_agent_md" ]]; then
                local _agent_fm _agent_domain _agent_desc
                _agent_fm=$(_fm_extract "$_agent_md")
                _agent_domain=$(_fm_list domain <<< "$_agent_fm" | head -1)
                _agent_desc=$(_fm_description "$_agent_md")
                _agent_desc="${_agent_desc:0:60}"
                if [[ -n "$_agent_domain" || -n "$_agent_desc" ]]; then
                    _agent_extra=$(printf '  %-20s  "%s"' "${_agent_domain:-}" "$_agent_desc")
                fi
            fi
        fi
        printf '  %-24s %-8s%s%s\n' "$n" "${AGENTS_M[$n]}" "$_agent_extra" "$note"
    done
    printf '\n'

    local command_total=${#COMMANDS_ORDER[@]}
    printf 'Slash commands (%d effective, Claude-only):\n' "$command_total"
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
    local lines="" range lb le
    if range=$(block_line_range AGENTS.md ROLES); then
        read -r lb le <<< "$range"
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

    # Notes section — emit only when the active stack has at least one
    # cross-namespace collision actually realised on disk (>=2 kinds of the
    # same name mounted by this stack). Annotated against the CI allowlist
    # so the user can tell intentional family overlaps from accidents.
    _status_render_notes
}

# _status_render_notes
# Emits a Notes: section when the active stack has actionable signal:
#   (1) drift — a manifest-pinned role whose .yaml no longer exists in the catalog
#   (2) collisions — a name mounted as >=2 kinds (skill, agent, command)
# Prints nothing when the stack is clean. Both signals can co-occur.
_status_render_notes() {
    local -a _active_lines=()

    # --- (1) Drift detection (no external helper required) ---
    # Check manifest-pinned base role.
    if [[ -n "$_STK_BASE" ]] && ! role_exists "$_STK_BASE"; then
        _active_lines+=("  - stack role '$_STK_BASE' is no longer in the toolkit catalog (roles/$_STK_BASE.yaml missing); run 'sciagent deactivate' to reset")
    fi
    # Check overlay, skipping synthetic _injected value.
    if [[ -n "$_STK_OVERLAY" && "$_STK_OVERLAY" != "_injected" ]] && ! role_exists "$_STK_OVERLAY"; then
        _active_lines+=("  - overlay role '$_STK_OVERLAY' is no longer in the toolkit catalog (roles/$_STK_OVERLAY.yaml missing); run 'sciagent deactivate' to reset")
    fi

    # --- (2) Collision detection (requires optional collisions_enumerate helper) ---
    if declare -F collisions_enumerate >/dev/null 2>&1; then
        # Build the active mount set: kind -> set of names.
        declare -A _mounted_skill=() _mounted_agent=() _mounted_command=()
        local n
        for n in "${SKILL_ORDER[@]:-}";     do [[ -n "$n" ]] && _mounted_skill[$n]=1;   done
        for n in "${MANIFEST_SKILLS[@]:-}"; do [[ -n "$n" ]] && _mounted_skill[$n]=1;   done
        for n in "${AGENTS_ORDER[@]:-}";    do [[ -n "$n" ]] && _mounted_agent[$n]=1;   done
        for n in "${COMMANDS_ORDER[@]:-}";  do [[ -n "$n" ]] && _mounted_command[$n]=1; done

        # Walk every collision; intersect with what this stack actually mounted.
        local col_name col_kinds
        while IFS=$'\t' read -r col_name col_kinds; do
            [[ -z "$col_name" ]] && continue
            local -a _active_kinds=()
            case ",$col_kinds," in *,skill,*)
                [[ -n "${_mounted_skill[$col_name]:-}" ]] && _active_kinds+=("skill") ;;
            esac
            case ",$col_kinds," in *,agent,*)
                [[ -n "${_mounted_agent[$col_name]:-}" ]] && _active_kinds+=("agent") ;;
            esac
            case ",$col_kinds," in *,command,*)
                [[ -n "${_mounted_command[$col_name]:-}" ]] && _active_kinds+=("command") ;;
            esac
            (( ${#_active_kinds[@]} >= 2 )) || continue
            # Render the kinds list as "a and b" / "a, b, and c".
            local kinds_phrase
            case ${#_active_kinds[@]} in
                2) kinds_phrase="${_active_kinds[0]} and ${_active_kinds[1]}" ;;
                3) kinds_phrase="${_active_kinds[0]}, ${_active_kinds[1]}, and ${_active_kinds[2]}" ;;
                *) kinds_phrase="${_active_kinds[*]}" ;;
            esac
            # Annotate against the allowlist. Status leans on the same file the
            # CI test reads; format mirrors the line shape there (`<name> <csv>`).
            local annotation
            if _status_collision_in_allowlist "$col_name" "$col_kinds"; then
                # Qualify the path: status runs in the user's project dir, so a
                # bare "tests/collision-allowlist.txt" looks like a sibling file
                # that doesn't exist. The file lives in the toolkit submodule.
                annotation="intentional family overlap per the toolkit's tests/collision-allowlist.txt"
            else
                annotation="UNEXPECTED — run 'sciagent validate' for details"
            fi
            _active_lines+=("  - '$col_name' appears as both $kinds_phrase ($annotation)")
        done < <(collisions_enumerate)
    fi


    (( ${#_active_lines[@]} == 0 )) && return 0

    printf '\nNotes:\n'
    local line
    for line in "${_active_lines[@]}"; do
        printf '%s\n' "$line"
    done
}

# _status_collision_in_allowlist <name> <kinds-csv>
# Returns 0 when the toolkit's allowlist records this exact (name, kinds)
# pair. Returns 1 when the name is absent OR the kinds csv differs — the
# CI test treats a kinds-csv drift as a hard-fail, so status mirrors it
# as "UNEXPECTED" rather than "intentional".
_status_collision_in_allowlist() {
    local want_name="$1" want_kinds="$2"
    local allowlist="${SCIAGENT_TOOLKIT:-}/tests/collision-allowlist.txt"
    [[ -f "$allowlist" ]] || return 1
    local line name kinds
    while IFS= read -r line || [[ -n "$line" ]]; do
        line="${line%%#*}"
        line="${line#"${line%%[![:space:]]*}"}"
        line="${line%"${line##*[![:space:]]}"}"
        [[ -z "$line" ]] && continue
        name="${line%%[[:space:]]*}"
        kinds="${line#"$name"}"
        kinds="${kinds#"${kinds%%[![:space:]]*}"}"
        [[ "$name" == "$want_name" && "$kinds" == "$want_kinds" ]] && return 0
    done < "$allowlist"
    return 1
}

_status_render_effective() {
    local n
    for n in "${SKILL_ORDER[@]:-}";      do [[ -n "$n" ]] && echo "$n"; done
    for n in "${INHERITED_SKILLS[@]:-}"; do [[ -n "$n" ]] && echo "$n"; done
    for n in "${AGENTS_ORDER[@]:-}";     do [[ -n "$n" ]] && echo "$n"; done
    for n in "${COMMANDS_ORDER[@]:-}";   do [[ -n "$n" ]] && echo "$n"; done
    return 0
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
    # On disk per the manifest, declared by nothing.
    for n in "${INHERITED_SKILLS[@]:-}"; do
        [[ "$n" == "$q" ]] || continue
        echo "on disk, not declared"
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

    # skills (with role embedded)
    printf '"skills":['
    first=1
    local n
    for n in "${SKILL_ORDER[@]:-}"; do
        [[ -z "$n" ]] && continue
        if (( first )); then first=0; else printf ','; fi
        printf '{"name":"%s","source":"%s"}' "$(_json_esc "$n")" "$(_json_esc "${SKILLS[$n]}")"
    done
    for n in "${INHERITED_SKILLS[@]:-}"; do
        [[ -z "$n" ]] && continue
        if (( first )); then first=0; else printf ','; fi
        printf '{"name":"%s","source":"inherited"}' "$(_json_esc "$n")"
    done
    printf '],'

    # agents — enriched with domain and description_brief when frontmatter.sh is loaded.
    _status_load_frontmatter
    printf '"agents":['
    first=1
    for n in "${AGENTS_ORDER[@]:-}"; do
        [[ -z "$n" ]] && continue
        if (( first )); then first=0; else printf ','; fi
        local _aj_domain="" _aj_desc=""
        if declare -F _fm_extract >/dev/null 2>&1; then
            local _aj_md=".claude/agents/${n}.md"
            if [[ -f "$_aj_md" ]]; then
                local _aj_fm
                _aj_fm=$(_fm_extract "$_aj_md")
                _aj_domain=$(_fm_list domain <<< "$_aj_fm" | paste -sd, -)
                _aj_desc=$(_fm_description "$_aj_md")
            fi
        fi
        printf '{"name":"%s","source":"%s","domain":"%s","description_brief":"%s"}' \
            "$(_json_esc "$n")" "$(_json_esc "${AGENTS_M[$n]}")" \
            "$(_json_esc "$_aj_domain")" "$(_json_esc "$_aj_desc")"
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

    # stale_roles: manifest-pinned roles whose .yaml no longer exists in catalog.
    printf '"stale_roles":['
    first=1
    local role
    for role in "$_STK_BASE" "$_STK_OVERLAY"; do
        [[ -z "$role" || "$role" == "_injected" ]] && continue
        role_exists "$role" && continue
        if (( first )); then first=0; else printf ','; fi
        printf '"%s"' "$(_json_esc "$role")"
    done
    printf '],'

    printf '"harness":{"claude":%s,"pi":%s}' "$claude" "$pi"
    printf '}\n'
}

# cmd_list [roles|skills|agents|commands|role <name>]
#
# Subcommands:
#   (no arg)             — list everything (roles, skills, agents, commands)
#   roles                — list all roles with one-line description + counts
#   role <name>          — detailed view of a specific role (skills/agents/commands)
#   skills               — list all available skills
#   agents               — list all available agents
#   commands             — list all available slash commands
cmd_list() {
    local what="${1:-all}"
    case "$what" in
        roles)    _list_roles ;;
        role)
            shift || true
            if [[ -z "${1:-}" ]]; then
                echo "sciagent list role: missing <name> argument" >&2
                return 1
            fi
            _list_role_detail "$1"
            ;;
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
            echo "sciagent list: unknown category '$what' (use: roles|role <name>|skills|agents|commands)" >&2
            return 1 ;;
    esac
}

_list_roles() {
    local f name desc skill_count agent_count cmd_count
    for f in "$SCIAGENT_TOOLKIT"/roles/*.yaml; do
        [[ -f "$f" ]] || continue
        name=$(basename "$f" .yaml)
        desc=$(role_scalar "$f" description)
        skill_count=$(role_array "$f" skills | wc -l | tr -d ' ')
        agent_count=$(role_array "$f" agents | wc -l | tr -d ' ')
        cmd_count=$(role_array "$f" commands | wc -l | tr -d ' ')
        # Suppress zero counts when array is empty (wc -l of empty = 0, still ok)
        printf '  %-26s  skills:%-3s agents:%-3s cmds:%-3s  %s\n' \
            "$name" "$skill_count" "$agent_count" "$cmd_count" "$desc"
    done
}

# _list_role_detail <name>
# Print a detailed RPG-hero view of a role: description, skills, agents, commands.
_list_role_detail() {
    local name="$1"
    local rfile
    rfile="$SCIAGENT_TOOLKIT/roles/${name}.yaml"
    if [[ ! -f "$rfile" ]]; then
        echo "sciagent list role: role '$name' not found (looked for roles/${name}.yaml)" >&2
        return 1
    fi
    local desc
    desc=$(role_scalar "$rfile" description)
    printf 'Role: %s — %s\n\n' "$name" "$desc"

    local -a skills=() agents=() cmds=()
    while IFS= read -r item; do [[ -n "$item" ]] && skills+=("$item"); done \
        < <(role_array "$rfile" skills)
    while IFS= read -r item; do [[ -n "$item" ]] && agents+=("$item"); done \
        < <(role_array "$rfile" agents)
    while IFS= read -r item; do [[ -n "$item" ]] && cmds+=("$item"); done \
        < <(role_array "$rfile" commands)

    printf 'Skills  (%d): %s\n' "${#skills[@]}" "$(IFS=', '; printf '%s' "${skills[*]:-}")"
    printf 'Agents  (%d): %s\n' "${#agents[@]}" "$(IFS=', '; printf '%s' "${agents[*]:-}")"
    printf 'Commands(%d): %s\n' "${#cmds[@]}" "$(IFS=', '; printf '%s' "${cmds[*]:-}")"
    return 0
}

_list_skills() {
    local d name
    # Active skills: flat dirs at skills/<name>/. Underscore-prefixed dirs
    # (_TEMPLATE, _attic, _archive) are scaffolding, not active skills.
    for d in "$SCIAGENT_TOOLKIT"/skills/*/; do
        [[ -f "$d/SKILL.md" ]] || continue
        name=$(basename "$d")
        [[ "$name" == _* ]] && continue
        printf '  %s\n' "$name"
    done

    # Attic: retired, reference-only skills (skills/_attic/<name>/). Listed
    # separately so they are never mistaken for active/available skills.
    local attic_dir="$SCIAGENT_TOOLKIT/skills/_attic"
    if [[ -d "$attic_dir" ]]; then
        local a aname
        local -a attic=()
        for a in "$attic_dir"/*/; do
            [[ -f "$a/SKILL.md" ]] || continue
            attic+=("$(basename "$a")")
        done
        if (( ${#attic[@]} > 0 )); then
            printf '\n  Attic (retired, reference-only — not installed by any role):\n'
            for aname in "${attic[@]}"; do
                printf '    _attic/%s\n' "$aname"
            done
        fi
    fi
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
