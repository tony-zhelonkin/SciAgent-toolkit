# lib/sciagent/stack.sh — shared stack-walking and block-body rendering.
#
# Consumers: activate.sh, inject.sh, status.sh.
# The dispatcher (bin/sciagent) does not source this directly; verb modules do.

# shellcheck shell=bash

# stack_walk <role1> [role2]
# Emits the resolved effective table in tagged form, one entry per line:
#   SKILL    <name>  <providing-role>  <shadowed-roles-csv-or-empty>
#   AGENT    <name>  <providing-role>  <shadowed>
#   COMMAND  <name>  <providing-role>  <shadowed>
#   STYLE    <name>  <providing-role>  <shadowed>
#
# Last-wins resolution is applied across the stack. Output order:
#   SKILL entries (insertion order), then AGENT, then COMMAND, then STYLE.
# The _injected synthetic role is silently skipped (no YAML file).
stack_walk() {
    local base="$1"
    local overlay="${2:-}"

    declare -A _sw_skills=() _sw_agents=() _sw_commands=()
    declare -A _sw_skill_shadows=() _sw_agent_shadows=() _sw_command_shadows=()
    local -a _sw_skill_order=() _sw_agent_order=() _sw_command_order=()
    local _sw_style="" _sw_style_role="" _sw_style_shadow=""

    local role line kind name
    for role in "$base" ${overlay:+"$overlay"}; do
        [[ "$role" == "_injected" ]] && continue
        while IFS= read -r line; do
            kind="${line%% *}"
            name="${line#* }"
            case "$kind" in
                SKILL)
                    if [[ -n "${_sw_skills[$name]:-}" ]]; then
                        # Already present: record the previous provider as shadowed.
                        local prev="${_sw_skills[$name]}"
                        if [[ -n "${_sw_skill_shadows[$name]:-}" ]]; then
                            _sw_skill_shadows[$name]="${_sw_skill_shadows[$name]},$prev"
                        else
                            _sw_skill_shadows[$name]="$prev"
                        fi
                    else
                        _sw_skill_order+=("$name")
                    fi
                    _sw_skills[$name]="$role"
                    ;;
                AGENT)
                    if [[ -n "${_sw_agents[$name]:-}" ]]; then
                        local prev="${_sw_agents[$name]}"
                        if [[ -n "${_sw_agent_shadows[$name]:-}" ]]; then
                            _sw_agent_shadows[$name]="${_sw_agent_shadows[$name]},$prev"
                        else
                            _sw_agent_shadows[$name]="$prev"
                        fi
                    else
                        _sw_agent_order+=("$name")
                    fi
                    _sw_agents[$name]="$role"
                    ;;
                COMMAND)
                    if [[ -n "${_sw_commands[$name]:-}" ]]; then
                        local prev="${_sw_commands[$name]}"
                        if [[ -n "${_sw_command_shadows[$name]:-}" ]]; then
                            _sw_command_shadows[$name]="${_sw_command_shadows[$name]},$prev"
                        else
                            _sw_command_shadows[$name]="$prev"
                        fi
                    else
                        _sw_command_order+=("$name")
                    fi
                    _sw_commands[$name]="$role"
                    ;;
                OUTPUT_STYLE)
                    if [[ -n "$_sw_style" ]]; then
                        _sw_style_shadow="$_sw_style_role"
                    fi
                    _sw_style="$name"
                    _sw_style_role="$role"
                    ;;
            esac
        done < <(role_load "$role")
    done

    # Emit: SKILL entries first, then AGENT, COMMAND, STYLE.
    local n
    for n in "${_sw_skill_order[@]+"${_sw_skill_order[@]}"}"; do
        printf 'SKILL\t%s\t%s\t%s\n' "$n" "${_sw_skills[$n]}" "${_sw_skill_shadows[$n]:-}"
    done
    for n in "${_sw_agent_order[@]+"${_sw_agent_order[@]}"}"; do
        printf 'AGENT\t%s\t%s\t%s\n' "$n" "${_sw_agents[$n]}" "${_sw_agent_shadows[$n]:-}"
    done
    for n in "${_sw_command_order[@]+"${_sw_command_order[@]}"}"; do
        printf 'COMMAND\t%s\t%s\t%s\n' "$n" "${_sw_commands[$n]}" "${_sw_command_shadows[$n]:-}"
    done
    if [[ -n "$_sw_style" ]]; then
        printf 'STYLE\t%s\t%s\t%s\n' "$_sw_style" "$_sw_style_role" "${_sw_style_shadow:-}"
    fi
}

# render_block_body <base> <overlay> [injected-skill ...]
# Produces the body text that goes between the BEGIN/END managed-block markers
# per architecture §4. Reads the role YAMLs via stack_walk.
# Injected skills are listed in a separate "## Injected (overlay)" subsection.
render_block_body() {
    local base="$1"; shift
    local overlay="$1"; shift
    local -a injected=( "$@" )

    # Run stack walk and load results into arrays.
    declare -A _rb_skills=() _rb_agents=() _rb_commands=()
    declare -A _rb_skill_shadows=() _rb_agent_shadows=() _rb_command_shadows=()
    local -a _rb_skill_order=() _rb_agent_order=() _rb_command_order=()
    local _rb_style="" _rb_style_role=""

    local kind name provider shadow
    while IFS=$'\t' read -r kind name provider shadow; do
        case "$kind" in
            SKILL)
                _rb_skills[$name]="$provider"
                _rb_skill_shadows[$name]="$shadow"
                _rb_skill_order+=("$name")
                ;;
            AGENT)
                _rb_agents[$name]="$provider"
                _rb_agent_shadows[$name]="$shadow"
                _rb_agent_order+=("$name")
                ;;
            COMMAND)
                _rb_commands[$name]="$provider"
                _rb_command_shadows[$name]="$shadow"
                _rb_command_order+=("$name")
                ;;
            STYLE)
                _rb_style="$name"
                _rb_style_role="$provider"
                ;;
        esac
    done < <(stack_walk "$base" "$overlay")

    {
        printf '# Active roles\n\n'
        printf 'Stack (in order, last-wins on name collisions):\n'
        local desc
        desc=$(role_description "$base")
        printf '1. **%s** — `roles/%s.yaml` — %s\n' "$base" "$base" "$desc"
        if [[ -n "$overlay" ]]; then
            if [[ "$overlay" == "_injected" ]]; then
                printf '2. **_injected** — (synthetic overlay for injected skills)  *(overlay)*\n'
            else
                desc=$(role_description "$overlay")
                printf '2. **%s** — `roles/%s.yaml` — %s  *(overlay)*\n' "$overlay" "$overlay" "$desc"
            fi
        fi
        printf '\n'

        local n s
        if [[ ${#_rb_skill_order[@]} -gt 0 ]]; then
            printf '## Skills (effective)\n'
            for n in "${_rb_skill_order[@]}"; do
                s=""
                [[ -n "${_rb_skill_shadows[$n]:-}" ]] && s=" [shadows ${_rb_skill_shadows[$n]}]"
                printf -- '- `%s` — (%s)%s\n' "$n" "${_rb_skills[$n]}" "$s"
            done
            printf '\n'
        fi

        if [[ ${#injected[@]} -gt 0 ]]; then
            printf '## Injected (overlay)\n'
            for n in "${injected[@]}"; do
                printf -- '- `%s` — (injected)\n' "$n"
            done
            printf '\n'
        fi

        # Inherited via `requires:` — populated by activate.sh via the
        # SCIAGENT_INHERITED env var ("<skill>=<parent>;..." pairs).
        if [[ -n "${SCIAGENT_INHERITED:-}" ]]; then
            printf '## Skills (inherited via requires:)\n'
            local pair sk parent
            IFS=';' read -ra _rb_inh <<< "$SCIAGENT_INHERITED"
            for pair in "${_rb_inh[@]+"${_rb_inh[@]}"}"; do
                [[ -z "$pair" ]] && continue
                sk="${pair%%=*}"
                parent="${pair#*=}"
                printf -- '- `%s` — (via `%s`)\n' "$sk" "$parent"
            done
            printf '\n'
        fi

        if [[ ${#_rb_agent_order[@]} -gt 0 ]]; then
            printf '## Sub-agents (effective, Claude-only)\n'
            for n in "${_rb_agent_order[@]}"; do
                s=""
                [[ -n "${_rb_agent_shadows[$n]:-}" ]] && s=" [shadows ${_rb_agent_shadows[$n]}]"
                printf -- '- `%s` — (%s)%s\n' "$n" "${_rb_agents[$n]}" "$s"
            done
            printf '\n'
        fi

        if [[ ${#_rb_command_order[@]} -gt 0 ]]; then
            printf '## Slash commands (effective, Claude-only)\n'
            for n in "${_rb_command_order[@]}"; do
                s=""
                [[ -n "${_rb_command_shadows[$n]:-}" ]] && s=" [shadows ${_rb_command_shadows[$n]}]"
                printf -- '- `/%s` — (%s)%s\n' "$n" "${_rb_commands[$n]}" "$s"
            done
            printf '\n'
        fi

        if [[ -n "$_rb_style" ]]; then
            printf '## Output style\n'
            printf -- '- `%s` — (%s)\n' "$_rb_style" "$_rb_style_role"
        fi
    }
}
