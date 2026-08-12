# lib/sciagent/stack.sh — shared stack-walking and block-body rendering.
#
# Consumers: activate.sh, status.sh.
# The dispatcher (bin/sciagent) does not source this directly; verb modules do.

# shellcheck shell=bash

# _sw_record <map-name> <shadow-map-name> <order-array-name> <name> <role>
# Records <name>=<role> into the provider map with last-wins semantics.
# On a repeat name, the *previous* provider is appended to the shadow CSV
# (preserving the existing comma-joined accumulation order). On a first
# sighting, <name> is pushed onto the order array.
# Locals are `_`-prefixed so a nameref can never alias one of this helper's
# own locals (a known bash nameref footgun); the caller arrays are `_sw_*`.
_sw_record() {
    local -n _map="$1"
    local -n _shadow="$2"
    local -n _order="$3"
    local name="$4" role="$5"

    if [[ -n "${_map[$name]:-}" ]]; then
        local prev="${_map[$name]}"
        if [[ -n "${_shadow[$name]:-}" ]]; then
            _shadow[$name]="${_shadow[$name]},$prev"
        else
            _shadow[$name]="$prev"
        fi
    else
        _order+=("$name")
    fi
    _map[$name]="$role"
}

# stack_walk <role1> [role2]
# Emits the resolved effective table in tagged form, one entry per line:
#   SKILL    <name>  <providing-role>  <shadowed-roles-csv-or-empty>
#   AGENT    <name>  <providing-role>  <shadowed>
#   COMMAND  <name>  <providing-role>  <shadowed>
#
# Last-wins resolution is applied across the stack. Output order:
#   SKILL entries (insertion order), then AGENT, then COMMAND.
# The _injected synthetic role is silently skipped (no YAML file).
#
# NOTHING IS FILTERED BY ROLE. The whole catalog — skills, agents, commands —
# is always resolved; a role's `skills:`/`agents:`/`commands:` lists only
# affect PROVENANCE (which role an item is attributed to in the managed
# block), never visibility. Roles are pure provenance labels (Phase 5d).
#
# Why (skills, original rationale): subsetting the catalog saved ~6k tokens
# (~3% of a 200k context) while breaking cross-references — 20 of 43 mounted
# skills in Meta-Aging/14616-DM pointed at skills that were not mounted, so
# the agent was routed to skills it could not see. Harnesses also load only
# `name` + `description` per skill (~100 tokens) until a skill is activated,
# so the platform already does the token optimization the filter was invented
# for, and does it without severing the routing graph. See docs 07a-07d.
#
# Why (agents, commands — Phase 5d): unlike skills, agents/commands had no
# fallback and were still genuinely gated by role — `base` mounted 8/21
# agents and 3/24 commands. The preload cost of mounting the rest is ~250
# tokens for all 21 agent descriptions combined; commands carry no
# `description:` frontmatter at all, so their preload cost is 0 tokens. That
# is negligible next to the ~6k tokens the skill filter saved (and which was
# judged not worth severing the routing graph for), so there is no token
# argument for keeping agents/commands gated either. The owner's call: roles
# are obsolete and were never actively used as a curation mechanism, so
# de-gate them the same way skills already are.
stack_walk() {
    local base="$1"
    local overlay="${2:-}"

    declare -A _sw_skills=() _sw_agents=() _sw_commands=()
    declare -A _sw_skill_shadows=() _sw_agent_shadows=() _sw_command_shadows=()
    local -a _sw_skill_order=() _sw_agent_order=() _sw_command_order=()

    local role line kind name
    for role in "$base" ${overlay:+"$overlay"}; do
        [[ "$role" == "_injected" ]] && continue
        while IFS= read -r line; do
            kind="${line%% *}"
            name="${line#* }"
            case "$kind" in
                SKILL)
                    _sw_record _sw_skills   _sw_skill_shadows   _sw_skill_order   "$name" "$role" ;;
                AGENT)
                    _sw_record _sw_agents   _sw_agent_shadows   _sw_agent_order   "$name" "$role" ;;
                COMMAND)
                    _sw_record _sw_commands _sw_command_shadows _sw_command_order "$name" "$role" ;;
            esac
        done < <(role_load "$role")
    done

    # The catalog IS the skill set. Roles above supplied provenance for the
    # skills they name; every remaining skill is mounted too, attributed to
    # "catalog". Sorted so mount order is deterministic across hosts.
    local _sw_sk _sw_name
    while IFS= read -r _sw_sk; do
        _sw_name=$(basename "$(dirname "$_sw_sk")")
        case "$_sw_name" in
            _*) continue ;;   # _TEMPLATE, and any other underscore-prefixed scaffold
        esac
        [[ -n "${_sw_skills[$_sw_name]:-}" ]] && continue   # already attributed to a role
        _sw_record _sw_skills _sw_skill_shadows _sw_skill_order "$_sw_name" catalog
    done < <(find "$SCIAGENT_TOOLKIT/skills" -mindepth 2 -maxdepth 2 -name SKILL.md 2>/dev/null | sort)

    # AGENT/COMMAND catalog-fallback (Phase 5d): same treatment as SKILL above.
    # agents/**/*.md and commands/**/*.md are nested under a subdirectory
    # (e.g. agents/analysis-base/captions.md); the mounted NAME is the
    # filename without extension (matches resolve_canonical's by-filename
    # search and _list_agents/_list_commands), not any frontmatter field.
    # README.md files and any path with an underscore-prefixed segment
    # (scaffolding/retired, e.g. commands/architect/_superseded/) are skipped.
    while IFS= read -r _sw_sk; do
        _sw_name=$(basename "$_sw_sk" .md)
        [[ "$_sw_name" == "README" ]] && continue
        case "$_sw_sk" in */_*) continue ;; esac
        [[ -n "${_sw_agents[$_sw_name]:-}" ]] && continue   # already attributed to a role
        _sw_record _sw_agents _sw_agent_shadows _sw_agent_order "$_sw_name" catalog
    done < <(find "$SCIAGENT_TOOLKIT/agents" -type f -name '*.md' 2>/dev/null | sort)

    while IFS= read -r _sw_sk; do
        _sw_name=$(basename "$_sw_sk" .md)
        [[ "$_sw_name" == "README" ]] && continue
        case "$_sw_sk" in */_*) continue ;; esac
        [[ -n "${_sw_commands[$_sw_name]:-}" ]] && continue   # already attributed to a role
        _sw_record _sw_commands _sw_command_shadows _sw_command_order "$_sw_name" catalog
    done < <(find "$SCIAGENT_TOOLKIT/commands" -type f -name '*.md' 2>/dev/null | sort)

    # Emit: SKILL entries first, then AGENT, then COMMAND.
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
}

# render_block_body <base> <overlay>
# Produces the body text that goes between the BEGIN/END managed-block markers
# per architecture §4. Reads the role YAMLs via stack_walk.
render_block_body() {
    local base="$1"
    local overlay="${2:-}"

    # Run stack walk and load results into arrays.
    declare -A _rb_skills=() _rb_agents=() _rb_commands=()
    declare -A _rb_skill_shadows=() _rb_agent_shadows=() _rb_command_shadows=()
    local -a _rb_skill_order=() _rb_agent_order=() _rb_command_order=()

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
    }
}
