# lib/sciagent/collisions.sh — cross-namespace name-collision detection.
#
# A "collision" is a single basename that exists in two or more of the four
# canonical namespaces (skills, agents, commands, roles). Mounting both is
# always legal — the symlink trees never overlap — but the human-facing
# disambiguation (`/foo` the command vs. `@foo` the agent vs. `foo` the
# skill vs. `foo` the role) becomes load-bearing. Most collisions in
# practice are deliberate "family overlaps" (a slash command that dispatches
# the same-named agent or invokes the same-named skill); the CI allowlist
# at tests/collision-allowlist.txt is the source of truth for which ones
# are blessed.
#
# Two callers consume this helper:
#   - validate.sh — emits a soft-warn per collision (exit 0 either way)
#   - status.sh   — surfaces collisions for the *active stack* in a Notes
#                   section, annotated against the allowlist
#
# The enumeration mirrors the recursive walks in inject.sh (agents/ and
# commands/ allow nested subdirs; skills/ is always flat at depth 1; roles/
# is always flat at depth 1).

# shellcheck shell=bash

# collisions_toolkit_root
# Resolve toolkit root the same way validate.sh does (env var first,
# this-file-relative fallback). Kept private to avoid bleeding helper state.
_collisions_toolkit_root() {
    if [[ -n "${SCIAGENT_TOOLKIT:-}" ]]; then
        printf '%s\n' "$SCIAGENT_TOOLKIT"
        return
    fi
    local self_dir
    self_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    (cd "$self_dir/../.." && pwd)
}

# collisions_enumerate
# Walks the four canonical dirs and prints one collision per line:
#   <name>\t<kind1>,<kind2>[,<kind3>[,<kind4>]]
# Kinds are emitted in fixed order: skill,agent,command,role. Names without
# at least two namespace matches are not emitted. Output is sorted by name.
#
# The four enumeration patterns:
#   skills    — skills/<name>/SKILL.md (flat; _TEMPLATE excluded)
#   agents    — agents/**/<name>.md (recursive; README.md and hidden dirs excluded)
#   commands  — commands/**/<name>.md (recursive; README.md and hidden dirs excluded)
#   roles     — roles/<name>.yaml (flat)
collisions_enumerate() {
    local tk_root
    tk_root="$(_collisions_toolkit_root)"
    [[ -d "$tk_root" ]] || return 0

    local tmp
    tmp=$(mktemp)
    # Emit "<name> <kind>" lines; sort + group below.
    {
        # Skills: flat directories at depth 1, excluding the scaffold template.
        if [[ -d "$tk_root/skills" ]]; then
            local d name
            for d in "$tk_root"/skills/*/; do
                [[ -d "$d" ]] || continue
                name="$(basename "$d")"
                [[ "$name" == "_TEMPLATE" ]] && continue
                [[ -f "$d/SKILL.md" ]] || continue
                printf '%s skill\n' "$name"
            done
        fi
        # Agents: recursive, but skip README.md and hidden dirs.
        if [[ -d "$tk_root/agents" ]]; then
            local f bn
            while IFS= read -r f; do
                [[ -z "$f" ]] && continue
                bn="$(basename "$f" .md)"
                [[ "$bn" == "README" ]] && continue
                printf '%s agent\n' "$bn"
            done < <(find "$tk_root/agents" -type f -name '*.md' -not -path '*/.*' 2>/dev/null)
        fi
        # Commands: recursive, same exclusions as agents.
        if [[ -d "$tk_root/commands" ]]; then
            local f bn
            while IFS= read -r f; do
                [[ -z "$f" ]] && continue
                bn="$(basename "$f" .md)"
                [[ "$bn" == "README" ]] && continue
                printf '%s command\n' "$bn"
            done < <(find "$tk_root/commands" -type f -name '*.md' -not -path '*/.*' 2>/dev/null)
        fi
        # Roles: flat *.yaml at depth 1.
        if [[ -d "$tk_root/roles" ]]; then
            local f bn
            for f in "$tk_root"/roles/*.yaml; do
                [[ -f "$f" ]] || continue
                bn="$(basename "$f" .yaml)"
                printf '%s role\n' "$bn"
            done
        fi
    } > "$tmp"

    # Group by name, keep names with >=2 entries, emit fixed-order kind csv.
    sort "$tmp" | awk '
        {
            name = $1
            kind = $2
            if (name != current) {
                if (current != "" && count >= 2) emit()
                current = name
                count = 0
                has_skill = 0; has_agent = 0; has_command = 0; has_role = 0
            }
            if (kind == "skill"   && !has_skill)   { has_skill   = 1; count++ }
            if (kind == "agent"   && !has_agent)   { has_agent   = 1; count++ }
            if (kind == "command" && !has_command) { has_command = 1; count++ }
            if (kind == "role"    && !has_role)    { has_role    = 1; count++ }
        }
        END { if (current != "" && count >= 2) emit() }
        function emit(   parts, n) {
            n = 0
            if (has_skill)   { n++; parts[n] = "skill"   }
            if (has_agent)   { n++; parts[n] = "agent"   }
            if (has_command) { n++; parts[n] = "command" }
            if (has_role)    { n++; parts[n] = "role"    }
            out = parts[1]
            for (i = 2; i <= n; i++) out = out "," parts[i]
            printf "%s\t%s\n", current, out
        }
    '
    rm -f "$tmp"
}
