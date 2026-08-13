# lib/sciagent/collisions.sh — cross-namespace name-collision detection.
#
# A "collision" is a single basename that exists in two or more of the four
# canonical namespaces (skills, agents, commands, roles). Mounting both is
# always legal — the symlink trees never overlap — but the human-facing
# disambiguation (`/foo` the command vs. `@foo` the agent vs. `foo` the
# skill vs. `foo` the role) becomes load-bearing. Most collisions in
# my practice are deliberate "family overlaps" (a slash command that dispatches
# the same-named agent or invokes the same-named skill); the CI allowlist
# at tests/collision-allowlist.txt is the source of truth for which ones
# are approved.
#
# `validate.sh` emits a soft warning per collision.
#
# The source namespaces are flat: direct skill directories and direct files
# under agents/, commands/, and roles/.

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
#   skills    — skills/<name>/SKILL.md (flat)
#   agents    — agents/<name>.md (flat)
#   commands  — commands/<name>.md (flat)
#   roles     — roles/<name>.yaml (flat)
collisions_enumerate() {
    local tk_root
    tk_root="$(_collisions_toolkit_root)"
    [[ -d "$tk_root" ]] || return 0

    local tmp
    tmp=$(mktemp)
    # Emit "<name> <kind>" lines; sort + group below.
    {
        # Skills: flat directories at depth 1.
        if [[ -d "$tk_root/skills" ]]; then
            local d name
            for d in "$tk_root"/skills/*/; do
                [[ -d "$d" ]] || continue
                name="$(basename "$d")"
                [[ "$name" == _* ]] && continue
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

# collisions_for_name <name>
# Emits the kind-csv (e.g. "skill,command") iff <name> collides across
# ≥2 namespaces in the toolkit; empty output + return 1 otherwise.
# Reuses collisions_enumerate rather than re-walking the four namespaces.
collisions_for_name() {
    local name="$1" nm csv
    while IFS=$'\t' read -r nm csv; do
        if [[ "$nm" == "$name" ]]; then
            printf '%s\n' "$csv"
            return 0
        fi
    done < <(collisions_enumerate)
    return 1
}

# collisions_is_allowlisted <name> <kinds-csv>
# Returns 0 when tests/collision-allowlist.txt records this exact (name, kinds)
# pair, 1 otherwise. Single parser shared by inject (the inline check) and
# validate (the soft-warn). A kinds-csv that differs from the allowlist line is
# treated as NOT allowlisted (callers surface it as an unexpected collision).
collisions_is_allowlisted() {
    local want_name="$1" want_kinds="$2"
    local allowlist
    allowlist="$(_collisions_toolkit_root)/tests/collision-allowlist.txt"
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
