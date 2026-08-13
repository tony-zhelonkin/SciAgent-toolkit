# lib/sciagent/activate.sh — sciagent activate <base> [overlay]
# Computes the effective merged stack (last-wins on name collisions), creates
# dual symlinks, rewrites the AGENTS.md managed block, writes manifest.
# Stack-walking and block-body rendering are delegated to stack.sh.

# shellcheck shell=bash

# ---------------------------------------------------------------------------
# ensure_claude_md_shim
# Guarantee the project CLAUDE.md carries the `@AGENTS.md` import line so Claude
# Code — which does not read AGENTS.md natively — picks up the canonical
# context. Idempotent + non-clobbering:
#   - absent      → create CLAUDE.md containing just `@AGENTS.md`
#   - present, has import → no-op
#   - present, no import  → prepend the import, preserving existing content
#
# Ownership record ($_CLAUDE_MD_STATE, plain text):
#   line 1: "created" (we wrote the whole file) or "prepended" (file
#           pre-existed; we added the import header in front of it)
#   line 2: sha1 — of the whole file as written (created), or of the
#           pre-existing content BEFORE we prepended (prepended)
# claude_md_shim_teardown (below) consults this exactly once to reverse
# precisely what we did, or to cede ownership with a warning if the file has
# since been edited/replaced. No state is written for the already-has-import
# no-op — there is nothing for teardown to reverse.
# ---------------------------------------------------------------------------
_CLAUDE_MD_STATE=".sciagent/claude_md.state"

ensure_claude_md_shim() {
    local f="CLAUDE.md"
    local import="@AGENTS.md"
    if [[ ! -f "$f" ]]; then
        printf '%s\n' "$import" > "$f" || {
            echo "sciagent activate: failed to write $f" >&2
            return 1
        }
        echo "wrote: $f (@AGENTS.md import shim)"
        mkdir -p "$(dirname "$_CLAUDE_MD_STATE")"
        printf 'created\n%s\n' "$(sciagent_sha1_file "$f")" > "$_CLAUDE_MD_STATE"
        return 0
    fi
    # Already imports AGENTS.md (a bare `@AGENTS.md` line) → nothing to do.
    if grep -qE '^[[:space:]]*@AGENTS\.md[[:space:]]*$' "$f"; then
        return 0
    fi
    # Present but missing the import — prepend it, preserving existing bytes.
    local orig_hash
    orig_hash=$(sciagent_sha1_file "$f")
    local tmp
    tmp=$(mktemp)
    printf '%s\n\n' "$import" > "$tmp"
    cat "$f" >> "$tmp"
    cp "$tmp" "$f"
    rm -f "$tmp"
    echo "updated: $f (added @AGENTS.md import shim)"
    mkdir -p "$(dirname "$_CLAUDE_MD_STATE")"
    printf 'prepended\n%s\n' "$orig_hash" > "$_CLAUDE_MD_STATE"
}

# claude_md_shim_teardown
# Reverse ensure_claude_md_shim using the saved state (see above):
#   created    → remove CLAUDE.md iff its content is unchanged since we wrote it.
#   prepended  → strip exactly the "@AGENTS.md\n\n" header we added, iff the
#                remainder still hashes to the pre-existing content we snapshotted;
#                otherwise leave the file untouched.
# Either way, mismatched content is left in place with a warning — never
# silently discarded — and the state file is removed regardless, so a second
# deactivate is a silent no-op (idempotent) rather than re-warning forever.
claude_md_shim_teardown() {
    [[ -f "$_CLAUDE_MD_STATE" ]] || return 0
    local tag stored f="CLAUDE.md"
    tag=$(sed -n '1p' "$_CLAUDE_MD_STATE")
    stored=$(sed -n '2p' "$_CLAUDE_MD_STATE")

    case "$tag" in
        created)
            if [[ -f "$f" ]]; then
                local cur
                cur=$(sciagent_sha1_file "$f")
                if [[ "$cur" == "$stored" ]]; then
                    rm -f "$f"
                else
                    echo "sciagent: warning — $f was modified since sciagent created it; leaving it in place" >&2
                fi
            fi
            ;;
        prepended)
            if [[ -f "$f" ]]; then
                local prefix=$'@AGENTS.md\n\n'
                local prefix_len=${#prefix}
                # Compare as byte streams via cmp, not `$(...)` — command
                # substitution strips trailing newlines, which would make
                # a 2-trailing-newline prefix ("@AGENTS.md\n\n") spuriously
                # mismatch the freshly-written file every time.
                if head -c "$prefix_len" "$f" | cmp -s - <(printf '%s' "$prefix"); then
                    local remainder_hash tmp
                    remainder_hash=$(tail -c "+$((prefix_len + 1))" "$f" | sciagent_sha1_stream)
                    if [[ "$remainder_hash" == "$stored" ]]; then
                        tmp=$(mktemp)
                        tail -c "+$((prefix_len + 1))" "$f" > "$tmp"
                        cp "$tmp" "$f"
                        rm -f "$tmp"
                    else
                        echo "sciagent: warning — $f content was modified since sciagent added the @AGENTS.md import; leaving it in place" >&2
                    fi
                else
                    echo "sciagent: warning — $f's @AGENTS.md import was modified since sciagent added it; leaving it in place" >&2
                fi
            fi
            ;;
    esac
    rm -f "$_CLAUDE_MD_STATE"
}

cmd_activate() {
    # Parse help anywhere in the args before any mutation; everything else is
    # positional (base [overlay]).
    local -a _pos=()
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -h|--help)
                # Deliberately the FIRST case, and it returns before the arg
                # loop can finish: everything cmd_activate mutates happens
                # after this loop, so `activate --help` mounts nothing. Without
                # this branch `-h` fell through to `_pos` and was read as a
                # role name ("role not found: -h").
                cat <<'USAGE'
sciagent activate <base> [overlay]
  Mount the toolkit's whole catalog of skills, agents and commands into this
  project's .claude/ and .agents/, write the SCIAGENT:ROLES and SCIAGENT:CRAFT
  blocks in AGENTS.md, and record the result in .sciagent/manifest.json.

  Roles do NOT filter what is mounted — the entire catalog is mounted either
  way. A role decides provenance only: which role a mounted name is credited
  to in `sciagent status`, and which of the two stack slots wins a name
  collision (last wins). Anything no role names is credited to `catalog`.

  <base>             Base role (roles/<base>.yaml). Required.
  [overlay]          Optional second role. Stack depth is capped at 2.

  Runs `sciagent validate --quiet` as a pre-flight and aborts before any
  mutation if it fails. Re-running is idempotent.

Exit code: 0 on success; 1 on a usage error, an unknown role, a stack deeper
than 2, a failed pre-flight, or a refused external-toolkit mutation.
USAGE
                return 0 ;;
            *)
                _pos+=("$1")
                shift ;;
        esac
    done
    set -- "${_pos[@]+"${_pos[@]}"}"

    if [[ $# -lt 1 ]]; then
        echo "usage: sciagent activate <base> [overlay]" >&2
        return 1
    fi
    if [[ $# -gt 2 ]]; then
        echo "Maximum stack depth is 2 (base + overlay)" >&2
        return 1
    fi

    local base="$1"
    local overlay="${2:-}"

    if ! role_exists "$base"; then
        echo "role not found: $base ($(role_path "$base"))" >&2
        return 1
    fi
    if [[ -n "$overlay" ]] && ! role_exists "$overlay"; then
        echo "role not found: $overlay ($(role_path "$overlay"))" >&2
        return 1
    fi

    # Pre-flight: check every skill's frontmatter shape before any mutation.
    # validate.sh must be sourced by the caller (bin/sciagent sources it for
    # the activate verb). Called with --quiet so success produces no output.
    if ! cmd_validate --quiet; then
        echo "sciagent: aborting activation; validate checks failed (see above)" >&2
        return 1
    fi

    # Materialize and register the two project guardrails.
    claude_settings_ensure_hooks || return 1

    # Phase A: gather direct entries via stack_walk, preserving insertion order.
    # No mutation yet — all validation/resolution must succeed before we touch
    # the filesystem.
    local -a SKILL_ORDER=() AGENT_ORDER=() COMMAND_ORDER=()
    declare -A SKILLS=() AGENTS_M=() COMMANDS_M=()

    local kind name provider _shadow
    while IFS=$'\t' read -r kind name provider _shadow; do
        case "$kind" in
            SKILL)   SKILLS[$name]="$provider";   SKILL_ORDER+=("$name") ;;
            AGENT)   AGENTS_M[$name]="$provider";  AGENT_ORDER+=("$name") ;;
            COMMAND) COMMANDS_M[$name]="$provider"; COMMAND_ORDER+=("$name") ;;
        esac
    done < <(stack_walk "$base" "$overlay")

    # Auto-deactivate if a stack is already active.
    # Deferred until AFTER skill_resolve_transitive succeeds, so a failed
    # resolution leaves the previous stack intact.
    if manifest_exists; then
        symlink_teardown_all
        # ROLES explicitly: this is the teardown-before-rewrite half of
        # re-activation, not a generic "clear the file" — the block_write
        # below re-renders the same id. CRAFT is torn down separately via
        # craft_remove (its own id).
        block_remove AGENTS.md ROLES 2>/dev/null || true
        craft_remove AGENTS.md 2>/dev/null || true
    fi

    # Start the manifest staging buffer before creating any symlinks.
    local stack="$base"
    [[ -n "$overlay" ]] && stack="$base $overlay"
    manifest_begin "$stack"

    # Create symlinks from the flat canonical catalog into both harness trees.
    local n src
    for n in "${SKILL_ORDER[@]:-}"; do
        [[ -z "$n" ]] && continue
        src=$(resolve_canonical skills "$n") || return 1
        symlink_create_dual skills "$n" "$src"
    done
    for n in "${AGENT_ORDER[@]:-}"; do
        [[ -z "$n" ]] && continue
        src=$(resolve_canonical agents "$n") || return 1
        symlink_create_dual agents "$n" "$src"
    done
    for n in "${COMMAND_ORDER[@]:-}"; do
        [[ -z "$n" ]] && continue
        src=$(resolve_canonical commands "$n") || return 1
        symlink_create_dual commands "$n" "$src"
    done
    # Create the helper-lib symlink for analysis-type repos (no-op if no 02_analysis/).
    symlink_create_helper_lib
    # ...and the importable shim modules beside those mounts. Same 02_analysis/
    # gate. Until 2026-08-12 only `new project` wrote these, so a project that
    # was already provisioned got the hyphenated mounts and nothing it could
    # import — the third instance of "written at scaffold time, never reaching
    # the field" after the hook bodies and the status line. Ownership discipline:
    # ownership.sh.
    helper_shims_ensure

    # Render and write the managed block.
    local body
    body=$(render_block_body "$base" "$overlay")
    block_write AGENTS.md "$body" ROLES || {
        echo "sciagent activate: failed to write managed block" >&2
        return 1
    }
    # Render the toolkit-owned CRAFT block (standing craft conventions) next to
    # ROLES. No-op when the toolkit ships no craft.yaml.
    craft_render_and_write AGENTS.md || {
        echo "sciagent activate: failed to write CRAFT block" >&2
        return 1
    }
    # AGENTS.md is the canonical context surface; Claude Code does not read it
    # natively, so guarantee the project CLAUDE.md = @AGENTS.md import shim.
    ensure_claude_md_shim
    manifest_finalize || {
        echo "sciagent activate: failed to finalize manifest" >&2
        return 1
    }

    # Summary.
    echo "Activated stack: $stack"
    echo "  skills:   ${#SKILL_ORDER[@]}"
    echo "  agents:   ${#AGENT_ORDER[@]}"
    echo "  commands: ${#COMMAND_ORDER[@]}"
}
