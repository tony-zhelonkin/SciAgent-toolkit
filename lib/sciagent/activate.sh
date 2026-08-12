# lib/sciagent/activate.sh — sciagent activate <base> [overlay] [--output-style <name>]
# Computes the effective merged stack (last-wins on name collisions), creates
# dual symlinks, rewrites the AGENTS.md managed block, writes manifest.
# Stack-walking and block-body rendering are delegated to stack.sh.
#
# Output style (Phase 5d) is no longer role-scoped — it is a SELECTION (one
# style exists on disk today), not a filter, so roles cannot express it via
# provenance the way skills/agents/commands do. Resolution precedence:
#   1. --output-style <name> flag (this invocation)
#   2. craft.yaml's `output_style:` key (toolkit-wide default)
#   3. none (no style mounted, no outputStyle written to settings.local.json)

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
    # Parse --output-style anywhere in the args; everything else is positional
    # (base [overlay]). --output-style wins over craft.yaml's output_style:
    # key when both are given; omitting both mounts no style at all.
    local -a _pos=()
    local OUTPUT_STYLE_FLAG=""
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -h|--help)
                # Deliberately the FIRST case, and it returns before the arg
                # loop can finish: everything cmd_activate mutates happens
                # after this loop, so `activate --help` mounts nothing. Without
                # this branch `-h` fell through to `_pos` and was read as a
                # role name ("role not found: -h").
                cat <<'USAGE'
sciagent activate <base> [overlay] [--output-style <name>]
  Mount the toolkit's whole catalog of skills, agents and commands into this
  project's .claude/ and .agents/, write the SCIAGENT:ROLES and SCIAGENT:CRAFT
  blocks in AGENTS.md, and record the result in .sciagent/manifest.json.

  Roles do NOT filter what is mounted — the entire catalog is mounted either
  way. A role decides provenance only: which role a mounted name is credited
  to in `sciagent status`, and which of the two stack slots wins a name
  collision (last wins). Anything no role names is credited to `catalog`.

  <base>             Base role (roles/<base>.yaml). Required.
  [overlay]          Optional second role. Stack depth is capped at 2.
  --output-style <n> Claude system prompt to mount/apply. Not role-scoped:
                     this flag wins, else craft.yaml's output_style:, else
                     none is mounted.

  Runs `sciagent validate --quiet` as a pre-flight and aborts before any
  mutation if it fails. Re-running is idempotent.

Exit code: 0 on success; 1 on a usage error, an unknown role, a stack deeper
than 2, a failed pre-flight, or a refused external-toolkit mutation.
USAGE
                return 0 ;;
            --output-style)
                if [[ -z "${2:-}" ]]; then
                    echo "usage: sciagent activate <base> [overlay] [--output-style <name>]" >&2
                    return 1
                fi
                OUTPUT_STYLE_FLAG="$2"
                shift 2 ;;
            *)
                _pos+=("$1")
                shift ;;
        esac
    done
    set -- "${_pos[@]+"${_pos[@]}"}"

    if [[ $# -lt 1 ]]; then
        echo "usage: sciagent activate <base> [overlay] [--output-style <name>]" >&2
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

    # Ensure the project-level Claude Code harness files exist, independent of
    # which role/stack is being activated. This is the fix for the toolkit
    # being vendored everywhere but statusline.sh/settings.json landing almost
    # nowhere: `new project` only rendered these at scaffold time, while
    # `activate` — the verb actually run against pre-existing / already-vendored
    # projects — never touched them. Non-clobbering: see claude_settings.sh.
    claude_settings_ensure_statusline
    claude_settings_ensure_project_defaults
    # settings.json registers hooks by path; materialize the bodies too, or the
    # enforcement tier is registered-but-absent (see claude_settings.sh).
    claude_settings_ensure_hooks

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

    # Resolve the output style (Phase 5d: no longer role-scoped). Precedence:
    # --output-style flag > craft.yaml's `output_style:` key > none.
    local OUTPUT_STYLE="" OUTPUT_STYLE_SOURCE=""
    if [[ -n "$OUTPUT_STYLE_FLAG" ]]; then
        OUTPUT_STYLE="$OUTPUT_STYLE_FLAG"
        OUTPUT_STYLE_SOURCE="--output-style flag"
    else
        local _craft_yaml _craft_style
        _craft_yaml=$(_craft_yaml_path)
        if [[ -f "$_craft_yaml" ]]; then
            _craft_style=$(role_scalar "$_craft_yaml" output_style)
            if [[ -n "$_craft_style" ]]; then
                OUTPUT_STYLE="$_craft_style"
                OUTPUT_STYLE_SOURCE="craft.yaml"
            fi
        fi
    fi

    # Resolve output_style → system-prompts/<file>.md by frontmatter name.
    # Validate before any filesystem mutation so a drifted request leaves
    # the project untouched (no half-state).
    local STYLE_SRC=""
    if [[ -n "$OUTPUT_STYLE" ]]; then
        STYLE_SRC=$(system_prompt_path "$OUTPUT_STYLE") || true
        if [[ -z "$STYLE_SRC" ]]; then
            echo "sciagent: $OUTPUT_STYLE_SOURCE requests output_style '$OUTPUT_STYLE'," >&2
            echo "  but no file in system-prompts/ has frontmatter 'name: $OUTPUT_STYLE'." >&2
            echo "  Available styles (frontmatter name → file):" >&2
            local pname pfile
            while IFS=$'\t' read -r pname pfile; do
                printf '    %s\t(%s)\n' "$pname" "$pfile" >&2
            done < <(system_prompt_inventory)
            return 1
        fi
    fi

    # Auto-deactivate if a stack is already active. claude_settings_teardown
    # must run BEFORE symlink_teardown_all rmdirs .sciagent (the state file
    # lives there) and BEFORE the new symlinks/block land so that the
    # settings.local.json revert is observable for the duration of the
    # re-activation rather than racing against the new apply.
    # Deferred until AFTER skill_resolve_transitive succeeds, so a failed
    # resolution leaves the previous stack intact.
    if manifest_exists; then
        claude_settings_teardown
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

    # Create symlinks. Canonical source paths are resolved via
    # resolve_canonical so that subfolders under agents/ and commands/ are
    # transparent to the consumer (symlinks stay flat in .claude/ and .agents/).
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

    local STYLE_APPLIED_TAG=""
    if [[ -n "$OUTPUT_STYLE" ]]; then
        symlink_create_dual output-styles "$OUTPUT_STYLE" "$STYLE_SRC"
        # Make the style Claude's active one by setting outputStyle in
        # .claude/settings.local.json. Symlinking alone makes the file
        # visible but does not select it as the active style.
        STYLE_APPLIED_TAG=$(claude_settings_apply "$OUTPUT_STYLE")
    fi

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
    if [[ -n "$OUTPUT_STYLE" ]]; then
        echo "  output_style: $OUTPUT_STYLE ($OUTPUT_STYLE_SOURCE) [settings.local.json: $STYLE_APPLIED_TAG]"
    fi
}
