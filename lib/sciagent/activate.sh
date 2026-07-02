# lib/sciagent/activate.sh — sciagent activate <base> [overlay]
# Computes the effective merged stack (last-wins on name collisions), creates
# dual symlinks, rewrites the AGENTS.md managed block, writes manifest.
# Stack-walking and block-body rendering are delegated to stack.sh.
#
# Complementary-skills warning surface:
#   Missing complementary-skills references are soft-warn only. Buffered
#   during the walk and emitted as one STDERR summary block after the
#   normal activation output. Exit code remains 0.

# shellcheck shell=bash

# ---------------------------------------------------------------------------
# _activate_read_complementary <skill-name>
# Emit complementary-skills entries, one per line, from frontmatter.
# Handles the same block-list format used by skill_read_requires.
# ---------------------------------------------------------------------------
_activate_read_complementary() {
    local name="$1"
    local file
    file=$(skill_frontmatter_path "$name") || return 0   # missing SKILL.md: skip
    _skill_extract_frontmatter "$file" | awk '
        BEGIN { in_meta=0; in_comp=0 }
        /^[A-Za-z_][A-Za-z0-9_-]*:/ {
            in_meta = ($0 ~ /^metadata:[ \t]*(#.*)?$/)
            in_comp = 0
            next
        }
        in_meta && /^[ \t]+complementary-skills:/ {
            line = $0
            sub(/^[ \t]+complementary-skills:[ \t]*/, "", line)
            sub(/[ \t]*#.*$/, "", line)
            sub(/[ \t]+$/, "", line)
            if (line ~ /^\[.*\]$/) {
                gsub(/^\[[ \t]*/, "", line)
                gsub(/[ \t]*\]$/, "", line)
                n = split(line, parts, /[ \t]*,[ \t]*/)
                for (i = 1; i <= n; i++) {
                    item = parts[i]
                    gsub(/^["'"'"']|["'"'"']$/, "", item)
                    if (item != "") print item
                }
                in_comp = 0
                next
            }
            in_comp = 1
            next
        }
        in_comp && /^[ \t]+[A-Za-z_][A-Za-z0-9_-]*:/ { in_comp=0; next }
        in_comp && /^[ \t]+-[ \t]+/ {
            item = $0
            sub(/^[ \t]+-[ \t]+/, "", item)
            sub(/[ \t]*#.*$/, "", item)
            sub(/[ \t]+$/, "", item)
            gsub(/^["'"'"']|["'"'"']$/, "", item)
            if (item != "") print item
        }
    '
}

cmd_activate() {
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

    # Pre-flight: run tag-vocab + requires-graph checks before any mutation.
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

    # Phase A: gather direct entries via stack_walk, preserving insertion order.
    # No mutation yet — all validation/resolution must succeed before we touch
    # the filesystem.
    local -a SKILL_ORDER=() AGENT_ORDER=() COMMAND_ORDER=()
    declare -A SKILLS=() AGENTS_M=() COMMANDS_M=()
    local OUTPUT_STYLE="" OUTPUT_STYLE_ROLE=""

    local kind name provider _shadow
    while IFS=$'\t' read -r kind name provider _shadow; do
        case "$kind" in
            SKILL)   SKILLS[$name]="$provider";   SKILL_ORDER+=("$name") ;;
            AGENT)   AGENTS_M[$name]="$provider";  AGENT_ORDER+=("$name") ;;
            COMMAND) COMMANDS_M[$name]="$provider"; COMMAND_ORDER+=("$name") ;;
            STYLE)   OUTPUT_STYLE="$name"; OUTPUT_STYLE_ROLE="$provider" ;;
        esac
    done < <(stack_walk "$base" "$overlay")

    # Phase B: transitive `requires:` resolution. For each direct skill, fold
    # its closure into SKILL_ORDER; new entries get provider=":requires:<parent>"
    # so the manifest is self-describing.
    # Resolver failures (cycles, missing targets) abort before any mutation.
    local direct
    for direct in "${SKILL_ORDER[@]+"${SKILL_ORDER[@]}"}"; do
        # Skip skills with no SKILL.md frontmatter (e.g. test fixtures).
        if ! skill_frontmatter_path "$direct" >/dev/null 2>&1; then
            continue
        fi
        # Capture closure via command-substitution so the resolver's exit
        # status propagates (process substitution + `done` would swallow it).
        local closure_out
        if ! closure_out=$(skill_resolve_transitive "$direct"); then
            echo "sciagent: aborting activation; requires resolution failed for '$direct'" >&2
            return 1
        fi
        local dep
        while IFS= read -r dep; do
            [[ -z "$dep" ]] && continue
            # Skip self and skills already in the effective set.
            [[ "$dep" == "$direct" ]] && continue
            [[ -n "${SKILLS[$dep]:-}" ]] && continue
            SKILLS[$dep]=":requires:$direct"
            SKILL_ORDER+=("$dep")
        done <<< "$closure_out"
    done

    # Phase B.5: scan complementary-skills for unresolvable references.
    # complementary-skills misses are soft-warn only — buffer here, 
    # emit at end-of-activate on STDERR.
    # Exit code stays 0 regardless of how many warnings accumulate.
    local -a _comp_warnings=()
    local sk
    for sk in "${SKILL_ORDER[@]+"${SKILL_ORDER[@]}"}"; do
        local cs
        while IFS= read -r cs; do
            [[ -z "$cs" ]] && continue
            # A complementary-skill is "missing" when it has no directory
            # under skills/. Use skill_frontmatter_path as the canonical check.
            if ! skill_frontmatter_path "$cs" >/dev/null 2>&1; then
                _comp_warnings+=("$sk references complementary-skill '$cs' — not present in skills/")
            fi
        done < <(_activate_read_complementary "$sk")
    done

    # Resolve output_style → system-prompts/<file>.md by frontmatter name.
    # Validate before any filesystem mutation so a drifted role spec leaves
    # the project untouched (no half-state).
    local STYLE_SRC=""
    if [[ -n "$OUTPUT_STYLE" ]]; then
        STYLE_SRC=$(system_prompt_path "$OUTPUT_STYLE") || true
        if [[ -z "$STYLE_SRC" ]]; then
            echo "sciagent: role '$OUTPUT_STYLE_ROLE' requests output_style '$OUTPUT_STYLE'," >&2
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
    #
    # Inject lifecycle note: activate is clean-slate w.r.t. injected entries
    # The previous manifest's injected rows are about to be torn down with 
    # the rest of the stack. Surface them on STDERR before teardown so 
    # the loss is attributed to this activate, not buried under the 
    # post-activation summary line. Exit code stays 0.
    if manifest_exists; then
        local -a _dropped_injected=()
        local _ov _sk _via _kind
        while IFS='|' read -r _ov _sk _via _kind; do
            [[ -z "$_sk" ]] && continue
            _dropped_injected+=("${_kind:-skill} $_sk")
        done < <(manifest_injected 2>/dev/null || true)
        if [[ "${#_dropped_injected[@]}" -gt 0 ]]; then
            {
                echo "sciagent: warning — activate is a clean-slate operation; dropping injected entries:"
                local _d
                for _d in "${_dropped_injected[@]}"; do
                    echo "  - $_d"
                done
                echo "  to preserve, run 'sciagent deactivate' first and re-inject after."
            } >&2
        fi
        claude_settings_teardown
        symlink_teardown_all
        block_remove AGENTS.md 2>/dev/null || true
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

    # Render and write the managed block. Inherited (`:requires:`) skills are
    # passed via the SCIAGENT_INHERITED env var so render_block_body can put
    # them under a dedicated subsection.
    local SCIAGENT_INHERITED=""
    local sk
    for sk in "${SKILL_ORDER[@]+"${SKILL_ORDER[@]}"}"; do
        if [[ "${SKILLS[$sk]:-}" == :requires:* ]]; then
            SCIAGENT_INHERITED+="${sk}=${SKILLS[$sk]#:requires:};"
        fi
    done
    export SCIAGENT_INHERITED
    local body
    body=$(render_block_body "$base" "$overlay")
    unset SCIAGENT_INHERITED
    block_write AGENTS.md "$body" || {
        echo "sciagent activate: failed to write managed block" >&2
        return 1
    }
    # Render the toolkit-owned CRAFT block (standing craft conventions) next to
    # ROLES. No-op when the toolkit ships no craft.yaml.
    craft_render_and_write AGENTS.md || {
        echo "sciagent activate: failed to write CRAFT block" >&2
        return 1
    }
    manifest_finalize "$(block_stored_hash AGENTS.md)" || {
        echo "sciagent activate: failed to finalize manifest" >&2
        return 1
    }

    # Summary.
    echo "Activated stack: $stack"
    echo "  skills:   ${#SKILL_ORDER[@]}"
    echo "  agents:   ${#AGENT_ORDER[@]}"
    echo "  commands: ${#COMMAND_ORDER[@]}"
    if [[ -n "$OUTPUT_STYLE" ]]; then
        echo "  output_style: $OUTPUT_STYLE ($OUTPUT_STYLE_ROLE) [settings.local.json: $STYLE_APPLIED_TAG]"
    fi

    # Emit complementary-skills warning block to STDERR after the activation
    # summary. One block, one place to look. Exit code stays 0.
    if [[ "${#_comp_warnings[@]}" -gt 0 ]]; then
        {
            echo ""
            echo "sciagent: activation completed with warnings:"
            echo ""
            echo "  complementary-skills (not installed):"
            local w
            for w in "${_comp_warnings[@]}"; do
                echo "    - $w"
            done
            echo ""
            echo "  These skills are optional companions. Install them via 'sciagent inject'"
            echo "  or add them to a role if the gap affects your current work."
        } >&2
    fi
}
