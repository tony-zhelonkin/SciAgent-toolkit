# lib/sciagent/harness.sh — reusable AI-harness detection + path map.
#
# Five target harnesses: claude, pi, codex, agy, opencode. A harness is
# "present" when its binary is on PATH OR its config home exists on disk.
# This generalises the ad-hoc `.pi/` / `.claude/` probes in status.sh into a
# pure, side-effect-free API that `provision` (and future adapters) consume.
#
# Accessors are case-based (not a source-time assoc array) so they re-read
# $HOME / $CLAUDE_CONFIG_DIR / $CODEX_HOME at call time — a scratch-HOME
# self-test or a mid-session env change resolves correctly.
#
# Config homes:
#   claude   ${CLAUDE_CONFIG_DIR:-$HOME/.claude}
#   pi       $HOME/.pi
#   codex    ${CODEX_HOME:-$HOME/.codex}
#   agy      $HOME/.gemini/antigravity-cli   (also accepts $HOME/.local/bin/agy)
#   opencode $HOME/.config/opencode
#
# Global AGENTS.md-class context paths (used by `provision --context`):
#   claude   ${CLAUDE_CONFIG_DIR:-$HOME/.claude}/CLAUDE.md   (shim; imports @AGENTS.md-class)
#   pi       $HOME/.pi/agent/AGENTS.md
#   codex    ${CODEX_HOME:-$HOME/.codex}/AGENTS.md
#   agy      $HOME/.gemini/antigravity-cli/AGENTS.md
#   opencode $HOME/.config/opencode/AGENTS.md

# shellcheck shell=bash

# Canonical enumeration order for the five harnesses.
_HARNESS_ALL="claude pi codex agy opencode"

# harness_bin <name> — print the binary name probed on PATH. Non-zero for an
# unknown harness (also the membership test used to validate --harness input).
harness_bin() {
    case "$1" in
        claude)   printf 'claude' ;;
        pi)       printf 'pi' ;;
        codex)    printf 'codex' ;;
        agy)      printf 'agy' ;;
        opencode) printf 'opencode' ;;
        *) return 1 ;;
    esac
}

# harness_config_dir <name> — print the harness's config home.
harness_config_dir() {
    case "$1" in
        claude)   printf '%s' "${CLAUDE_CONFIG_DIR:-$HOME/.claude}" ;;
        pi)       printf '%s' "$HOME/.pi" ;;
        codex)    printf '%s' "${CODEX_HOME:-$HOME/.codex}" ;;
        agy)      printf '%s' "$HOME/.gemini/antigravity-cli" ;;
        opencode) printf '%s' "$HOME/.config/opencode" ;;
        *) return 1 ;;
    esac
}

# harness_global_context_path <name> — print the global AGENTS.md-class file
# the harness reads at user scope.
harness_global_context_path() {
    case "$1" in
        claude)   printf '%s' "${CLAUDE_CONFIG_DIR:-$HOME/.claude}/CLAUDE.md" ;;
        pi)       printf '%s' "$HOME/.pi/agent/AGENTS.md" ;;
        codex)    printf '%s' "${CODEX_HOME:-$HOME/.codex}/AGENTS.md" ;;
        agy)      printf '%s' "$HOME/.gemini/antigravity-cli/AGENTS.md" ;;
        opencode) printf '%s' "$HOME/.config/opencode/AGENTS.md" ;;
        *) return 1 ;;
    esac
}

# harness_is_present <name>
#   0 = present (binary on PATH or config home exists), 1 = absent,
#   2 = unknown harness name. Pure: no filesystem writes.
harness_is_present() {
    local name="$1" bin dir
    bin=$(harness_bin "$name") || return 2
    command -v "$bin" >/dev/null 2>&1 && return 0
    dir=$(harness_config_dir "$name")
    [[ -d "$dir" ]] && return 0
    # agy also ships as a user-local binary outside its config home.
    if [[ "$name" == agy && -x "$HOME/.local/bin/agy" ]]; then
        return 0
    fi
    return 1
}

# harness_detect — print each present harness name, one per line, in the
# canonical order. Empty output means none of the five are installed.
harness_detect() {
    local name
    for name in $_HARNESS_ALL; do
        harness_is_present "$name" && printf '%s\n' "$name"
    done
}
