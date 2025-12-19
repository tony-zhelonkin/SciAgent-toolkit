#!/usr/bin/env bash
#
# validate-secrets.sh - Validates no API keys are exposed to git
#
# This script checks that all files potentially containing API keys are:
# 1. Properly gitignored
# 2. Not staged for commit
#
# Usage:
#   ./scripts/validate-secrets.sh [PROJECT_DIR]
#
# Exit codes:
#   0 - All checks passed
#   1 - Security issues found
#
# Can be used as a pre-commit hook:
#   ln -s ../../01_modules/SciAgent-toolkit/scripts/validate-secrets.sh .git/hooks/pre-commit

set -euo pipefail

# Source common utilities
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -f "${SCRIPT_DIR}/common.sh" ]; then
    source "${SCRIPT_DIR}/common.sh"
else
    # Fallback if common.sh not available
    log_info()  { echo "[INFO] $*"; }
    log_ok()    { echo "[OK] $*"; }
    log_warn()  { echo "[WARN] $*"; }
    log_error() { echo "[ERROR] $*"; }
fi

# Configuration
PROJECT_DIR="${1:-$(pwd)}"

# Files that may contain API keys and MUST be gitignored
FILES_TO_CHECK=(
    ".mcp.json"
    ".gemini/settings.json"
    ".claude/settings.local.json"
    ".env"
    ".devcontainer/.env"
    ".devcontainer/.env.local"
)

# Regex patterns for detecting API keys
# These patterns match common API key formats
KEY_PATTERNS=(
    'sk-proj-[A-Za-z0-9_-]{20,}'           # OpenAI project keys
    'sk-[A-Za-z0-9_-]{40,}'                 # OpenAI legacy keys
    'AIzaSy[A-Za-z0-9_-]{30,}'              # Google/Gemini API keys
    'sk-or-v1-[a-f0-9]{64}'                 # OpenRouter keys
    'xai-[A-Za-z0-9_-]{20,}'                # X.AI/Grok keys
    'ghp_[A-Za-z0-9]{36}'                   # GitHub personal access tokens
    'gho_[A-Za-z0-9]{36}'                   # GitHub OAuth tokens
)

# Combine patterns for grep
COMBINED_PATTERN=$(IFS='|'; echo "${KEY_PATTERNS[*]}")

errors=0
warnings=0

log_info "Validating secrets protection in: ${PROJECT_DIR}"
echo ""

# Check if we're in a git repository
if ! git -C "${PROJECT_DIR}" rev-parse --git-dir &>/dev/null; then
    log_warn "Not a git repository. Skipping gitignore checks."
    exit 0
fi

# --------------------------
# Check 1: Gitignore Status
# --------------------------
log_info "Checking gitignore status for sensitive files..."

for file in "${FILES_TO_CHECK[@]}"; do
    filepath="${PROJECT_DIR}/${file}"

    if [ -f "$filepath" ]; then
        # Check if file is gitignored
        if git -C "${PROJECT_DIR}" check-ignore -q "$filepath" 2>/dev/null; then
            log_ok "$file is gitignored"
        else
            log_error "$file exists but is NOT gitignored!"
            log_error "  -> Add '$file' to .gitignore immediately"
            errors=$((errors + 1))
        fi
    fi
done

echo ""

# --------------------------
# Check 2: Staged Files
# --------------------------
log_info "Checking for API keys in staged files..."

# Get list of staged files
staged_files=$(git -C "${PROJECT_DIR}" diff --cached --name-only 2>/dev/null || true)

if [ -n "$staged_files" ]; then
    while IFS= read -r staged_file; do
        staged_path="${PROJECT_DIR}/${staged_file}"
        if [ -f "$staged_path" ]; then
            # Check for API key patterns in staged content
            if git -C "${PROJECT_DIR}" diff --cached "$staged_file" 2>/dev/null | grep -qE "$COMBINED_PATTERN"; then
                log_error "Staged file '$staged_file' appears to contain API keys!"
                log_error "  -> Unstage with: git reset HEAD $staged_file"
                errors=$((errors + 1))
            fi
        fi
    done <<< "$staged_files"
fi

log_ok "No API keys detected in staged changes"
echo ""

# --------------------------
# Check 3: Content Scan
# --------------------------
log_info "Scanning sensitive files for API key patterns..."

for file in "${FILES_TO_CHECK[@]}"; do
    filepath="${PROJECT_DIR}/${file}"

    if [ -f "$filepath" ]; then
        # Check for actual key patterns (not placeholders)
        if grep -qE "$COMBINED_PATTERN" "$filepath" 2>/dev/null; then
            # Verify it's not a placeholder like ${GEMINI_API_KEY}
            if grep -E "$COMBINED_PATTERN" "$filepath" 2>/dev/null | grep -qvE '\$\{[^}]+\}'; then
                log_warn "$file contains what appears to be actual API keys"
                warnings=$((warnings + 1))

                # Check if gitignored (this is acceptable if gitignored)
                if git -C "${PROJECT_DIR}" check-ignore -q "$filepath" 2>/dev/null; then
                    log_info "  -> File is gitignored (acceptable)"
                else
                    log_error "  -> File is NOT gitignored! This is a security risk!"
                    errors=$((errors + 1))
                fi
            fi
        fi
    fi
done

echo ""

# --------------------------
# Check 4: Template Validation
# --------------------------
log_info "Validating templates use placeholders (not actual keys)..."

TEMPLATE_DIRS=(
    "${PROJECT_DIR}/01_modules/SciAgent-toolkit/templates"
    "${SCRIPT_DIR}/../templates"
)

for template_dir in "${TEMPLATE_DIRS[@]}"; do
    if [ -d "$template_dir" ]; then
        # Find all template files
        while IFS= read -r -d '' template_file; do
            if grep -qE "$COMBINED_PATTERN" "$template_file" 2>/dev/null; then
                # Check if these are placeholders
                if grep -E "$COMBINED_PATTERN" "$template_file" 2>/dev/null | grep -qvE '\$\{[^}]+\}'; then
                    log_error "Template '$template_file' contains actual API keys!"
                    log_error "  -> Templates should only contain \${VAR} placeholders"
                    errors=$((errors + 1))
                fi
            fi
        done < <(find "$template_dir" -type f \( -name "*.json" -o -name "*.template" -o -name "*.yaml" \) -print0 2>/dev/null)
    fi
done

log_ok "Templates use proper placeholders"
echo ""

# --------------------------
# Summary
# --------------------------
echo "========================================"
if [ $errors -gt 0 ]; then
    log_error "FAILED: Found $errors security issue(s)"
    if [ $warnings -gt 0 ]; then
        log_warn "Also found $warnings warning(s)"
    fi
    echo ""
    log_info "To fix gitignore issues, ensure your .gitignore contains:"
    echo "  # AI CLI configurations (contain API keys)"
    echo "  .gemini/"
    echo "  .mcp.json"
    echo "  .claude/settings.local.json"
    echo "  .env"
    echo "  .devcontainer/.env"
    exit 1
elif [ $warnings -gt 0 ]; then
    log_warn "PASSED with $warnings warning(s)"
    log_info "Warnings are acceptable if files are properly gitignored"
    exit 0
else
    log_ok "PASSED: All security checks passed"
    exit 0
fi
