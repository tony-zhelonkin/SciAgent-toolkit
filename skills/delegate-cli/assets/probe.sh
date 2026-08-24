#!/usr/bin/env bash

set -eu

usage() {
    echo "usage: probe.sh --workdir DIR [--model MODEL] [--sandbox MODE] [--web] [--file-read-check]" >&2
}

fail() {
    code="$1"
    shift
    printf 'probe.sh: %s\n' "$*" >&2
    exit "$code"
}

workdir=
model=
sandbox=
want_web=0
file_read_check=0

while [ "$#" -gt 0 ]; do
    case "$1" in
        --workdir)
            [ "$#" -ge 2 ] || { usage; exit 2; }
            workdir="$2"
            shift 2
            ;;
        --model)
            [ "$#" -ge 2 ] || { usage; exit 2; }
            model="$2"
            shift 2
            ;;
        --sandbox)
            [ "$#" -ge 2 ] || { usage; exit 2; }
            sandbox="$2"
            shift 2
            ;;
        --web)
            want_web=1
            shift
            ;;
        --file-read-check)
            file_read_check=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            usage
            exit 2
            ;;
    esac
done

[ -n "$workdir" ] || { usage; exit 2; }
[ -d "$workdir" ] || fail 2 "workdir is not a directory: $workdir"

command -v codex >/dev/null 2>&1 || fail 10 "codex is not available on PATH"

set +e
codex_version=$(codex --version 2>/dev/null)
version_status=$?
codex_help=$(codex exec --help 2>&1)
help_status=$?
set -e

[ "$version_status" -eq 0 ] && [ -n "$codex_version" ] \
    || fail 11 "codex --version failed"
[ "$help_status" -eq 0 ] && [ -n "$codex_help" ] \
    || fail 11 "codex exec --help failed"

has_output_last_message=0
has_stdin_prompt=0
has_config=0
has_enable=0
has_bypass=0
has_search_flag=0
has_workdir=0
has_skip_git_check=0
has_model=0
has_sandbox=0

printf '%s\n' "$codex_help" | grep -Eq -- '--output-last-message([=[:space:]<,]|$)' \
    && has_output_last_message=1
printf '%s\n' "$codex_help" | grep -Eiq -- 'stdin' \
    && has_stdin_prompt=1
printf '%s\n' "$codex_help" | grep -Eq -- '(^|[[:space:],])-c([[:space:],<]|$)|--config([=[:space:]<,]|$)' \
    && has_config=1
printf '%s\n' "$codex_help" | grep -Eq -- '--enable([=[:space:]<,]|$)' \
    && has_enable=1
printf '%s\n' "$codex_help" | grep -Eq -- '--dangerously-bypass-approvals-and-sandbox([[:space:],]|$)' \
    && has_bypass=1
printf '%s\n' "$codex_help" | grep -Eq -- '(^|[[:space:],])--search([=[:space:]<,]|$)' \
    && has_search_flag=1
printf '%s\n' "$codex_help" | grep -Eq -- '(^|[[:space:],])-C([[:space:],<]|$)|--cd([=[:space:]<,]|$)' \
    && has_workdir=1
printf '%s\n' "$codex_help" | grep -Eq -- '--skip-git-repo-check([[:space:],]|$)' \
    && has_skip_git_check=1
printf '%s\n' "$codex_help" | grep -Eq -- '(^|[[:space:],])-m([[:space:],<]|$)|--model([=[:space:]<,]|$)' \
    && has_model=1
printf '%s\n' "$codex_help" | grep -Eq -- '(^|[[:space:],])-s([[:space:],<]|$)|--sandbox([=[:space:]<,]|$)' \
    && has_sandbox=1

web_mode=unsupported
if [ "$has_config" -eq 1 ] && [ "$has_enable" -eq 1 ]; then
    web_mode=config_enable
fi

printf 'CODEX_VERSION=%q\n' "$codex_version"
printf 'HAS_OUTPUT_LAST_MESSAGE=%s\n' "$has_output_last_message"
printf 'HAS_STDIN_PROMPT=%s\n' "$has_stdin_prompt"
printf 'HAS_CONFIG=%s\n' "$has_config"
printf 'HAS_ENABLE=%s\n' "$has_enable"
printf 'HAS_BYPASS=%s\n' "$has_bypass"
printf 'HAS_SEARCH_FLAG=%s\n' "$has_search_flag"
printf 'WEB_MODE=%s\n' "$web_mode"

[ "$has_output_last_message" -eq 1 ] \
    || fail 12 "codex exec lacks --output-last-message"
[ "$has_stdin_prompt" -eq 1 ] \
    || fail 12 "codex exec help does not document stdin prompts"
[ "$has_workdir" -eq 1 ] \
    || fail 12 "codex exec lacks -C/--cd"
[ "$has_skip_git_check" -eq 1 ] \
    || fail 12 "codex exec lacks --skip-git-repo-check"
if [ -n "$model" ] && [ "$has_model" -ne 1 ]; then
    fail 12 "codex exec lacks -m/--model"
fi
if [ -n "$sandbox" ] && [ "$has_sandbox" -ne 1 ]; then
    fail 12 "codex exec lacks -s/--sandbox"
fi
if [ "$want_web" -eq 1 ] && [ "$web_mode" != config_enable ]; then
    fail 13 "web search is unsupported by this codex exec"
fi

if [ "$file_read_check" -eq 1 ]; then
    check_dir=$(mktemp -d "${TMPDIR:-/tmp}/scio-delegate-probe.XXXXXX") \
        || fail 14 "could not create the file-read check directory"
    trap 'rm -rf "$check_dir"' EXIT

    nonce="SCIO_DELEGATE_FILE_READ_${$}_$(date +%s)"
    nonce_file="$check_dir/nonce.txt"
    prompt_file="$check_dir/prompt.md"
    final_file="$check_dir/final.md"
    stream_file="$check_dir/stream.log"
    printf '%s\n' "$nonce" > "$nonce_file"
    printf 'Read %s with a file tool and reply with only its exact contents.\n' "$nonce_file" \
        > "$prompt_file"

    codex_args=(exec -C "$workdir" --skip-git-repo-check -o "$final_file")
    [ -z "$model" ] || codex_args+=(-m "$model")
    [ -z "$sandbox" ] || codex_args+=(-s "$sandbox")

    if codex "${codex_args[@]}" - < "$prompt_file" > "$stream_file" 2>&1; then
        check_status=0
    else
        check_status=$?
    fi
    # A sandboxed codex cannot read files at all when the container blocks the
    # namespace or filter operations its sandbox is built on. That failure looks
    # like a model refusing to cooperate, so name it: the remedy is a container
    # policy change or an explicitly authorized --bypass, and no prompt wording
    # substitutes for either.
    if [ "$check_status" -ne 0 ]; then
        if grep -Eqi 'bwrap|landlock|seccomp|unshare|namespace|Operation not permitted' \
                "$stream_file" 2>/dev/null; then
            printf 'SANDBOX_FILE_READ=blocked\n'
            printf 'probe.sh: codex could not inspect a file with sandbox %s.\n' \
                "${sandbox:-default}" >&2
            printf 'probe.sh: the container blocks the operation its sandbox needs. Either relax the\n' >&2
            printf 'probe.sh: container policy (see scbio-docker docs/ai-integration.md) or obtain\n' >&2
            printf 'probe.sh: explicit authorization for --bypass, which removes the sandbox.\n' >&2
            exit 14
        fi
        fail 14 "codex could not inspect the nonce file"
    fi
    grep -Fq -- "$nonce" "$final_file" 2>/dev/null \
        || fail 14 "codex final output did not contain the nonce"
    printf 'SANDBOX_FILE_READ=ok\n'
fi
