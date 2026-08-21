#!/usr/bin/env bash
# tests/test_delegate_cli_assets.sh — probe and launch delegation without the real CLI.

set -u
. "$(dirname "$0")/_lib.sh"

PROBE="$TOOLKIT_ROOT/skills/delegate-cli/assets/probe.sh"
LAUNCH="$TOOLKIT_ROOT/skills/delegate-cli/assets/launch.sh"

assert_file_exists "$PROBE" "delegate-cli probe asset missing"
assert_file_exists "$LAUNCH" "delegate-cli launch asset missing"
[ -x "$PROBE" ] || { echo "FAIL [$_TEST_NAME] probe.sh is not executable" >&2; exit 1; }
[ -x "$LAUNCH" ] || { echo "FAIL [$_TEST_NAME] launch.sh is not executable" >&2; exit 1; }
bash -n "$PROBE"
bash -n "$LAUNCH"

setup_tmpdir
mkdir -p "$TMPDIR_TEST/bin" "$TMPDIR_TEST/no-bin" "$TMPDIR_TEST/work" "$TMPDIR_TEST/runtime"

cat > "$TMPDIR_TEST/bin/codex" <<'EOF'
#!/usr/bin/env bash
set -u

if [ "${1:-}" = "--version" ]; then
    echo "codex-cli-exec 9.9.9-stub"
    exit 0
fi

if [ "${1:-}" = "exec" ] && [ "${2:-}" = "--help" ]; then
    cat <<'HELP'
Usage: codex exec [OPTIONS] [PROMPT]
Arguments:
  [PROMPT]  Use - to read the prompt from stdin
Options:
  -m, --model <MODEL>
  -s, --sandbox <MODE>
  -C, --cd <DIR>
      --skip-git-repo-check
  -o, --output-last-message <FILE>
      --dangerously-bypass-approvals-and-sandbox
HELP
    if [ "${SCIO_STUB_NO_WEB:-0}" -ne 1 ]; then
        cat <<'HELP'
  -c, --config <key=value>
      --enable <FEATURE>
HELP
    fi
    exit 0
fi

[ "${1:-}" = "exec" ] || exit 90
shift
output=
last_arg=
while [ "$#" -gt 0 ]; do
    case "$1" in
        -o|--output-last-message)
            output="$2"
            shift 2
            ;;
        -m|--model|-s|--sandbox|-C|--cd|-c|--config|--enable)
            shift 2
            ;;
        --skip-git-repo-check|--dangerously-bypass-approvals-and-sandbox)
            shift
            ;;
        *)
            last_arg="$1"
            shift
            ;;
    esac
done
[ "$last_arg" = "-" ] || exit 91
prompt_text=$(cat)

if [ -n "${SCIO_STUB_BLOCK:-}" ]; then
    while [ -e "$SCIO_STUB_BLOCK" ]; do
        sleep 0.02
    done
fi

[ -n "$output" ] || exit 92
if printf '%s\n' "$prompt_text" | grep -q '^Read .*nonce.txt with a file tool'; then
    nonce_path=$(printf '%s\n' "$prompt_text" | sed -n 's/^Read \(.*nonce.txt\) with a file tool.*$/\1/p')
    cat "$nonce_path" > "$output"
elif [ "${SCIO_STUB_EMPTY_FINAL:-0}" -ne 1 ]; then
    printf 'stub final for: %s\n' "$prompt_text" > "$output"
else
    : > "$output"
fi
printf 'stub stream\n'
EOF
chmod +x "$TMPDIR_TEST/bin/codex"

STUB_PATH="$TMPDIR_TEST/bin:/usr/bin:/bin"

probe_out=$(PATH="$STUB_PATH" "$PROBE" \
    --workdir "$TMPDIR_TEST/work" --model stub-model --sandbox read-only)
printf '%s\n' "$probe_out" | grep -q '^HAS_SEARCH_FLAG=0$' || {
    echo "FAIL [$_TEST_NAME] probe reported a search flag the stub does not expose" >&2
    exit 1
}
printf '%s\n' "$probe_out" | grep -q '^WEB_MODE=config_enable$' || {
    echo "FAIL [$_TEST_NAME] probe did not select config_enable for config + enable" >&2
    exit 1
}
PATH="$STUB_PATH" "$PROBE" --workdir "$TMPDIR_TEST/work" --model stub-model \
    --sandbox read-only --file-read-check > "$TMPDIR_TEST/file-read.out"

set +e
PATH="$TMPDIR_TEST/no-bin" /bin/bash "$PROBE" --workdir "$TMPDIR_TEST/work" \
    > "$TMPDIR_TEST/absent.out" 2>&1
rc=$?
set -e
assert_eq "$rc" "10" "probe must exit 10 when codex is absent"

set +e
SCIO_STUB_NO_WEB=1 PATH="$STUB_PATH" "$PROBE" --workdir "$TMPDIR_TEST/work" --web \
    > "$TMPDIR_TEST/no-web.out" 2>&1
rc=$?
set -e
assert_eq "$rc" "13" "probe must exit 13 when requested web support is absent"
assert_grep '^WEB_MODE=unsupported$' "$TMPDIR_TEST/no-web.out" \
    "unsupported web mode was not reported"

set +e
PATH="$STUB_PATH" "$LAUNCH" --unit missing --workdir "$TMPDIR_TEST/work" \
    --prompt "$TMPDIR_TEST/missing.md" > "$TMPDIR_TEST/missing.out" 2>&1
rc=$?
set -e
[ "$rc" -ne 0 ] || { echo "FAIL [$_TEST_NAME] launch accepted a missing prompt" >&2; exit 1; }

: > "$TMPDIR_TEST/empty.md"
set +e
PATH="$STUB_PATH" "$LAUNCH" --unit empty --workdir "$TMPDIR_TEST/work" \
    --prompt "$TMPDIR_TEST/empty.md" > "$TMPDIR_TEST/empty.out" 2>&1
rc=$?
set -e
[ "$rc" -ne 0 ] || { echo "FAIL [$_TEST_NAME] launch accepted an empty prompt" >&2; exit 1; }

printf 'Implement the bounded fixture task.\n' > "$TMPDIR_TEST/prompt.md"
: > "$TMPDIR_TEST/block"

set +e
held_out=$(SCIO_STUB_BLOCK="$TMPDIR_TEST/block" TMPDIR="$TMPDIR_TEST/runtime" PATH="$STUB_PATH" \
    timeout 2 "$LAUNCH" --unit held --workdir "$TMPDIR_TEST/work" \
    --prompt "$TMPDIR_TEST/prompt.md" --bg 2>&1)
rc=$?
set -e
assert_eq "$rc" "0" "--bg must return before the held stub finishes"
held_pid_file=$(printf '%s\n' "$held_out" | sed -n 's/^PID_FILE=//p')
[ -s "$held_pid_file" ] || { echo "FAIL [$_TEST_NAME] --bg did not write its pid file" >&2; exit 1; }

set +e
SCIO_STUB_BLOCK="$TMPDIR_TEST/block" TMPDIR="$TMPDIR_TEST/runtime" PATH="$STUB_PATH" \
    "$LAUNCH" --unit held --workdir "$TMPDIR_TEST/work" --prompt "$TMPDIR_TEST/prompt.md" --bg \
    > "$TMPDIR_TEST/duplicate.out" 2>&1
rc=$?
set -e
[ "$rc" -ne 0 ] || { echo "FAIL [$_TEST_NAME] a second run acquired the same unit lock" >&2; exit 1; }

set +e
parallel_out=$(SCIO_STUB_BLOCK="$TMPDIR_TEST/block" TMPDIR="$TMPDIR_TEST/runtime" PATH="$STUB_PATH" \
    "$LAUNCH" --unit disjoint --workdir "$TMPDIR_TEST/work" \
    --prompt "$TMPDIR_TEST/prompt.md" --parallel-ok --bg 2>&1)
rc=$?
set -e
assert_eq "$rc" "0" "a different unit with --parallel-ok must be accepted"
parallel_pid_file=$(printf '%s\n' "$parallel_out" | sed -n 's/^PID_FILE=//p')
[ -s "$parallel_pid_file" ] || {
    echo "FAIL [$_TEST_NAME] parallel background run did not write its pid file" >&2
    exit 1
}

rm "$TMPDIR_TEST/block"
for status_path in "${held_pid_file%/pid}/status.tsv" "${parallel_pid_file%/pid}/status.tsv"; do
    attempts=0
    while [ ! -s "$status_path" ] && [ "$attempts" -lt 200 ]; do
        sleep 0.02
        attempts=$((attempts + 1))
    done
    [ -s "$status_path" ] || {
        echo "FAIL [$_TEST_NAME] background child did not record status: $status_path" >&2
        exit 1
    }
done

set +e
TMPDIR="$TMPDIR_TEST/runtime" PATH="$STUB_PATH" "$LAUNCH" --unit expected \
    --workdir "$TMPDIR_TEST/work" --prompt "$TMPDIR_TEST/prompt.md" \
    --expect absent-result.txt > "$TMPDIR_TEST/expect.out" 2>&1
rc=$?
set -e
[ "$rc" -ne 0 ] || { echo "FAIL [$_TEST_NAME] launch accepted an absent --expect path" >&2; exit 1; }

pass
