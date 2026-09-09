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

# The lock assertions below need a launcher that holds its unit while another
# tries to take it. --bg used to provide that; it is gone, so the test does what
# a harness does — backgrounds the FOREGROUND launcher and waits for the status
# file it writes before spawning the child.
_hold() {
    SCIO_STUB_BLOCK="$TMPDIR_TEST/block" TMPDIR="$TMPDIR_TEST/runtime" PATH="$STUB_PATH" \
        "$LAUNCH" --unit "$1" --workdir "$TMPDIR_TEST/work" \
        --prompt "$TMPDIR_TEST/prompt.md" ${2:+"$2"} \
        > "$TMPDIR_TEST/$1.out" 2>&1 &
    HOLD_PID=$!
}

# Wait for state=running, not merely for the file: the launcher writes
# `starting` before it spawns codex, so the file exists before a child does.
_await_status() {
    local path="$TMPDIR_TEST/runtime/scio-delegate/$1/current/status.tsv" tries=0 state=
    while [ "$tries" -lt 300 ]; do
        if [ -s "$path" ]; then
            state=$(awk -F'\t' '$1 == "state" { print $2 }' "$path" | tail -n 1)
            [ "$state" = "running" ] && return 0
        fi
        sleep 0.02
        tries=$(( tries + 1 ))
    done
    echo "FAIL [$_TEST_NAME] unit $1 never reached state=running (last: ${state:-none})" >&2
    exit 1
}

_hold held
held_pid=$HOLD_PID
_await_status held

# Status exists while the child is still blocked. Before the rebuild this file
# appeared only after the child was reaped, so an in-flight question had no
# answer at all.
held_state=$(awk -F'\t' '$1 == "state" { print $2 }' \
    "$TMPDIR_TEST/runtime/scio-delegate/held/current/status.tsv" | tail -n 1)
assert_eq "$held_state" "running" "a held run must report state=running"

set +e
SCIO_STUB_BLOCK="$TMPDIR_TEST/block" TMPDIR="$TMPDIR_TEST/runtime" PATH="$STUB_PATH" \
    "$LAUNCH" --unit held --workdir "$TMPDIR_TEST/work" --prompt "$TMPDIR_TEST/prompt.md" \
    > "$TMPDIR_TEST/duplicate.out" 2>&1
rc=$?
set -e
[ "$rc" -ne 0 ] || { echo "FAIL [$_TEST_NAME] a second run acquired the same unit lock" >&2; exit 1; }

_hold disjoint --parallel-ok
parallel_pid=$HOLD_PID
_await_status disjoint

# A different active unit WITHOUT --parallel-ok is refused: independence is an
# affirmation the caller makes, never an inference the launcher draws.
set +e
SCIO_STUB_BLOCK="$TMPDIR_TEST/block" TMPDIR="$TMPDIR_TEST/runtime" PATH="$STUB_PATH" \
    "$LAUNCH" --unit third --workdir "$TMPDIR_TEST/work" --prompt "$TMPDIR_TEST/prompt.md" \
    > "$TMPDIR_TEST/third.out" 2>&1
rc=$?
set -e
[ "$rc" -ne 0 ] || {
    echo "FAIL [$_TEST_NAME] a concurrent unit was accepted without --parallel-ok" >&2
    exit 1
}

rm "$TMPDIR_TEST/block"
set +e
wait "$held_pid"; wait "$parallel_pid"
set -e

for unit_name in held disjoint; do
    status_path="$TMPDIR_TEST/runtime/scio-delegate/$unit_name/current/status.tsv"
    final_state=$(awk -F'\t' '$1 == "state" { print $2 }' "$status_path" | tail -n 1)
    case "$final_state" in
        succeeded|failed) ;;
        *)
            echo "FAIL [$_TEST_NAME] unit $unit_name ended in state '$final_state'" >&2
            exit 1
            ;;
    esac
done

set +e
TMPDIR="$TMPDIR_TEST/runtime" PATH="$STUB_PATH" "$LAUNCH" --unit expected \
    --workdir "$TMPDIR_TEST/work" --prompt "$TMPDIR_TEST/prompt.md" \
    --expect absent-result.txt > "$TMPDIR_TEST/expect.out" 2>&1
rc=$?
set -e
[ "$rc" -ne 0 ] || { echo "FAIL [$_TEST_NAME] launch accepted an absent --expect path" >&2; exit 1; }

pass
