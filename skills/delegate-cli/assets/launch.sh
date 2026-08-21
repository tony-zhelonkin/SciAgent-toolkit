#!/usr/bin/env bash

set -eu

usage() {
    echo "usage: launch.sh --unit NAME --workdir DIR --prompt FILE [options]" >&2
    echo "options: --rules FILE --model MODEL --effort LEVEL --sandbox MODE --bypass --web" >&2
    echo "         --expect PATH [--expect PATH ...] --parallel-ok --bg" >&2
}

fail() {
    code="$1"
    shift
    printf 'launch.sh: %s\n' "$*" >&2
    exit "$code"
}

shell_value() {
    printf '%q' "$1"
}

unit=
workdir=
prompt=
rules=
model=
effort=
sandbox=
bypass=0
want_web=0
parallel_ok=0
background=0
expects=()

while [ "$#" -gt 0 ]; do
    case "$1" in
        --unit)
            [ "$#" -ge 2 ] || { usage; exit 2; }
            unit="$2"
            shift 2
            ;;
        --workdir)
            [ "$#" -ge 2 ] || { usage; exit 2; }
            workdir="$2"
            shift 2
            ;;
        --prompt)
            [ "$#" -ge 2 ] || { usage; exit 2; }
            prompt="$2"
            shift 2
            ;;
        --rules)
            [ "$#" -ge 2 ] || { usage; exit 2; }
            rules="$2"
            shift 2
            ;;
        --model)
            [ "$#" -ge 2 ] || { usage; exit 2; }
            model="$2"
            shift 2
            ;;
        --effort)
            [ "$#" -ge 2 ] || { usage; exit 2; }
            effort="$2"
            shift 2
            ;;
        --sandbox)
            [ "$#" -ge 2 ] || { usage; exit 2; }
            sandbox="$2"
            shift 2
            ;;
        --bypass)
            bypass=1
            shift
            ;;
        --web)
            want_web=1
            shift
            ;;
        --expect)
            [ "$#" -ge 2 ] || { usage; exit 2; }
            expects+=("$2")
            shift 2
            ;;
        --parallel-ok)
            parallel_ok=1
            shift
            ;;
        --bg)
            background=1
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

[ -n "$unit" ] || { usage; exit 2; }
case "$unit" in
    .|..|*[![:alnum:]_.-]*) fail 2 "unit must contain only letters, numbers, dot, underscore, or hyphen" ;;
esac
[ -n "$workdir" ] || { usage; exit 2; }
[ -d "$workdir" ] || fail 2 "workdir is not a directory: $workdir"
[ -n "$prompt" ] || { usage; exit 2; }
[ -f "$prompt" ] || fail 2 "prompt file is missing: $prompt"
[ -s "$prompt" ] || fail 2 "prompt file is empty: $prompt"
if [ -n "$rules" ] && [ ! -f "$rules" ]; then
    fail 2 "rules file is missing: $rules"
fi
if [ "$bypass" -eq 1 ] && [ -n "$sandbox" ]; then
    fail 2 "--bypass and --sandbox are mutually exclusive"
fi

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
if [ "$want_web" -eq 1 ]; then
    probe_args=(--workdir "$workdir" --web)
    [ -z "$model" ] || probe_args+=(--model "$model")
    [ -z "$sandbox" ] || probe_args+=(--sandbox "$sandbox")
    set +e
    probe_output=$("$script_dir/probe.sh" "${probe_args[@]}" 2>&1)
    probe_status=$?
    set -e
    if [ "$probe_status" -ne 0 ]; then
        printf '%s\n' "$probe_output" >&2
        exit "$probe_status"
    fi
    printf '%s\n' "$probe_output" | grep -q '^WEB_MODE=config_enable$' \
        || fail 13 "probe did not authorize the configured web-search mechanism"
fi

run_root="${TMPDIR:-/tmp}/scio-delegate"
lock_root="$run_root/locks"
mkdir -p "$lock_root"
lock="$lock_root/$unit.lock"

lock_pid_live() {
    candidate_lock="$1"
    candidate_pid=
    if [ -r "$candidate_lock/pid" ]; then
        IFS= read -r candidate_pid < "$candidate_lock/pid" || true
    fi
    case "$candidate_pid" in
        ''|*[!0-9]*) return 1 ;;
    esac
    kill -0 "$candidate_pid" 2>/dev/null
}

clear_dead_lock() {
    dead_lock="$1"
    rm -f "$dead_lock/pid" 2>/dev/null || true
    rmdir "$dead_lock" 2>/dev/null || true
}

if ! mkdir "$lock" 2>/dev/null; then
    if lock_pid_live "$lock"; then
        fail 20 "unit is already active: $unit"
    fi
    clear_dead_lock "$lock"
    mkdir "$lock" 2>/dev/null || fail 20 "unit lock cannot be acquired: $unit"
fi
printf '%s\n' "$$" > "$lock/pid"

release_lock() {
    rm -f "$lock/pid" 2>/dev/null || true
    rmdir "$lock" 2>/dev/null || true
}
trap 'release_lock' EXIT HUP INT TERM

for candidate_lock in "$lock_root"/*.lock; do
    [ -d "$candidate_lock" ] || continue
    [ "$candidate_lock" != "$lock" ] || continue
    if lock_pid_live "$candidate_lock"; then
        if [ "$parallel_ok" -ne 1 ]; then
            fail 21 "another delegation unit is active; use --parallel-ok only for disjoint work"
        fi
    else
        clear_dead_lock "$candidate_lock"
        if [ -d "$candidate_lock" ] && [ "$parallel_ok" -ne 1 ]; then
            fail 21 "another delegation lock has no verifiable live pid: $candidate_lock"
        fi
    fi
done

run_id="$(date -u +%Y%m%dT%H%M%SZ)-$$"
run_dir="$run_root/$unit/$run_id"
mkdir -p "$run_dir"
rendered="$run_dir/prompt.md"
final="$run_dir/final.md"
stream="$run_dir/stream.log"
status_file="$run_dir/status.tsv"
pid_file="$run_dir/pid"

if [ -n "$rules" ]; then
    cat "$rules" "$prompt" > "$rendered"
else
    cat "$prompt" > "$rendered"
fi

expect_paths=()
for expected in "${expects[@]}"; do
    case "$expected" in
        /*) expect_paths+=("$expected") ;;
        *) expect_paths+=("$workdir/$expected") ;;
    esac
done

print_paths() {
    printf 'RUN_DIR=%s\n' "$(shell_value "$run_dir")"
    printf 'PROMPT=%s\n' "$(shell_value "$rendered")"
    printf 'FINAL=%s\n' "$(shell_value "$final")"
    printf 'STREAM=%s\n' "$(shell_value "$stream")"
    printf 'STATUS=%s\n' "$(shell_value "$status_file")"
    printf 'PID_FILE=%s\n' "$(shell_value "$pid_file")"
}

path_bytes() {
    measured_path="$1"
    if [ -f "$measured_path" ]; then
        wc -c < "$measured_path" | tr -d '[:space:]'
    else
        printf '0'
    fi
}

record_status() {
    codex_status="$1"
    final_status="$2"
    {
        printf 'kind\tpath\tvalue\n'
        printf 'codex_status\t-\t%s\n' "$codex_status"
        printf 'status\t-\t%s\n' "$final_status"
        printf 'bytes\t%s\t%s\n' "$rendered" "$(path_bytes "$rendered")"
        printf 'bytes\t%s\t%s\n' "$final" "$(path_bytes "$final")"
        printf 'bytes\t%s\t%s\n' "$stream" "$(path_bytes "$stream")"
        for measured_path in "${expect_paths[@]}"; do
            printf 'bytes\t%s\t%s\n' "$measured_path" "$(path_bytes "$measured_path")"
        done
    } > "$status_file"
}

print_bytes() {
    printf 'BYTES\t%s\t%s\n' "$(path_bytes "$rendered")" "$rendered"
    printf 'BYTES\t%s\t%s\n' "$(path_bytes "$final")" "$final"
    printf 'BYTES\t%s\t%s\n' "$(path_bytes "$stream")" "$stream"
    for measured_path in "${expect_paths[@]}"; do
        printf 'BYTES\t%s\t%s\n' "$(path_bytes "$measured_path")" "$measured_path"
    done
}

run_codex() {
    emit_counts="$1"
    codex_args=(exec -C "$workdir" --skip-git-repo-check -o "$final")
    [ -z "$model" ] || codex_args+=(-m "$model")
    [ -z "$effort" ] || codex_args+=(-c "model_reasoning_effort=$effort")
    [ -z "$sandbox" ] || codex_args+=(-s "$sandbox")
    [ "$bypass" -ne 1 ] || codex_args+=(--dangerously-bypass-approvals-and-sandbox)
    if [ "$want_web" -eq 1 ]; then
        codex_args+=(-c tools.web_search=true --enable web_search_request)
    fi

    if codex "${codex_args[@]}" - < "$rendered" > "$stream" 2>&1; then
        codex_status=0
    else
        codex_status=$?
    fi

    final_status="$codex_status"
    if [ ! -s "$final" ] && [ "$final_status" -eq 0 ]; then
        final_status=30
    fi
    for measured_path in "${expect_paths[@]}"; do
        if [ ! -e "$measured_path" ] && [ "$final_status" -eq 0 ]; then
            final_status=31
        fi
    done

    record_status "$codex_status" "$final_status"
    [ "$emit_counts" -ne 1 ] || print_bytes
    return "$final_status"
}

print_paths
printf 'PROMPT_BYTES=%s\n' "$(path_bytes "$rendered")"

if [ "$background" -eq 1 ]; then
    (
        trap 'release_lock' EXIT HUP INT TERM
        printf '%s\n' "$BASHPID" > "$pid_file"
        printf '%s\n' "$BASHPID" > "$lock/pid"
        if run_codex 0; then
            exit 0
        else
            exit $?
        fi
    ) >/dev/null 2>&1 &
    child_pid=$!
    trap - EXIT HUP INT TERM

    attempts=0
    while [ ! -s "$pid_file" ] && [ "$attempts" -lt 100 ]; do
        sleep 0.01
        attempts=$((attempts + 1))
    done
    if [ ! -s "$pid_file" ]; then
        kill "$child_pid" 2>/dev/null || true
        wait "$child_pid" 2>/dev/null || true
        fail 32 "background child did not write its pid"
    fi
    printf 'PID=%s\n' "$(sed -n '1p' "$pid_file")"
    exit 0
fi

printf '%s\n' "$$" > "$pid_file"
if run_codex 1; then
    run_status=0
else
    run_status=$?
fi
exit "$run_status"
