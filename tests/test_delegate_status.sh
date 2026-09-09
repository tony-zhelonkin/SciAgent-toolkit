#!/usr/bin/env bash
# tests/test_delegate_status.sh — delegate-cli's two observability contracts.
#
# The forensic trace of the JR-MC session found status.tsv appearing only after
# the child was reaped, so it could not answer the question a caller actually
# has: how far along is this. And --bg detached the launcher with no
# harness-owned completion signal, so every long dead-time launch in that
# session was a run nobody was waiting on. Disk observability and orchestrator
# wake-up are separate contracts and both are needed.
#
# codex is stubbed. The real binary must not run from a test.
#
# Tests:
#   1. status appears while the child runs, state=running, stream_bytes climbing
#   2. the per-unit `current` symlink resolves to the run
#   3. a clean run ends succeeded, codex_exit=0, result_exit=0
#   4. an absent expected artifact is result_exit=31 with codex_exit=0
#   5. --wait blocks to a terminal state and exits with the run's result
#   6. --unit for a unit that never ran exits 3
#   7. a status file with no schema row exits 4
#   8. --bg is refused with an actionable message
#   9. status.sh invents no phase or percentage
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
ASSETS="$TOOLKIT_ROOT/skills/delegate-cli/assets"

# A codex stub: streams for <secs>, then writes whatever -o names.
# _stub <dir> <secs> [--no-final]
_stub() {
    local bin="$1/bin" secs="$2" nofinal="${3:-}"
    mkdir -p "$bin"
    cat > "$bin/codex" <<EOF
#!/usr/bin/env bash
out=; prev=
for a in "\$@"; do [ "\$prev" = "-o" ] && out="\$a"; prev="\$a"; done
cat > /dev/null
i=0
while [ "\$i" -lt $secs ]; do echo "thinking \$i"; sleep 1; i=\$(( i + 1 )); done
[ "$nofinal" = "--no-final" ] || printf 'the answer\n' > "\$out"
EOF
    chmod +x "$bin/codex"
}

# _env <name> — an isolated runtime root, workdir and prompt. Sets TD/W/P.
_env() {
    TD="$TMPDIR_TEST/$1"
    mkdir -p "$TD/tmp" "$TD/work"
    W="$TD/work"
    P="$TD/prompt.md"
    printf 'do a thing\n' > "$P"
}

# --- 1+2+4. in-flight status, the current symlink, and a missing artifact ----
_env e1
_stub "$TD" 5
(
    export TMPDIR="$TD/tmp" PATH="$TD/bin:$PATH"
    "$ASSETS/launch.sh" --unit u_run --workdir "$W" --prompt "$P" \
        --expect "$W/never.txt" >"$TD/out" 2>"$TD/err" &
    lp=$!
    sleep 3

    mid=$("$ASSETS/status.sh" --unit u_run 2>&1) || {
        echo "FAIL [$_TEST_NAME] case1: status.sh failed mid-run" >&2
        printf '%s\n' "$mid" >&2; exit 1
    }
    printf '%s\n' "$mid" | grep -q 'state=running' || {
        echo "FAIL [$_TEST_NAME] case1: expected state=running mid-flight" >&2
        printf '%s\n' "$mid" >&2; exit 1
    }
    # The point of the rebuild: bytes are observable BEFORE the child is reaped.
    sb=$(printf '%s\n' "$mid" | sed -n 's/.*stream_bytes=\([0-9]*\).*/\1/p')
    [ "${sb:-0}" -gt 0 ] || {
        echo "FAIL [$_TEST_NAME] case1: stream_bytes was $sb mid-run" >&2
        printf '%s\n' "$mid" >&2; exit 1
    }

    [ -L "$TD/tmp/scio-delegate/u_run/current" ] || {
        echo "FAIL [$_TEST_NAME] case2: no current symlink" >&2; exit 1
    }
    [ -f "$TD/tmp/scio-delegate/u_run/current/status.tsv" ] || {
        echo "FAIL [$_TEST_NAME] case2: current/ does not resolve to the run" >&2; exit 1
    }

    set +e; wait "$lp"; lrc=$?; set -e
    [ "$lrc" -eq 31 ] || {
        echo "FAIL [$_TEST_NAME] case4: launcher rc was $lrc, expected 31" >&2; exit 1
    }
    fin=$("$ASSETS/status.sh" --unit u_run)
    printf '%s\n' "$fin" | grep -q 'state=failed' || {
        echo "FAIL [$_TEST_NAME] case4: expected state=failed" >&2
        printf '%s\n' "$fin" >&2; exit 1
    }
    printf '%s\n' "$fin" | grep -q 'codex_exit=0.*result_exit=31' || {
        echo "FAIL [$_TEST_NAME] case4: codex succeeded but the artifact was absent — expected codex_exit=0 result_exit=31" >&2
        printf '%s\n' "$fin" >&2; exit 1
    }
) || exit 1

# --- 3+5. a clean run, and --wait shadowing it ------------------------------
_env e2
_stub "$TD" 3
(
    export TMPDIR="$TD/tmp" PATH="$TD/bin:$PATH"
    "$ASSETS/launch.sh" --unit u_ok --workdir "$W" --prompt "$P" \
        >"$TD/out" 2>"$TD/err" &
    lp=$!
    # Wait for the status file to exist, then shadow the run as a harness task
    # would: this process's exit IS the wake-up signal.
    tries=0
    while [ ! -f "$TD/tmp/scio-delegate/u_ok/current/status.tsv" ] && [ "$tries" -lt 200 ]; do
        sleep 0.05; tries=$(( tries + 1 ))
    done
    set +e
    shadow=$("$ASSETS/status.sh" --unit u_ok --wait 2>&1); src=$?
    wait "$lp"; lrc=$?
    set -e
    [ "$lrc" -eq 0 ] || {
        echo "FAIL [$_TEST_NAME] case3: clean launcher rc was $lrc" >&2
        cat "$TD/err" >&2; exit 1
    }
    printf '%s\n' "$shadow" | grep -q 'state=succeeded' || {
        echo "FAIL [$_TEST_NAME] case3: expected state=succeeded" >&2
        printf '%s\n' "$shadow" >&2; exit 1
    }
    printf '%s\n' "$shadow" | grep -q 'result_exit=0' || {
        echo "FAIL [$_TEST_NAME] case3: expected result_exit=0" >&2
        printf '%s\n' "$shadow" >&2; exit 1
    }
    [ "$src" -eq 0 ] || {
        echo "FAIL [$_TEST_NAME] case5: --wait exited $src for a successful run" >&2; exit 1
    }
) || exit 1

# --- 6. a unit that never ran ----------------------------------------------
_env e3
set +e
out6=$(TMPDIR="$TD/tmp" "$ASSETS/status.sh" --unit u_absent 2>&1); rc6=$?
set -e
[ "$rc6" -eq 3 ] || {
    echo "FAIL [$_TEST_NAME] case6: rc was $rc6, expected 3" >&2
    printf '%s\n' "$out6" >&2; exit 1
}

# --- 7. a broken record is not the same failure as an absent one -----------
_env e4
printf 'state\trunning\n' > "$TD/broken.tsv"
set +e
out7=$("$ASSETS/status.sh" --file "$TD/broken.tsv" 2>&1); rc7=$?
set -e
[ "$rc7" -eq 4 ] || {
    echo "FAIL [$_TEST_NAME] case7: rc was $rc7, expected 4 for a schema-less record" >&2
    printf '%s\n' "$out7" >&2; exit 1
}

# --- 8. --bg is refused, not ignored --------------------------------------
_env e5
set +e
out8=$("$ASSETS/launch.sh" --unit u_bg --workdir "$W" --prompt "$P" --bg 2>&1); rc8=$?
set -e
[ "$rc8" -ne 0 ] || {
    echo "FAIL [$_TEST_NAME] case8: --bg was accepted" >&2; exit 1
}
printf '%s\n' "$out8" | grep -q 'background-task mechanism' || {
    echo "FAIL [$_TEST_NAME] case8: --bg refusal is not actionable" >&2
    printf '%s\n' "$out8" >&2; exit 1
}

# --- 9. no manufactured progress -----------------------------------------
# Nothing on disk supports a phase or a percentage, so the OUTPUT must not carry
# one. Asserted against emitted text, not the source, which discusses the rule.
# Hung versus slow stays the reader's judgement from the ages reported.
_env e6
printf 'schema\t1\nstate\trunning\nunit\tu\nrun_id\tr\n' > "$TD/min.tsv"
out9=$("$ASSETS/status.sh" --file "$TD/min.tsv")
if printf '%s\n' "$out9" | grep -qiE 'percent|progress=|phase=|[0-9]+%|eta='; then
    echo "FAIL [$_TEST_NAME] case9: output manufactures a phase, percentage or ETA" >&2
    printf '%s\n' "$out9" >&2
    exit 1
fi
# The quantities it does report are the ones the trace found useful.
for field in elapsed_s heartbeat_age_s activity_age_s stream_bytes expected_present; do
    printf '%s\n' "$out9" | grep -q "$field=" || {
        echo "FAIL [$_TEST_NAME] case9: output omits $field" >&2
        printf '%s\n' "$out9" >&2
        exit 1
    }
done

pass
