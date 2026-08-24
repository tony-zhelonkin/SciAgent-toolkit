#!/usr/bin/env bash
# tests/test_delegate_environments.sh — the launcher's behaviour is predictable
# across the environments it actually meets.
#
# codex does not behave the same everywhere. In scbio-docker dev containers as
# shipped, a SANDBOXED codex cannot read files at all: the container policy
# blocks the namespace or filter operation the sandbox is built on, so the run
# fails in a way that reads like a model refusing to cooperate. It launches
# cleanly only with an authorized bypass. These tests pin what the launcher and
# probe do in each of those conditions, so the answer is the same next time
# rather than rediscovered.
#
# codex is stubbed throughout: no test spends real tokens or depends on a
# network. Set SCIO_DELEGATE_LIVE=1 to add one real-binary smoke test at the end.
#
# Tests:
#   1. codex absent from PATH             -> probe exits 10, names it
#   2. sandbox blocks file reads          -> probe exits 14, SANDBOX_FILE_READ=blocked
#   3. sandbox works                      -> probe reports SANDBOX_FILE_READ=ok
#   4. codex fails under launch           -> state=failed, the REAL exit code recorded
#   5. codex writes no final              -> result_exit=30, codex_exit=0
#   6. --bypass reaches codex             -> the flag appears in the invocation
#   7. --sandbox MODE reaches codex       -> -s MODE appears in the invocation
#   8. the run root follows TMPDIR        -> nothing is written outside it
#   9. codex killed mid-run               -> terminal state, not a frozen 'running'
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
ASSETS="$TOOLKIT_ROOT/skills/delegate-cli/assets"

# _mkenv <name> — isolated runtime root, workdir, prompt. Sets TD, W, P, BIN.
_mkenv() {
    TD="$TMPDIR_TEST/$1"
    BIN="$TD/bin"
    W="$TD/work"
    P="$TD/prompt.md"
    mkdir -p "$TD/tmp" "$W" "$BIN"
    printf 'do a bounded thing\n' > "$P"
}

# Every stub records its own argv so a test can assert what the launcher passed.
_argv_recorder() {
    printf 'printf "%%s\\n" "$@" > "%s/argv"\n' "$TD"
}

# --- 1. codex absent -------------------------------------------------------
# The system directories stay on PATH so the interpreter still resolves; only
# the directory holding codex is dropped. Skipped if codex is installed into a
# system directory on this host, where the condition cannot be constructed.
_mkenv e1
if [ -x /usr/bin/codex ] || [ -x /bin/codex ]; then
    echo "SKIP [$_TEST_NAME] case1: codex is in a system directory here" >&2
else
    set +e
    out1=$(PATH="$BIN:/usr/bin:/bin" "$ASSETS/probe.sh" --workdir "$W" 2>&1); rc1=$?
    set -e
    [ "$rc1" -eq 10 ] || {
        echo "FAIL [$_TEST_NAME] case1: rc was $rc1, expected 10 for absent codex" >&2
        printf '%s\n' "$out1" >&2; exit 1
    }
    printf '%s\n' "$out1" | grep -q 'not available on PATH' || {
        echo "FAIL [$_TEST_NAME] case1: message does not name the cause" >&2
        printf '%s\n' "$out1" >&2; exit 1
    }
fi

# A codex stub whose exec help advertises the flags probe.sh requires, so the
# capability gate passes and the test isolates the condition under study.
_capable_stub() {
    local body="$1"
    cat > "$BIN/codex" <<EOF
#!/usr/bin/env bash
case "\${1:-}" in
  --version) echo "codex-cli-exec 9.9.9"; exit 0 ;;
esac
if [ "\${1:-}" = exec ] && [ "\${2:-}" = --help ]; then
  cat <<'HELP'
Usage: codex exec [OPTIONS] [PROMPT]
  -C, --cd <DIR>
      --skip-git-repo-check
  -o, --output-last-message <FILE>
  -m, --model <MODEL>
  -s, --sandbox <MODE>
  -c, --config <KEY=VALUE>
      --enable <FEATURE>
      --dangerously-bypass-approvals-and-sandbox
  Reads the prompt from stdin when PROMPT is '-'.
HELP
  exit 0
fi
$body
EOF
    chmod +x "$BIN/codex"
}

# --- 2. the sandbox blocks file reads -------------------------------------
# This is the scbio-docker condition. The bwrap denial goes to the stream, which
# is where probe.sh looks, and the exit is nonzero.
_mkenv e2
_capable_stub 'cat > /dev/null
echo "bwrap: setting up uid map: Operation not permitted" >&2
exit 1'
set +e
out2=$(TMPDIR="$TD/tmp" PATH="$BIN:$PATH" \
    "$ASSETS/probe.sh" --workdir "$W" --sandbox workspace-write --file-read-check 2>&1); rc2=$?
set -e
[ "$rc2" -eq 14 ] || {
    echo "FAIL [$_TEST_NAME] case2: rc was $rc2, expected 14" >&2
    printf '%s\n' "$out2" >&2; exit 1
}
printf '%s\n' "$out2" | grep -q 'SANDBOX_FILE_READ=blocked' || {
    echo "FAIL [$_TEST_NAME] case2: the blocked condition is not reported as a fact" >&2
    printf '%s\n' "$out2" >&2; exit 1
}
# The remedy must be stated, because rewording the prompt is not one.
printf '%s\n' "$out2" | grep -q 'bypass' || {
    echo "FAIL [$_TEST_NAME] case2: no remedy named" >&2
    printf '%s\n' "$out2" >&2; exit 1
}

# --- 3. the sandbox works -------------------------------------------------
_mkenv e3
_capable_stub 'out=; prev=; pfile=
for a in "$@"; do [ "$prev" = "-o" ] && out="$a"; prev="$a"; done
prompt=$(cat)
# Echo back the nonce the probe asked us to read.
nfile=$(printf "%s" "$prompt" | sed -n "s#.*Read \\(/[^ ]*nonce.txt\\).*#\\1#p" | head -n 1)
cat "$nfile" > "$out"
exit 0'
set +e
out3=$(TMPDIR="$TD/tmp" PATH="$BIN:$PATH" \
    "$ASSETS/probe.sh" --workdir "$W" --file-read-check 2>&1); rc3=$?
set -e
[ "$rc3" -eq 0 ] || {
    echo "FAIL [$_TEST_NAME] case3: rc was $rc3 for a working sandbox" >&2
    printf '%s\n' "$out3" >&2; exit 1
}
printf '%s\n' "$out3" | grep -q 'SANDBOX_FILE_READ=ok' || {
    echo "FAIL [$_TEST_NAME] case3: a working file read is not reported" >&2
    printf '%s\n' "$out3" >&2; exit 1
}

# _state <runtime-root> <unit>
_state() {
    awk -F'\t' '$1 == "state" { print $2 }' \
        "$1/scio-delegate/$2/current/status.tsv" 2>/dev/null | tail -n 1
}
_field() {
    awk -F'\t' -v k="$3" '$1 == k { print $2 }' \
        "$1/scio-delegate/$2/current/status.tsv" 2>/dev/null | tail -n 1
}

# --- 4. codex fails: the real exit code survives -------------------------
# The old wrapper's printed exit code was decorative. The recorded one must be
# the child's actual status, because that is what a caller acts on.
_mkenv e4
_capable_stub 'cat > /dev/null; echo "sandbox denied" >&2; exit 7'
set +e
TMPDIR="$TD/tmp" PATH="$BIN:$PATH" "$ASSETS/launch.sh" \
    --unit u_fail --workdir "$W" --prompt "$P" --sandbox workspace-write \
    > "$TD/out" 2>&1
rc4=$?
set -e
[ "$rc4" -eq 7 ] || {
    echo "FAIL [$_TEST_NAME] case4: launcher rc was $rc4, expected codex's 7" >&2
    cat "$TD/out" >&2; exit 1
}
[ "$(_state "$TD/tmp" u_fail)" = "failed" ] || {
    echo "FAIL [$_TEST_NAME] case4: state was '$(_state "$TD/tmp" u_fail)'" >&2; exit 1
}
[ "$(_field "$TD/tmp" u_fail codex_exit)" = "7" ] || {
    echo "FAIL [$_TEST_NAME] case4: codex_exit was '$(_field "$TD/tmp" u_fail codex_exit)'" >&2
    exit 1
}

# --- 5. codex succeeds but writes nothing --------------------------------
_mkenv e5
_capable_stub 'cat > /dev/null; echo "done talking"; exit 0'
set +e
TMPDIR="$TD/tmp" PATH="$BIN:$PATH" "$ASSETS/launch.sh" \
    --unit u_nofinal --workdir "$W" --prompt "$P" > "$TD/out" 2>&1
rc5=$?
set -e
[ "$rc5" -eq 30 ] || {
    echo "FAIL [$_TEST_NAME] case5: rc was $rc5, expected 30 for an absent final" >&2
    cat "$TD/out" >&2; exit 1
}
[ "$(_field "$TD/tmp" u_nofinal codex_exit)" = "0" ] || {
    echo "FAIL [$_TEST_NAME] case5: codex succeeded, so codex_exit must be 0" >&2; exit 1
}

# --- 6+7. the authorization flags actually reach codex -------------------
# The container that blocks the sandbox is the reason --bypass exists, so the
# flag reaching the child is the thing worth pinning.
_mkenv e6
_capable_stub "$(_argv_recorder)
out=; prev=
for a in \"\$@\"; do [ \"\$prev\" = \"-o\" ] && out=\"\$a\"; prev=\"\$a\"; done
cat > /dev/null; printf 'ok\n' > \"\$out\"; exit 0"
TMPDIR="$TD/tmp" PATH="$BIN:$PATH" "$ASSETS/launch.sh" \
    --unit u_bypass --workdir "$W" --prompt "$P" --bypass > "$TD/out" 2>&1
grep -qx -- '--dangerously-bypass-approvals-and-sandbox' "$TD/argv" || {
    echo "FAIL [$_TEST_NAME] case6: --bypass did not reach codex" >&2
    cat "$TD/argv" >&2; exit 1
}

_mkenv e7
_capable_stub "$(_argv_recorder)
out=; prev=
for a in \"\$@\"; do [ \"\$prev\" = \"-o\" ] && out=\"\$a\"; prev=\"\$a\"; done
cat > /dev/null; printf 'ok\n' > \"\$out\"; exit 0"
TMPDIR="$TD/tmp" PATH="$BIN:$PATH" "$ASSETS/launch.sh" \
    --unit u_sbx --workdir "$W" --prompt "$P" --sandbox read-only > "$TD/out" 2>&1
grep -qx -- 'read-only' "$TD/argv" || {
    echo "FAIL [$_TEST_NAME] case7: --sandbox value did not reach codex" >&2
    cat "$TD/argv" >&2; exit 1
}
# --bypass and --sandbox together are refused rather than silently reconciled.
set +e
out7=$(TMPDIR="$TD/tmp" PATH="$BIN:$PATH" "$ASSETS/launch.sh" \
    --unit u_both --workdir "$W" --prompt "$P" --bypass --sandbox read-only 2>&1); rc7=$?
set -e
[ "$rc7" -ne 0 ] || {
    echo "FAIL [$_TEST_NAME] case7: --bypass with --sandbox was accepted" >&2; exit 1
}

# --- 8. the run root follows TMPDIR -------------------------------------
# A container gets its own TMPDIR; writing outside it would put run state
# somewhere the next session cannot find and the host may not allow.
_mkenv e8
_capable_stub 'out=; prev=
for a in "$@"; do [ "$prev" = "-o" ] && out="$a"; prev="$a"; done
cat > /dev/null; printf "ok\n" > "$out"; exit 0'
TMPDIR="$TD/tmp" PATH="$BIN:$PATH" "$ASSETS/launch.sh" \
    --unit u_tmp --workdir "$W" --prompt "$P" > "$TD/out" 2>&1
[ -d "$TD/tmp/scio-delegate/u_tmp" ] || {
    echo "FAIL [$_TEST_NAME] case8: run state is not under TMPDIR" >&2; exit 1
}
run_dir=$(sed -n 's/^RUN_DIR=//p' "$TD/out" | tr -d "'")
case "$run_dir" in
    "$TD/tmp"/*) ;;
    *) echo "FAIL [$_TEST_NAME] case8: RUN_DIR '$run_dir' escapes TMPDIR" >&2; exit 1 ;;
esac

# --- 9. a killed child still reaches a terminal state -------------------
# SIGKILL on the child is the case a status writer can still record, unlike a
# SIGKILL on the launcher. Anything left at 'running' would read as in-flight
# forever.
_mkenv e9
_capable_stub 'cat > /dev/null; kill -9 $$'
set +e
TMPDIR="$TD/tmp" PATH="$BIN:$PATH" "$ASSETS/launch.sh" \
    --unit u_killed --workdir "$W" --prompt "$P" > "$TD/out" 2>&1
rc9=$?
set -e
[ "$rc9" -ne 0 ] || {
    echo "FAIL [$_TEST_NAME] case9: a killed child reported success" >&2; exit 1
}
st9=$(_state "$TD/tmp" u_killed)
case "$st9" in
    failed|interrupted) ;;
    *)
        echo "FAIL [$_TEST_NAME] case9: state was '$st9', expected a terminal state" >&2
        exit 1
        ;;
esac

# --- 10. opt-in live smoke test ----------------------------------------
# Costs tokens and needs a configured codex, so it runs only on request.
if [ "${SCIO_DELEGATE_LIVE:-0}" = "1" ]; then
    _mkenv live
    if ! command -v codex >/dev/null 2>&1; then
        echo "SKIP [$_TEST_NAME] case10: SCIO_DELEGATE_LIVE=1 but codex is absent" >&2
    else
        printf 'Reply with exactly: SCIO_LIVE_OK\n' > "$P"
        set +e
        TMPDIR="$TD/tmp" "$ASSETS/launch.sh" --unit u_live --workdir "$W" \
            --prompt "$P" > "$TD/out" 2>&1
        rcl=$?
        set -e
        echo "INFO [$_TEST_NAME] case10: live launcher rc=$rcl state=$(_state "$TD/tmp" u_live)" >&2
        [ "$rcl" -eq 0 ] || {
            echo "FAIL [$_TEST_NAME] case10: live run failed; probe --file-read-check for the cause" >&2
            cat "$TD/out" >&2
            exit 1
        }
    fi
fi

pass
