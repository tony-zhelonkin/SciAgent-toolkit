#!/usr/bin/env bash
# Regression: real role YAMLs (e.g. roles/base.yaml) write skill lines as
# `  - <name>    # description`. The stack-mount guard in eject.sh must
# still recognise them; otherwise `sciagent eject <stack-mounted-skill>`
# silently exits 0 and treats the skill as merely "not injected".
#
# Synthesise a fixture role with trailing-comment skill lines, activate it,
# attempt to eject one of the listed skills, and assert that eject refuses
# with the `sciagent deactivate` pointer.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

# Rewrite base.yaml so every skill line carries a trailing comment, mirroring
# the real roles/base.yaml style.
cat > "$FAKE/roles/base.yaml" <<'EOF'
name: base
description: fixture base role with commented skill lines
skills:
  - s_a                                # short description for s_a
  - s_b   # tight spacing variant
agents:
  - ag_a
commands:
  - c_a
EOF

mkdir project && cd project
"$SCIAGENT" activate base >/dev/null

# Both s_a (wide gap) and s_b (tight gap) are stack-mounted with trailing
# comments. Each eject attempt must fail with the deactivate pointer.
for sk in s_a s_b; do
    out=$("$SCIAGENT" eject "$sk" 2>&1)
    rc=$?
    if [[ "$rc" -eq 0 ]]; then
        echo "FAIL [$_TEST_NAME] expected exit 1 ejecting commented stack-mounted '$sk', got 0" >&2
        echo "output: $out" >&2
        exit 1
    fi
    printf '%s\n' "$out" | grep -qi 'deactivate' || {
        echo "FAIL [$_TEST_NAME] eject of '$sk' must mention 'sciagent deactivate', got: $out" >&2
        exit 1
    }
done

pass
