#!/usr/bin/env bash
# tests/test_status_injected_agent_classified_correctly.sh
# Regression: prior to the kind-partition fix in status.sh, every injected
# entry (regardless of kind) was rendered under "Skills (N effective)" and
# counted toward the skill total — agents and commands disappeared from
# their own sections. Inject one agent and one command into a stack whose
# base already has a non-zero count per kind, then assert each lands under
# its matching section with the count incremented.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

mkdir project && cd project
"$SCIAGENT" activate base >/dev/null

# Baseline `base` role: skills s_a, s_b; agents ag_a; commands c_a.
# Inject ag_b (agent, overlay-only) and c_b (command, overlay-only).
"$SCIAGENT" inject ag_b >/dev/null
"$SCIAGENT" inject c_b  >/dev/null

out=$("$SCIAGENT" status 2>&1)

# Skills section count must be unchanged at 2 (no skill was injected).
if ! printf '%s\n' "$out" | grep -Eq '^Skills \(2 effective\):'; then
    echo "FAIL [$_TEST_NAME] Skills count must remain 2 (injected entries were agent + command, not skill)" >&2
    echo "--- status output ---" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# Skills section must NOT mention the injected agent or command.
skills_block=$(printf '%s\n' "$out" | awk '/^Skills \(/{flag=1; next} /^$/{flag=0} flag')
if printf '%s\n' "$skills_block" | grep -qw ag_b; then
    echo "FAIL [$_TEST_NAME] injected agent ag_b leaked into the Skills section" >&2
    echo "--- skills block ---" >&2
    printf '%s\n' "$skills_block" >&2
    exit 1
fi
if printf '%s\n' "$skills_block" | grep -qw c_b; then
    echo "FAIL [$_TEST_NAME] injected command c_b leaked into the Skills section" >&2
    echo "--- skills block ---" >&2
    printf '%s\n' "$skills_block" >&2
    exit 1
fi

# Sub-agents count: base contributes ag_a (1) + injected ag_b (1) = 2.
if ! printf '%s\n' "$out" | grep -Eq '^Sub-agents \(2 effective'; then
    echo "FAIL [$_TEST_NAME] Sub-agents count must include the injected agent (expected 2)" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi
agents_block=$(printf '%s\n' "$out" | awk '/^Sub-agents \(/{flag=1; next} /^$/{flag=0} flag')
if ! printf '%s\n' "$agents_block" | grep -qw ag_b; then
    echo "FAIL [$_TEST_NAME] injected agent ag_b missing from the Sub-agents section" >&2
    printf '%s\n' "$agents_block" >&2
    exit 1
fi
if ! printf '%s\n' "$agents_block" | grep -Eq 'ag_b[[:space:]]+injected'; then
    echo "FAIL [$_TEST_NAME] injected agent ag_b must be annotated 'injected'" >&2
    printf '%s\n' "$agents_block" >&2
    exit 1
fi

# Slash commands count: base contributes c_a (1) + injected c_b (1) = 2.
if ! printf '%s\n' "$out" | grep -Eq '^Slash commands \(2 effective'; then
    echo "FAIL [$_TEST_NAME] Slash commands count must include the injected command (expected 2)" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi
commands_block=$(printf '%s\n' "$out" | awk '/^Slash commands \(/{flag=1; next} /^$/{flag=0} flag')
if ! printf '%s\n' "$commands_block" | grep -qw c_b; then
    echo "FAIL [$_TEST_NAME] injected command c_b missing from the Slash commands section" >&2
    printf '%s\n' "$commands_block" >&2
    exit 1
fi
if ! printf '%s\n' "$commands_block" | grep -Eq '/c_b[[:space:]]+injected'; then
    echo "FAIL [$_TEST_NAME] injected command c_b must be annotated 'injected'" >&2
    printf '%s\n' "$commands_block" >&2
    exit 1
fi

pass
