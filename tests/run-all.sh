#!/usr/bin/env bash
# tests/run-all.sh — run every test_*.sh in this dir, summarize.

set -u
cd "$(dirname "$0")"

pass=0
fail=0
failed_tests=()

for t in test_*.sh; do
    [[ -f "$t" ]] || continue
    if bash "$t"; then
        pass=$((pass + 1))
    else
        fail=$((fail + 1))
        failed_tests+=("$t")
    fi
done

echo
echo "==== test summary ===="
echo "passed: $pass"
echo "failed: $fail"
if [[ $fail -gt 0 ]]; then
    echo "failures:"
    for t in "${failed_tests[@]}"; do
        echo "  - $t"
    done
    exit 1
fi
