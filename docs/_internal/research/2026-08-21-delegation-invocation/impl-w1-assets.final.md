Implemented Worker 1 only:

- [probe.sh](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/skills/delegate-cli/assets/probe.sh)
- [launch.sh](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/skills/delegate-cli/assets/launch.sh)
- [test_delegate_cli_assets.sh](/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/tests/test_delegate_cli_assets.sh)

Both assets are executable. Tests use a stub `codex`; the real CLI was never invoked. No commits or git write commands were run.

Gate outputs:

```text
$ bash tests/run-all.sh
==== test summary ====
passed: 63
failed: 0
```

```text
$ bin/scio lint --check toolkit --strict
# no output; exit 0
```

```text
$ bash -n skills/delegate-cli/assets/*.sh
# no output; exit 0
```

```text
$ grep -n '\-\-search\|run_in_background' skills/delegate-cli/SKILL.md
# no output; exit 1 (expected: no matches)
```

The baseline agreed with the plan: 62 runners before this test, 63 afterward.

A concurrent modification to `skills/delegate-cli/SKILL.md` appeared after my initial clean-worktree check. I did not touch or revert it. Final status:

```text
 M skills/delegate-cli/SKILL.md
?? skills/delegate-cli/assets/
?? tests/test_delegate_cli_assets.sh
```