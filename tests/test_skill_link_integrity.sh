#!/usr/bin/env bash
# tests/test_skill_link_integrity.sh
#
# Relative markdown links inside skills/ must resolve to a file that exists.
#
# The recurring mistake is writing a cross-skill reference as
# `[label](other-skill-name.md)`, which reads naturally and is always wrong:
# skills are DIRECTORIES, so the body lives at ../<name>/SKILL.md. There were
# 32 of these across 7 skills. Every one named a real skill, so the reference
# was semantically correct and only the path was dead — which is exactly why
# they survived: nothing renders these bodies, and an agent following the link
# just gets nothing back rather than an error.
#
# Cheap to check, so check it: any relative *.md link under skills/ must point
# at something on disk. External (http[s]) and absolute links are ignored —
# their liveness is not this repo's business.
set -u
. "$(dirname "$0")/_lib.sh"

cd "$TOOLKIT_ROOT"

report=$(python3 - <<'PY'
import os, re, sys

skills = {d for d in os.listdir('skills') if os.path.isdir(os.path.join('skills', d))}
# ](target.md) or ](target.md#anchor); skip http(s) and absolute paths.
pat = re.compile(r'\]\((?!https?:|/)([^)\s#]+\.md)(#[^)]*)?\)')

bad = []
for dirpath, dirs, files in os.walk('skills'):
    dirs[:] = [d for d in dirs if d not in ('.git', '.venv', 'node_modules', '__pycache__')]
    for fn in files:
        if not fn.endswith('.md'):
            continue
        p = os.path.join(dirpath, fn)
        try:
            text = open(p, encoding='utf-8').read()
        except (OSError, UnicodeDecodeError):
            continue
        for m in pat.finditer(text):
            tgt = m.group(1)
            if os.path.exists(os.path.normpath(os.path.join(dirpath, tgt))):
                continue
            name = os.path.basename(tgt)[:-3]
            hint = (f"  -> did you mean ../{name}/SKILL.md ?  "
                    f"(skills are directories)") if name in skills else ""
            bad.append(f"{p}: broken link ({tgt}){hint}")

if bad:
    print("\n".join(bad))
    sys.exit(1)
PY
) || {
    echo "FAIL [$_TEST_NAME] broken relative links under skills/:" >&2
    printf '%s\n' "$report" >&2
    exit 1
}

pass
