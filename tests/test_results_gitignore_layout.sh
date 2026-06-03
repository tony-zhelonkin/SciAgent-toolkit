#!/usr/bin/env bash
# Phase-based 03_results/ .gitignore: data artifacts ignored, .gitkeep tracked,
# objects/ checkpoints hard-ignored.

. "$(dirname "$0")/_lib.sh"

setup_tmpdir

proj="$TMPDIR_TEST/proj"
"$SCIAGENT_TOOLKIT/bin/sciagent" new project "$proj" --type analysis >/dev/null 2>&1 \
    || { echo "FAIL [$_TEST_NAME] scaffold failed" >&2; exit 1; }

# git check-ignore needs a repo with the seeded .gitignore.
git -C "$proj" init -q
git -C "$proj" add .gitignore >/dev/null 2>&1

# NOTE: figure/table gitignore rules were deliberately commented out in
# .gitignore-seed ("Decided to leave for now"), so we no longer assert
# that PDFs/PNGs/CSVs in phase figures/tables dirs are ignored.

# Checkpoint state objects must be hard-ignored.
touch "$proj/03_results/objects/data.h5ad"
if git -C "$proj" check-ignore -q "03_results/objects/data.h5ad"; then
    :
else
    echo "FAIL [$_TEST_NAME] objects/data.h5ad not hard-ignored" >&2
    exit 1
fi

# Scratch contents must be ignored, but the skeleton .gitkeep stays tracked.
touch "$proj/03_results/_scratch/throwaway.png"
if git -C "$proj" check-ignore -q "03_results/_scratch/throwaway.png"; then
    :
else
    echo "FAIL [$_TEST_NAME] _scratch/throwaway.png not gitignored" >&2
    exit 1
fi
if git -C "$proj" check-ignore -q "03_results/_scratch/.gitkeep"; then
    echo "FAIL [$_TEST_NAME] _scratch/.gitkeep was gitignored (should be tracked)" >&2
    exit 1
fi

pass
