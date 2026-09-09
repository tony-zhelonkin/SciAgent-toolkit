# lib/scio/common.sh — constants shared across verbs.
#
# Sourced first by every verb in bin/scio, which is where module order is
# decided. A module that needs one of these does not source this file itself:
# the entry point owns the ordering.

# shellcheck shell=bash

# The toolkit's directory name inside a consumer project. `scio` is the name
# ADR-D9 moves to; `SciAgent-toolkit` is where every vendored copy still is
# until the fleet pass renames them. Both are live during that migration, so
# discovery accepts either and prefers the new one — the order here IS the
# preference, and it is deterministic rather than filesystem-dependent.
_SCIO_TOOLKIT_DIRS=(scio SciAgent-toolkit)
