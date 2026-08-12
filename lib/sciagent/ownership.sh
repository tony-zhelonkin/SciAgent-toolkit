# lib/sciagent/ownership.sh — the ownership discipline for MATERIALIZED bodies.
#
# WHY THIS MODULE EXISTS (moved here 2026-08-11)
#
# `activate` writes a handful of toolkit-owned file bodies into a project with a
# plain `cp` — no placeholder substitution. Each one needs the same three
# questions answered before it may be written, refreshed, or removed:
#
#   1. Is the file absent?                              → write it, record it.
#   2. Are these bytes ones the toolkit ever shipped?    → ours; safe to refresh.
#   3. Anything else?                                    → the user's; cede it.
#
# That machinery was born inside claude_settings.sh, because its first two
# consumer was .claude/hooks/*.sh. It is not a
# Claude concern. Its second consumer set — the 02_analysis/helpers/*.{py,R}
# shims (symlinks.sh) — is analysis-repo LAYOUT, imported by R and Python and
# entirely unrelated to Claude Code. Calling a `claude_settings_*` function to
# write an R file would entrench a misnaming, so the generic half lives here
# and both callers depend on this module rather than on each other:
#
#   claude_settings.sh  →  ownership_ensure_body / ownership_teardown_body
#   symlinks.sh         →  ownership_ensure_body / ownership_teardown_body
#
# DEPENDENCY CLOSURE (bin/sciagent's VERB_MODULES). This module is required by
# every verb whose closure loads claude_settings.sh or symlinks.sh — activate,
# deactivate, status, update — whether or not the code path a given
# invocation takes actually calls into it. The closure is over the MODULES
# LOADED, not over the branches taken: a module that references a function
# nobody sourced is a "command not found" waiting for the first mutation that
# reaches that line (the lesson of 4b291b2). `validate` and `lint` load neither
# caller and are therefore correctly without it.
#
# Content hashing routes through block.sh's sciagent_sha1_file — block.sh is the
# single home for the hashing primitive in lib/ (the 5.3 invariant, enforced by
# tests/test_block_marker_boundary.sh, which greps for the tool name literally,
# so do not name it here either).

# shellcheck shell=bash

# ----- Ownership records ----------------------------------------------------
#
# Every artifact this module writes carries a SNAPSHOT of the content hash at
# the moment we wrote/adopted it, taken by the caller-supplied <state> path and
# consulted exactly once at teardown:
#
#   State file ABSENT  → we made no change here; teardown touches nothing.
#   State file PRESENT, current hash MATCHES → our content, untouched since →
#     teardown reverses it exactly.
#   State file PRESENT, hash MISMATCH → the user edited or replaced it since →
#     teardown leaves it alone and warns on stderr. It never re-checks: the
#     state file is removed either way, so the decision is made once and the
#     artifact is the user's from then on. That is what keeps teardown
#     idempotent — a second run finds no state file and is a silent no-op
#     rather than re-warning forever.
#
# A CEDED marker (the <ceded> path) records that we found a body we did not
# write, said so once, and stepped back. It holds the content hash at the moment
# we ceded, so the warning fires again if the user later replaces the file with
# something ELSE (a genuinely new situation) but not on every activate for the
# same file. Teardown deletes these markers without touching the files they name
# — the whole point is that those files are not ours to remove.

# ----- Content-provenance: "did WE write this file, and is it untouched?" ----
#
# The ownership records above answer that question only for artifacts created
# AFTER the records existed. Measured across the real fleet on 2026-08-11: not
# one consumer has a .sciagent/hook_state directory, because every one of them
# was provisioned before those records were introduced. So the record is absent
# in exactly the population that needs repairing, and "absent -> touch nothing"
# would make any record-based refresh a no-op everywhere it matters.
#
# The bytes are the remaining evidence, and they are conclusive. These bodies
# are materialized with a plain `cp` (no placeholder substitution — enforced by
# tests/test_template_provenance.sh), so a project's copy is byte-identical to
# whichever template version produced it. If a file's content equals SOME
# version of the template this toolkit has ever shipped, then the toolkit wrote
# it and nobody has edited it since. That is precisely the condition under
# which overwriting is safe, and it is decidable offline from a committed
# manifest of historical hashes.
#
# Anything else — content matching no version we ever shipped — is the user's,
# and is ceded permanently (see ownership_ensure_body).
_SCIAGENT_TEMPLATE_PROVENANCE="templates/PROVENANCE.sha1"

# ownership_template_hash_known <hash> <template-rel-path>
# True if <hash> is a content hash <template-rel-path> has ever had.
# Absent/unreadable manifest → false, so a stripped install degrades to the old
# create-if-absent behaviour rather than overwriting on a guess.
ownership_template_hash_known() {
    local hash="$1" rel="$2"
    local mf="$SCIAGENT_TOOLKIT/$_SCIAGENT_TEMPLATE_PROVENANCE"
    [[ -n "$hash" && -f "$mf" ]] || return 1
    local h r
    while read -r h r; do
        [[ "$h" == \#* || -z "$h" ]] && continue
        [[ "$h" == "$hash" && "$r" == "$rel" ]] && return 0
    done < "$mf"
    return 1
}

# ownership_ensure_body <src> <dst> <template-rel> <state-file> <ceded-file> [exec|plain]
#
# The single materialize-and-keep-current routine behind claude_settings.sh's
# ensure_hooks and symlinks.sh's helper_shims_ensure. Four
# outcomes, in order:
#
#   dst ABSENT            → write it, record the hash. (Original behaviour.)
#   dst == template       → nothing to write. Adopt it: record the hash if we
#                           had none, so `deactivate` can reverse a body that
#                           predates the ownership records.
#   dst is a KNOWN older  → REFRESH — this is the defect this exists to fix.
#   version (record hash    Verified against the fleet: all 8 consumers carrying
#   or manifest hash)       no_ephemeral.sh hold 12dafe8a3e39, the af20086
#                           template, and none of them ever saw de5719e's fix.
#   anything else         → the user's file. Warn ONCE, drop the record, write
#                           a ceded marker, and never touch it again. Warning
#                           on every activate is what the marker prevents; the
#                           teardown path already made this decision once and
#                           for all for the same reason.
#
# The 6th argument decides the execute bit, and every call site states it:
#   exec  (default) — chmod +x on every run, regardless of ownership. A hook
#                     that is not executable is silently dead, and
#                     that is never what the user meant by editing its contents.
#   plain           — modes are left exactly as `cp` produced them. Correct for
#                     the R/Python helper shims: they are imported, never run,
#                     and an executable data file is a lie about how it is used.
ownership_ensure_body() {
    local src="$1" dst="$2" rel="$3" state="$4" ceded="$5" mode="${6:-exec}"
    [[ -f "$src" ]] || return 0

    mkdir -p "$(dirname "$dst")"

    if [[ ! -f "$dst" ]]; then
        cp "$src" "$dst" || { echo "sciagent: failed to write $dst" >&2; return 1; }
        [[ "$mode" == exec ]] && chmod +x "$dst"
        echo "wrote: $dst"
        mkdir -p "$(dirname "$state")"
        sciagent_sha1_file "$dst" > "$state"
        rm -f "$ceded"
        return 0
    fi

    [[ "$mode" == exec ]] && { chmod +x "$dst" 2>/dev/null || true; }

    local cur tpl
    cur=$(sciagent_sha1_file "$dst")
    tpl=$(sciagent_sha1_file "$src")

    if [[ "$cur" == "$tpl" ]]; then
        # Already current. Adopt it if we have no record, so that a body
        # materialized before the ownership records existed is reversible.
        if [[ ! -f "$state" && ! -f "$ceded" ]]; then
            mkdir -p "$(dirname "$state")"
            printf '%s\n' "$cur" > "$state"
        fi
        return 0
    fi

    # Already ceded to the user at this exact content — stay silent.
    if [[ -f "$ceded" && "$(cat "$ceded" 2>/dev/null)" == "$cur" ]]; then
        return 0
    fi

    local ours=false
    if [[ -f "$state" && "$(cat "$state" 2>/dev/null)" == "$cur" ]]; then
        ours=true
    elif ownership_template_hash_known "$cur" "$rel"; then
        ours=true
    fi

    if [[ "$ours" == true ]]; then
        cp "$src" "$dst" || { echo "sciagent: failed to refresh $dst" >&2; return 1; }
        [[ "$mode" == exec ]] && chmod +x "$dst"
        echo "refreshed: $dst (was an older toolkit version)"
        mkdir -p "$(dirname "$state")"
        printf '%s\n' "$tpl" > "$state"
        rm -f "$ceded"
        return 0
    fi

    echo "sciagent: $dst differs from the toolkit's version and was not written by sciagent" >&2
    echo "  leaving it as yours; sciagent will not manage it from now on" >&2
    rm -f "$state"
    mkdir -p "$(dirname "$ceded")"
    printf '%s\n' "$cur" > "$ceded"
    return 0
}

# ownership_teardown_body <dst> <state-file> <ceded-file>
#
# Reverse ONE ownership_ensure_body, using the saved content-hash snapshot:
# remove the file iff it is still exactly what we wrote/adopted; otherwise warn
# and leave it (see the ownership-record note above). The ceded marker is always
# dropped — it names a file that is explicitly NOT ours, so there is nothing to
# reverse, but the marker itself is our bookkeeping. Silent no-op when no state
# file exists (we never wrote that body, or it has already been reversed).
#
# Directories are never removed here; callers rmdir their own (never rm -r), so
# a non-empty directory means real content remains and is left in place.
ownership_teardown_body() {
    local dst="$1" state="$2" ceded="$3"
    rm -f "$ceded"
    [[ -f "$state" ]] || return 0
    local stored cur
    stored=$(cat "$state")
    if [[ -f "$dst" ]]; then
        cur=$(sciagent_sha1_file "$dst")
        if [[ "$cur" == "$stored" ]]; then
            rm -f "$dst"
        else
            echo "sciagent: warning — $dst was modified since sciagent created it; leaving it in place" >&2
        fi
    fi
    rm -f "$state"
    return 0
}
