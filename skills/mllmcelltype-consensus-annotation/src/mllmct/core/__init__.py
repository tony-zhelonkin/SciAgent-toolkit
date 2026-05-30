"""mllmct.core — domain-agnostic remediation library.

Nothing in this package names a cell type, lineage, tissue, or program. It holds the
reusable machinery distilled from a real multi-LLM annotation run:

  - capture       : recover token usage the library discards + force determinism
  - normalize     : the shared token-usage shape
  - consensus     : Python-recomputed consensus_proportion + Shannon entropy (#1)
  - debug_capture : keep library DEBUG capture alive across the consensus call (#4)
  - prompt        : non-invasive custom-prompt-template install/restore seam (#5)
  - harmonize     : open/closed vocabulary harmonization with a pluggable guard (#8/#10)
  - guards        : build guard functions + post-hoc label guardrails
  - trace         : self-describing per-run trace tree + idempotent token non-clobber (#6/#10)
  - cost          : logging-only USD cost estimate + aggregate summary (#7)
  - evidence      : EvidenceProvider ABC — the "inject extra metadata" extension point
"""
