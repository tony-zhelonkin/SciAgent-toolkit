"""Byte-faithful custom-prompt-template rendering for the trace.

``interactive_consensus_annotation`` accepts ``prompt_template`` natively
(``mllmcelltype>=2.0.7``) and threads it through to ``create_prompt``, so the engine
passes the filled template straight into the call — no ``DEFAULT_PROMPT_TEMPLATE``
module-global monkeypatch (that was seam (i), retired at the 2.0.7 bump; see
``references/monkeypatch-internals.md``). This module keeps only the preview renderer,
which calls the library's own ``create_prompt`` so the trace records the exact prompt the
model receives.
"""

from __future__ import annotations


def render_prompt_preview(
    species: str,
    tissue: str,
    template: str,
    marker_dict: dict[str, list[str]],
) -> str:
    """Render the EXACT prompt the library will produce, for trace/verification.

    Uses the library's own ``create_prompt`` with ``prompt_template=template`` — the same
    template the engine threads into ``interactive_consensus_annotation`` — so the preview
    is byte-faithful to what the model receives and never mutates library state.
    """
    from mllmcelltype.prompts import create_prompt

    return create_prompt(
        marker_genes=marker_dict,
        species=species,
        tissue=tissue,
        prompt_template=template,
    )
