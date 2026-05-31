"""Non-invasive custom-prompt-template install seam (remediation fix #5).

``interactive_consensus_annotation`` does NOT thread ``prompt_template`` through to
``create_prompt`` (it hardcodes ``None`` → ``DEFAULT_PROMPT_TEMPLATE``). ``create_prompt``
reads ``mllmcelltype.prompts.DEFAULT_PROMPT_TEMPLATE`` from the module global AT CALL TIME,
so swapping that global is a clean, documented seam to inject a custom template without
editing the vendored library. This module is the SOLE owner of the snapshot sentinel so
install/restore is unambiguous. Ported from ``cellstate_llm``.
"""

from __future__ import annotations

_ORIGINAL_TEMPLATE: str | None = None


def install_prompt_template(template: str) -> None:
    """Monkeypatch the vendored ``DEFAULT_PROMPT_TEMPLATE`` with ``template``.

    Idempotent on the saved original (the first install snapshots the real default so
    ``restore_prompt_template`` can undo it). Call per run so each run gets its own baked-in
    template.
    """
    global _ORIGINAL_TEMPLATE
    import mllmcelltype.prompts as _prompts

    if _ORIGINAL_TEMPLATE is None:
        _ORIGINAL_TEMPLATE = _prompts.DEFAULT_PROMPT_TEMPLATE
    _prompts.DEFAULT_PROMPT_TEMPLATE = template


def restore_prompt_template() -> None:
    """Restore the vendored ``DEFAULT_PROMPT_TEMPLATE`` (best-effort, idempotent)."""
    global _ORIGINAL_TEMPLATE
    if _ORIGINAL_TEMPLATE is None:
        return
    import mllmcelltype.prompts as _prompts

    _prompts.DEFAULT_PROMPT_TEMPLATE = _ORIGINAL_TEMPLATE


def render_prompt_preview(
    species: str,
    tissue: str,
    template: str,
    marker_dict: dict[str, list[str]],
) -> str:
    """Render the EXACT prompt the library will produce, for trace/verification.

    Uses the library's own ``create_prompt`` so the preview is byte-faithful to what the
    model receives. Independent of ``install_prompt_template`` — it passes ``template``
    explicitly — so it is safe to call without mutating the global.
    """
    from mllmcelltype.prompts import create_prompt

    return create_prompt(
        marker_genes=marker_dict,
        species=species,
        tissue=tissue,
        prompt_template=template,
    )
