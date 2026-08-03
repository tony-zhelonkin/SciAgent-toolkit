"""#5 — byte-faithful custom-prompt render via the native prompt_template path."""

import mllmcelltype.prompts as prompts

from mllmct.core.prompt import render_prompt_preview


def test_render_is_byte_faithful():
    tmpl = "STATE for {species} from {tissue}.\n{markers}\n"
    out = render_prompt_preview("mouse", "synthetic culture", tmpl,
                                {"0": ["MKI67", "TOP2A"], "1": ["BAX", "CASP3"]})
    assert "mouse" in out and "synthetic culture" in out
    assert "MKI67" in out and "CASP3" in out


def test_render_does_not_mutate_library_default():
    """render_prompt_preview passes prompt_template explicitly, so the library default global
    is untouched (no install/restore seam anymore)."""
    before = prompts.DEFAULT_PROMPT_TEMPLATE
    render_prompt_preview("mouse", "lung", "CUSTOM {species} {tissue} {markers}",
                          {"0": ["CD3E"]})
    assert prompts.DEFAULT_PROMPT_TEMPLATE == before
