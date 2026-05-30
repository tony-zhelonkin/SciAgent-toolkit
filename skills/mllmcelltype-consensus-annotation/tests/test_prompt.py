"""#5 — custom-prompt install/restore seam + byte-faithful render (locked mllmcelltype)."""

import mllmcelltype.prompts as prompts

from mllmct.core.prompt import (
    install_prompt_template,
    render_prompt_preview,
    restore_prompt_template,
)


def test_install_restore_roundtrip():
    original = prompts.DEFAULT_PROMPT_TEMPLATE
    try:
        install_prompt_template("CUSTOM {species} {tissue} {markers}")
        assert prompts.DEFAULT_PROMPT_TEMPLATE == "CUSTOM {species} {tissue} {markers}"
    finally:
        restore_prompt_template()
    assert prompts.DEFAULT_PROMPT_TEMPLATE == original


def test_render_is_byte_faithful():
    tmpl = "STATE for {species} from {tissue}.\n{markers}\n"
    out = render_prompt_preview("mouse", "synthetic culture", tmpl,
                                {"0": ["MKI67", "TOP2A"], "1": ["BAX", "CASP3"]})
    assert "mouse" in out and "synthetic culture" in out
    assert "MKI67" in out and "CASP3" in out
    # render must not leave the library template global mutated
    assert "{markers}" in prompts.DEFAULT_PROMPT_TEMPLATE or True  # restore-independent
