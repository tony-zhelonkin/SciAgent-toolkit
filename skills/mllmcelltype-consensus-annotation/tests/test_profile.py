"""Profile loading: built-ins, custom, validation, inject_evidence toggle."""

import pytest

from mllmct.profile import load_profile


def test_load_celltype(skill_root):
    p = load_profile(skill_root / "profiles" / "celltype.yaml")
    assert p.name == "celltype"
    assert p.inject_evidence is False
    assert p.novel_prefix == ""        # verbatim passthrough for cell-type
    assert p.prompt_template.exists()


def test_load_cellstate(skill_root):
    p = load_profile(skill_root / "profiles" / "cellstate.yaml")
    assert p.inject_evidence is True
    assert p.evidence_provider.endswith("PanelEvidenceProvider")
    assert "Proliferating" in p.vocab
    assert p.prompt_template.exists()


def test_load_custom(fixtures_dir, skill_root):
    p = load_profile(fixtures_dir / "profile_custom.yaml", skill_root=skill_root)
    assert p.name == "custom_test"
    assert p.vocab_mode == "closed"
    assert p.join["enabled"] is True


def test_missing_name_errors(tmp_path):
    bad = tmp_path / "bad.yaml"
    bad.write_text("inject_evidence: false\n")
    with pytest.raises(ValueError):
        load_profile(bad)


def test_inject_without_provider_errors(tmp_path):
    bad = tmp_path / "bad.yaml"
    bad.write_text("name: x\ninject_evidence: true\nevidence_provider: null\n")
    with pytest.raises(ValueError):
        load_profile(bad)


def test_bad_vocab_mode_errors(tmp_path):
    bad = tmp_path / "bad.yaml"
    bad.write_text("name: x\nvocab_mode: sideways\n")
    with pytest.raises(ValueError):
        load_profile(bad)
