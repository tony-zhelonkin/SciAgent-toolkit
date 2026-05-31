"""#2/#3 — reference PanelEvidenceProvider: program decode, signature de-sat, axis cap, companion."""

import yaml

from mllmct.plugins.panel_evidence import PanelEvidenceProvider


def _provider(fixtures_dir, skill_root):
    opts = yaml.safe_load((skill_root / "profiles" / "cellstate.yaml").read_text())["evidence_options"]
    return PanelEvidenceProvider(options=opts, aux_path=str(fixtures_dir / "programs_synth.csv"))


def test_assemble_all_clusters(fixtures_dir, skill_root):
    ev = _provider(fixtures_dir, skill_root).assemble(str(fixtures_dir / "evidence_synth.csv"))
    assert set(ev) == {"0", "1", "2", "3", "4", "5"}


def test_program_decode_and_unmapped(fixtures_dir, skill_root):
    ev = _provider(fixtures_dir, skill_root).assemble(str(fixtures_dir / "evidence_synth.csv"))
    assert "program=P1[Biological:" in ev["0"]      # decoded, human-readable + genes
    assert "program=P9(unmapped)" in ev["3"]         # absent id surfaced, not silently dropped
    assert "Technical:P3(artifact)" in ev["2"]       # drop-category flagged, not shown as biology


def test_signature_desaturation(fixtures_dir, skill_root):
    ev = _provider(fixtures_dir, skill_root).assemble(str(fixtures_dir / "evidence_synth.csv"))
    # always-on signature excluded from ranked slot, surfaced as a flag
    assert "Housekeeping_sig" not in ev["0"].split("signatures=")[1].split("(")[0]
    assert "Housekeeping_sig+" in ev["0"]
    # cluster 2 has ONLY the always-on signature → none-discriminative + flag
    assert "none-discriminative" in ev["2"] and "Housekeeping_sig+" in ev["2"]


def test_axis_cap_and_companion(fixtures_dir, skill_root):
    ev = _provider(fixtures_dir, skill_root).assemble(str(fixtures_dir / "evidence_synth.csv"))
    # cluster 2 has 4 non-baseline axes but max_axes=3 → activation (lowest |z|) dropped
    assert "activation=" not in ev["2"]
    assert "proliferation=" in ev["2"] and "stress_response=" in ev["2"]
    # apoptosis is the companion axis → rendered with an inline companion bin
    assert "apoptosis=High companion=High" in ev["2"]
    # cluster 1: stress_response (trigger) fires → apoptosis force-rendered though it is Mid
    assert "apoptosis=Mid companion=High" in ev["1"]


def test_flag_token(fixtures_dir, skill_root):
    ev = _provider(fixtures_dir, skill_root).assemble(str(fixtures_dir / "evidence_synth.csv"))
    assert "proliferating" in ev["0"]   # flag_metric 0.5 > 0.05
    assert "resting" in ev["1"]         # flag_metric 0.01 <= 0.05
