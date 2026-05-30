"""End-to-end engine loop with NO API (FakeConsensus) — exercises the full sequence."""

import pandas as pd


def test_engine_cellstate_no_api(tmp_path, fixtures_dir, skill_root, cassette, fake_consensus):
    from mllmct.engine import AnnotationEngine
    from mllmct.profile import load_profile

    profile = load_profile(skill_root / "profiles" / "cellstate.yaml")
    engine = AnnotationEngine(profile, git_cwd=skill_root)

    with fake_consensus(cassette):
        res = engine.run(
            markers_csv=str(fixtures_dir / "markers_synth.csv"),
            evidence_csv=str(fixtures_dir / "evidence_synth.csv"),
            aux_path=str(fixtures_dir / "programs_synth.csv"),
            species="mouse", tissue="synthetic in-vitro culture",
            models=["openai/gpt-5"], api_keys={"openrouter": "x"},
            lens="lensA", out_dir=str(tmp_path),
        )

    assert res.status == "PASS"
    df = pd.read_csv(res.labels_csv).set_index("cluster")

    # #1: Python-recomputed metrics differ from the (deliberately wrong) llm-reported ones.
    assert df.loc[3, "py_consensus_proportion"] == 0.5
    assert df.loc[3, "llm_reported_proportion"] == 1.0

    # #8/#10: harmonize ran — synonym mapped, novel preserved.
    assert df.loc[1, "harmonized_label"] == "Stressed"
    assert df.loc[3, "harmonized_label"].startswith("Novel:")
    assert df.loc[0, "harmonized_label"] == "Proliferating"

    # #6: trace tree complete.
    td = res.out_dir / "trace" / "lensA"
    for f in ("prompt.txt", "model_responses.json", "discussion.json", "tokens.json", "meta.json"):
        assert (td / f).exists(), f"missing trace artifact {f}"

    # evidence was injected into the prompt (cell-state path).
    prompt = (td / "prompt.txt").read_text()
    assert "Per-cluster evidence" in prompt and "program=P1[Biological" in prompt


def test_engine_celltype_markers_only(tmp_path, fixtures_dir, skill_root, cassette, fake_consensus):
    from mllmct.engine import AnnotationEngine
    from mllmct.profile import load_profile

    profile = load_profile(skill_root / "profiles" / "celltype.yaml")
    engine = AnnotationEngine(profile, git_cwd=skill_root)

    with fake_consensus(cassette):
        res = engine.run(
            markers_csv=str(fixtures_dir / "markers_synth.csv"),
            evidence_csv=None, aux_path=None,
            species="human", tissue="peripheral blood",
            models=["gemini-2.5-flash"], api_keys={"gemini": "x"},
            lens="pbmc", out_dir=str(tmp_path),
        )

    assert res.status == "PASS"
    prompt = (res.out_dir / "trace" / "pbmc" / "prompt.txt").read_text()
    # markers-only: no evidence section, markers present
    assert "Per-cluster evidence" not in prompt
    assert "MKI67" in prompt
    # cell-type open-vocab passthrough: free-text label kept verbatim (no Novel: prefix)
    df = pd.read_csv(res.labels_csv).set_index("cluster")
    assert df.loc[3, "harmonized_label"] == "some weird state"
