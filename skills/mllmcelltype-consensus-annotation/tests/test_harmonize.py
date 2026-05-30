"""#8/#10 — open/closed vocabulary harmonization + guard ordering."""

from mllmct.core.harmonize import harmonize_label, is_novel_label, normalize_novel_token

VOCAB = ["ISG_High", "Apoptotic", "Ambiguous_LowSignal"]
SYN = {"interferon": "ISG_High", "apoptot": "Apoptotic"}


def test_exact_case_insensitive():
    assert harmonize_label("isg_high", VOCAB, SYN) == "ISG_High"


def test_longest_synonym_substring():
    assert harmonize_label("strong interferon response", VOCAB, SYN) == "ISG_High"


def test_open_mode_preserves_novel():
    out = harmonize_label("metabolic rewiring state", VOCAB, SYN, mode="open")
    assert is_novel_label(out)
    assert out == "Novel:metabolic_rewiring_state"


def test_closed_mode_snaps_to_fallback():
    out = harmonize_label("metabolic rewiring", VOCAB, SYN, mode="closed",
                          fallback_label="Ambiguous_LowSignal")
    assert out == "Ambiguous_LowSignal"


def test_blank_is_fallback():
    assert harmonize_label("   ", VOCAB, SYN, fallback_label="X") == "X"


def test_celltype_verbatim_passthrough():
    # cell-type config: empty vocab, no prefix, no cleaning → exact label kept
    assert harmonize_label("CD8+ cytotoxic T cells", [], {}, mode="open",
                           novel_prefix="", clean_novel=False) == "CD8+ cytotoxic T cells"


def test_guard_rejects_before_novel():
    from mllmct.core.guards import make_reject_guard
    guard = make_reject_guard([r"\bth1\b"])
    # a guarded label must NOT become Novel:* even in open mode
    out = harmonize_label("Th1 cells", VOCAB, SYN, guard_fn=guard, mode="open",
                          fallback_label="Ambiguous_LowSignal")
    assert out == "Ambiguous_LowSignal"


def test_normalize_novel_token():
    assert normalize_novel_token("DNA Damage Response cells") == "DNA_Damage_Response"
    assert normalize_novel_token("Stress / Recovery cells") == "Stress_Recovery"
