"""Guard hook + restricted-state guardrails."""

from mllmct.core.guards import check_restricted_states, make_reject_guard


def test_reject_guard_whole_word():
    g = make_reject_guard([r"\bth1\b", r"\bth17\b", r"\bcd4\b"])
    assert g("th1 cells")
    assert g("th17")
    assert g("cd4 t-helper")
    # a genuine state that merely ends in 'cells' must survive
    assert not g("dna damage response cells")


def test_empty_patterns_never_fire():
    g = make_reject_guard([])
    assert not g("anything at all")


def test_restricted_state_violation():
    restricted = {"StateX": ["GroupB"]}
    warns = check_restricted_states("StateX", "GroupA", restricted)
    assert warns and "RESTRICTED-STATE" in warns[0]


def test_restricted_state_allowed():
    restricted = {"StateX": ["GroupB"]}
    assert check_restricted_states("StateX", "GroupB", restricted) == []
    assert check_restricted_states("Unlisted", "GroupA", restricted) == []
