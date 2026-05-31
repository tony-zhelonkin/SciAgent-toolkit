"""Shared offline fixtures. No API, no AnnData, no network."""

from __future__ import annotations

import json
from contextlib import contextmanager
from pathlib import Path

import pytest

FIX = Path(__file__).parent / "fixtures"
SKILL_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def fixtures_dir() -> Path:
    return FIX


@pytest.fixture
def skill_root() -> Path:
    return SKILL_ROOT


@pytest.fixture
def cassette() -> dict:
    return json.loads((FIX / "cassette_response.json").read_text())


@pytest.fixture
def fake_consensus():
    """Patch the name AS IMPORTED BY mllmct.engine (so a `from mllmcelltype import ...`
    local binding cannot escape the patch) to return a canned response. No network."""
    @contextmanager
    def _fake(response):
        import mllmct.engine as eng
        orig = eng.interactive_consensus_annotation
        eng.interactive_consensus_annotation = lambda **kwargs: response
        try:
            yield
        finally:
            eng.interactive_consensus_annotation = orig
    return _fake
