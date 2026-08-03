"""Single version source-of-truth (kept in sync with pyproject [project].version)."""

__version__ = "0.3.0"

# Exact pins this build's monkeypatches are validated against. checks/smoke_check_versions.py
# asserts the installed versions match these AND that the structural seams still exist.
PINNED = {
    "mllmcelltype": "2.0.7",
    "google-genai": "2.6.0",
    "pydantic": "2.13.3",
    "requests": "2.33.1",
}
