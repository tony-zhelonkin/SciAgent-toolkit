"""mllmct — version-locked multi-LLM consensus cell-type / cell-state annotation.

Public surface is intentionally small. Most callers use the CLI (`mllmct ...`); the
Python API below is for embedding. Internals (token capture, determinism, trace,
metrics) are in ``mllmct.core`` and are documented in ``references/`` — you do not
need to understand them to use the tool.
"""

from ._version import __version__

__all__ = ["__version__"]
