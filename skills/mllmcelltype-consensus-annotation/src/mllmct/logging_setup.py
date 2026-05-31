"""CLI logging: human-readable to stderr (stdout stays parseable) + optional file."""

from __future__ import annotations

import logging
from pathlib import Path


def configure_cli_logging(level: str = "INFO", log_file: str | Path | None = None) -> None:
    root = logging.getLogger()
    root.setLevel(getattr(logging, level.upper(), logging.INFO))
    # Clear our own prior handlers (idempotent across re-invocations in one process).
    for h in list(root.handlers):
        root.removeHandler(h)

    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s", "%H:%M:%S")

    stderr = logging.StreamHandler()  # defaults to stderr
    stderr.setFormatter(fmt)
    root.addHandler(stderr)

    if log_file:
        fh = logging.FileHandler(Path(log_file))
        fh.setFormatter(fmt)
        root.addHandler(fh)
