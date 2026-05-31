"""Keep library DEBUG capture alive across the consensus call (remediation fix #4).

mLLMCelltype logs raw model responses only at DEBUG, and it re-initializes its logger to
INFO once PER MODEL inside the consensus call (``annotate.setup_logging`` → the
"just update level" fast path), which downgrades any DEBUG handler we attached BEFORE the
first provider call fires — so 0 raw-response lines land. ``llm_debug_capture`` wraps
``mllmcelltype.logger.setup_logging`` (and every module that imported it by value) so we
RE-ASSERT DEBUG after each library re-init. Reversible; the vendored package is not edited.

Decoupled from the original project paths: ``log_dir`` is an explicit argument here
(was ``PATHS.annotation/llm_trace/<lens>``). Ported from ``cellstate_obs``.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

log = logging.getLogger(__name__)

# Name of the vendored mllmcelltype logger (mllmcelltype/logger.py).
LLMCELLTYPE_LOGGER_NAME = "llmcelltype"


def _handler_name(lens: str) -> str:
    """Stable tag for the per-lens DEBUG FileHandler (identifiable/removable)."""
    return f"mllmct_debug::{lens}"


def _ensure_debug_handler(lens: str, log_file: Path) -> logging.Handler:
    """Idempotently attach our DEBUG ``FileHandler`` to the ``llmcelltype`` logger.

    Returns the (existing or newly-created) handler and forces it + the logger to DEBUG.
    Re-callable any number of times — neither duplicates a handler pointing at the same
    file nor leaves the level downgraded.
    """
    lib_logger = logging.getLogger(LLMCELLTYPE_LOGGER_NAME)
    lib_logger.setLevel(logging.DEBUG)

    target = str(log_file.resolve())
    for h in lib_logger.handlers:
        if isinstance(h, logging.FileHandler) and getattr(h, "baseFilename", None) == target:
            h.setLevel(logging.DEBUG)
            return h

    handler = logging.FileHandler(log_file)
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(
        logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s",
                          "%Y-%m-%d %H:%M:%S")
    )
    handler.set_name(_handler_name(lens))
    lib_logger.addHandler(handler)
    return handler


def configure_llm_logging(lens: str, log_dir: Path | str) -> Path:
    """Attach a DEBUG ``FileHandler`` to the vendored ``llmcelltype`` logger.

    Raises the library logger to DEBUG so its raw-response DEBUG lines are persisted to a
    per-lens file, without editing the library. Idempotent per (lens, file). NOTE — this
    alone is NOT enough during a consensus run; wrap the call in ``llm_debug_capture`` so
    the DEBUG level survives the library's mid-call re-inits.

    Args:
        lens: lens / run name used in the file name.
        log_dir: directory for ``llmcelltype_debug.log``.

    Returns:
        Path to the per-lens debug log file.
    """
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "llmcelltype_debug.log"
    _ensure_debug_handler(lens, log_file)
    log.info("LLM DEBUG capture for lens '%s' -> %s", lens, log_file)
    return log_file


class llm_debug_capture:
    """Context manager that keeps DEBUG capture alive across the consensus run.

    Wraps ``mllmcelltype.logger.setup_logging`` so that, AFTER the library's own
    ``setup_logging`` runs (per model), we re-add our FileHandler if removed and force the
    logger + handler back to DEBUG. The original is restored on ``__exit__``.

    Usage::

        configure_llm_logging(lens, log_dir)
        with llm_debug_capture(lens, log_dir):
            result = interactive_consensus_annotation(...)
    """

    def __init__(self, lens: str, log_dir: Path | str) -> None:
        self.lens = lens
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.log_file = self.log_dir / "llmcelltype_debug.log"
        self._patched: list[tuple[Any, str, Any]] = []  # (target_obj, attr, original)

    def _reassert(self) -> None:
        """Re-add our handler if dropped and force logger+handler back to DEBUG."""
        _ensure_debug_handler(self.lens, self.log_file)

    def __enter__(self) -> "llm_debug_capture":
        _ensure_debug_handler(self.lens, self.log_file)

        try:
            import mllmcelltype.logger as _logger_mod
        except Exception as exc:  # pragma: no cover - library absent
            log.info("llm_debug_capture: mllmcelltype not importable (%s) — "
                     "DEBUG re-assert disabled", exc)
            return self

        orig_setup = _logger_mod.setup_logging
        guard = self

        def _wrapped_setup(*args, **kwargs):  # noqa: ANN002, ANN003
            result = orig_setup(*args, **kwargs)
            try:
                guard._reassert()
            except Exception as exc:  # pragma: no cover - defensive
                log.debug("llm_debug_capture: re-assert failed: %s", exc)
            return result

        # Patch the canonical definition AND every module that imported the name by value
        # (``annotate`` does ``from .logger import setup_logging``; the package re-exports it).
        targets: list[tuple[Any, str]] = [(_logger_mod, "setup_logging")]
        try:
            import mllmcelltype.annotate as _annotate_mod
            if getattr(_annotate_mod, "setup_logging", None) is orig_setup:
                targets.append((_annotate_mod, "setup_logging"))
        except Exception:  # pragma: no cover - defensive
            pass
        try:
            import mllmcelltype as _pkg
            if getattr(_pkg, "setup_logging", None) is orig_setup:
                targets.append((_pkg, "setup_logging"))
        except Exception:  # pragma: no cover - defensive
            pass

        for obj, attr in targets:
            self._patched.append((obj, attr, getattr(obj, attr)))
            setattr(obj, attr, _wrapped_setup)
        return self

    def __exit__(self, *exc) -> None:
        for obj, attr, original in self._patched:
            try:
                setattr(obj, attr, original)
            except Exception:  # pragma: no cover - defensive
                pass
        self._patched = []
        try:
            self._reassert()
        except Exception:  # pragma: no cover - defensive
            pass
