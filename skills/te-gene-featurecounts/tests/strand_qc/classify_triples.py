#!/usr/bin/env python3
"""Back-compat shim — the regime classifier now lives in the runnable QC suite.

The single source of truth for the (s0,s2,s1) status-triple -> regime classifier is
``qc/tools/05_classify_triples.py`` (the SHARED module the runnable suite and this regression
both use). This shim re-exports ``classify`` and ``main`` from there so any caller that still
imports ``tests/strand_qc/classify_triples.py`` keeps working without a duplicated copy of the
regime logic. The regression (run_regression.sh) imports the suite module directly.
"""
import importlib.util
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_SHARED = os.path.normpath(os.path.join(_HERE, "..", "..", "qc", "tools", "05_classify_triples.py"))

_spec = importlib.util.spec_from_file_location("classify_triples_shared", _SHARED)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)

classify = _mod.classify
load_joint_counts = _mod.load_joint_counts
main = _mod.main

if __name__ == "__main__":
    sys.exit(main(sys.argv))
