#!/usr/bin/env bash
# Run the peak-atlas-multiome test scaffolds under testthat.
#
# Convention (matches the toolkit's packaged-skill runners): if R or testthat is
# unavailable, SKIP with exit 0 so CI / tests/run-all.sh stays green. The tests
# are SCAFFOLDS — they report as skipped (via skip_scaffold()) until their
# fixtures are implemented in an R + Bioconductor env. See tests/README.md.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if ! command -v Rscript >/dev/null 2>&1; then
  echo "SKIP: Rscript not found — peak-atlas tests need R + Bioconductor."
  exit 0
fi
if ! Rscript -e 'if (!requireNamespace("testthat", quietly=TRUE)) quit(status=1)' >/dev/null 2>&1; then
  echo "SKIP: R package 'testthat' is not installed."
  exit 0
fi

Rscript -e "testthat::test_dir('${here}/testthat', stop_on_failure = TRUE)"
