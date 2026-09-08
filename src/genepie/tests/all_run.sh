#!/bin/bash
# Run the whole genepie test suite with pytest (all arguments are forwarded).
#
#   bash src/genepie/tests/all_run.sh            # everything, incl. slow tests
#   bash src/genepie/tests/all_run.sh -m "not slow"
#
# Optional datasets self-skip when absent:
#   python -m genepie.tests.download_test_data    # chignolin (integration tests)
#   python -m genepie.tests.download_tremd_data   # T-REMD (MBAR resampling)
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
exec python -m pytest "$SCRIPT_DIR" "$@"
