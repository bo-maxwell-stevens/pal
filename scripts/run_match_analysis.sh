#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-/home/stevens/miniconda3/bin/python}"

"${PYTHON_BIN}" scripts/compare_fred_sources.py "$@"
