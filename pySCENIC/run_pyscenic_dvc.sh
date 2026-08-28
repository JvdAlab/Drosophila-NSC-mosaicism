#!/usr/bin/env bash
# Script to run the pySCENIC DVC pipeline end to end
set -euo pipefail

# Run from the repository root regardless of the caller's current directory.
cd "$(dirname "${BASH_SOURCE[0]}")/.."

echo "Reproducing the pySCENIC DVC pipeline from the submission directory"
DVC_BIN="${DVC_BIN:-dvc}"
if ! command -v "$DVC_BIN" >/dev/null 2>&1 && [ -x "${HOME}/bin/dvc" ]; then
  DVC_BIN="${HOME}/bin/dvc"
fi
"$DVC_BIN" repro "$@" pySCENIC/dvc.yaml
