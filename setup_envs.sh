#!/usr/bin/env bash
#
# Create every conda environment used in this repository.
#
# There is no single combined environment
# See the "Environments" table in README.md. Run this from the repository root.

set -euo pipefail

MODULES=(scRNA_pipeline edgeR hdWGCNA metabolomics statistical_analysis reports pySCENIC)

for m in "${MODULES[@]}"; do
    echo "=== conda env create -f ${m}/environment.yml ==="
    conda env create -f "${m}/environment.yml"
    echo
done

echo "All environments created. Activate one with 'conda activate <name>'"
echo "(see the Environments table in README.md)."
