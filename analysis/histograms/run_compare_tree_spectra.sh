#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

TREE_DIR="${1:-${SCRIPT_DIR}/../../trees}"
OUT_DIR="${2:-${SCRIPT_DIR}/comparison_tree_spectra}"
MACRO="${SCRIPT_DIR}/compare_tree_spectra.C"

echo "----------------------------------------"
echo "Comparing tree spectra"
echo "Script dir : ${SCRIPT_DIR}"
echo "Tree dir   : ${TREE_DIR}"
echo "Output dir : ${OUT_DIR}"
echo "Macro      : ${MACRO}"
echo "----------------------------------------"

root -l -b -q "${MACRO}(\"${TREE_DIR}\",\"${OUT_DIR}\")"

echo "----------------------------------------"
echo "Done."
echo "----------------------------------------"