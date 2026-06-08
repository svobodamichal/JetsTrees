#!/bin/bash

# Stop on errors
set -e

# Directory where *this* script lives
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default file names; can be overridden from command line:
#   ./run_make_hists_pthat_diagnostics.sh embedding_merged.root hists_pthat_diagnostics.root
IN_BASENAME="${1:-embedding_merged_NoVeto.root}"
OUT_BASENAME="${2:-hists_pthat_diagnostics.root}"

# Paths relative to the script location
INPUT="${SCRIPT_DIR}/../../trees/${IN_BASENAME}"
OUTPUT="${SCRIPT_DIR}/${OUT_BASENAME}"
MACRO="${SCRIPT_DIR}/make_hists_pthat_diagnostics.C"

echo "----------------------------------------"
echo "Running pThat diagnostic histogram production"
echo "Script dir : ${SCRIPT_DIR}"
echo "Input      : ${INPUT}"
echo "Output     : ${OUTPUT}"
echo "Macro      : ${MACRO}"
echo "----------------------------------------"

if [ ! -f "${INPUT}" ]; then
  echo "ERROR: input file does not exist:"
  echo "  ${INPUT}"
  exit 1
fi

if [ ! -f "${MACRO}" ]; then
  echo "ERROR: macro does not exist:"
  echo "  ${MACRO}"
  exit 1
fi

root -l -b -q "${MACRO}(\"${INPUT}\",\"${OUTPUT}\")"

echo "----------------------------------------"
echo "Done."
echo "ROOT output:"
echo "  ${OUTPUT}"
echo "PDF/PNG output should be under:"
echo "  ${SCRIPT_DIR}/pdf_pthat_diagnostics/"
echo "----------------------------------------"
