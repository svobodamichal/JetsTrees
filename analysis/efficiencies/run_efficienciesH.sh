#!/bin/bash

# Stop on errors
set -e

# Directory where *this* script lives
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default file names (can be overridden via arguments)
# 1st arg: input merged file (old, with histos)
# 2nd arg: output efficiencies file
IN_BASENAME="${1:-embedding_mergedH.root}"
OUT_BASENAME="${2:-efficienciesH.root}"

# Paths
INPUT="${SCRIPT_DIR}/../../trees/${IN_BASENAME}"
OUTPUT="${SCRIPT_DIR}/${OUT_BASENAME}"
MACRO="${SCRIPT_DIR}/make_efficienciesH.C"

echo "----------------------------------------"
echo "Running efficiency production (hist-based)"
echo "Script dir : $SCRIPT_DIR"
echo "Input      : $INPUT"
echo "Output     : $OUTPUT"
echo "Macro      : $MACRO"
echo "----------------------------------------"

root -l -b -q "$MACRO(\"$INPUT\",\"$OUTPUT\")"

echo "----------------------------------------"
echo "Done."
echo "----------------------------------------"