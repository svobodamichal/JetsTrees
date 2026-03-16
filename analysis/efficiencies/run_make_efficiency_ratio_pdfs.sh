#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Defaults match your layout from the screenshot
INDIR="${1:-$SCRIPT_DIR}"                      # where efficiencies_*.root live
OUTDIR="${2:-$SCRIPT_DIR/PDFs_Ratios}"         # output folder for PDFs
REFPAT="${3:-50_-1}"                           # substring identifying the reference file

MACRO="${SCRIPT_DIR}/make_efficiency_ratio_pdfs.C"

echo "----------------------------------------"
echo "Efficiency ratio PDFs"
echo "indir   : $INDIR"
echo "outdir  : $OUTDIR"
echo "refpat  : $REFPAT"
echo "macro   : $MACRO"
echo "----------------------------------------"

root -l -b -q "$MACRO(\"$INDIR\",\"$OUTDIR\",\"$REFPAT\")"

echo "----------------------------------------"
echo "Done."
echo "----------------------------------------"