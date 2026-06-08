#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Default: give full path explicitly, or relative to analysis/unfolding
INPUT="${1:-${BASE}/analysis/unfolding/out_embedding_BAYES_Inclusive_MCRC1p5/responses_embedding.root}"
OUT_DIR="${2:-${SCRIPT_DIR}/ptlead_comparison_plots}"

MACRO="${SCRIPT_DIR}/draw_ptlead_comparison.cxx"

echo "----------------------------------------"
echo "Drawing pTlead comparison"
echo "BASE    : ${BASE}"
echo "Input   : ${INPUT}"
echo "Macro   : ${MACRO}"
echo "Out dir : ${OUT_DIR}"
echo "----------------------------------------"

[[ -f "${INPUT}" ]] || { echo "ERROR: input not found: ${INPUT}"; exit 1; }
[[ -f "${MACRO}" ]] || { echo "ERROR: macro not found: ${MACRO}"; exit 1; }

mkdir -p "${OUT_DIR}"

root -l -b <<EOF
.L ${MACRO}+
draw_ptlead_comparison("${INPUT}","${OUT_DIR}");
gSystem->Exit(0);
EOF

echo "Done."