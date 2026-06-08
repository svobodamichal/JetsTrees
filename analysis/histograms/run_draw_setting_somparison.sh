#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$(cd "${SCRIPT_DIR}/../.." && pwd)"

UNFOLDING_DIR="${1:-${BASE}/analysis/unfolding}"
OUT_DIR="${2:-${SCRIPT_DIR}/setting_comparison_plots}"

MACRO="${SCRIPT_DIR}/draw_setting_comparison.cxx"
SIF="${BASE}/analysis/unfolding/roounfold.sif"

echo "----------------------------------------"
echo "Drawing setting comparison"
echo "BASE          : ${BASE}"
echo "Unfolding dir : ${UNFOLDING_DIR}"
echo "Macro         : ${MACRO}"
echo "Out dir       : ${OUT_DIR}"
echo "Reference     : out_embedding_BAYES_Inclusive_MCRC1p5"
echo "----------------------------------------"

[[ -d "${UNFOLDING_DIR}" ]] || { echo "ERROR: unfolding dir not found: ${UNFOLDING_DIR}"; exit 1; }
[[ -f "${MACRO}" ]] || { echo "ERROR: macro not found: ${MACRO}"; exit 1; }
[[ -f "${SIF}" ]] || { echo "ERROR: SIF file not found: ${SIF}"; exit 1; }

mkdir -p "${OUT_DIR}"

apptainer exec -e -B /gpfs01 \
  "${SIF}" \
  root -l -b <<EOF
gSystem->Load("libRooUnfold");
.L ${MACRO}+
draw_setting_comparison("${UNFOLDING_DIR}","${OUT_DIR}");
gSystem->Exit(0);
EOF

echo "----------------------------------------"
echo "Done."
echo "----------------------------------------"