#!/usr/bin/env bash
set -euo pipefail

########################
# Configuration
########################

# Directory where this script lives
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Base dir = two levels up from unfolding/
# (i.e. .../JetsTrees, assuming the same structure)
BASE="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Singularity image and macro, both in the same dir as this script
SIF="${SCRIPT_DIR}/roounfold.sif"
MACRO="${SCRIPT_DIR}/unfold_data.cxx"

########################
# Arguments
########################

# 1st arg: data input (either absolute path or basename under ${BASE}/trees)
if [[ $# -ge 1 ]]; then
  if [[ "$1" = /* ]]; then
    INPUT="$1"
  else
    INPUT="${BASE}/trees/$1"
  fi
else
  INPUT="${BASE}/trees/data_merged.root"
fi

# 2nd arg: RESPONSE ROOT FILE (single file with all tag directories)
RESP_FILE="${2:-${SCRIPT_DIR}/out_embedding/responses_embedding.root}"

# 3rd arg: EFFICIENCIES ROOT FILE
# default: analysis/efficiencies/efficiencies.root (based on your folder layout)
EFF_FILE="${3:-${BASE}/analysis/efficiencies/efficiencies.root}"

# 4th arg: output directory for unfolded data spectra
OUT_DIR="${4:-${SCRIPT_DIR}/out_data}"

# 5th arg: number of Bayes iterations
NITER="${5:-4}"

########################
# Checks
########################

echo "----------------------------------------"
echo "Running unfolding on REAL DATA"
echo "SCRIPT_DIR  : $SCRIPT_DIR"
echo "BASE        : $BASE"
echo "SIF         : $SIF"
echo "Macro       : $MACRO"
echo "Input data  : $INPUT"
echo "Resp. file  : $RESP_FILE"
echo "Eff. file   : $EFF_FILE"
echo "Output dir  : $OUT_DIR"
echo "Iterations  : $NITER"
echo "----------------------------------------"

[[ -f "$SIF"       ]] || { echo "ERROR: SIF not found:       $SIF";       exit 1; }
[[ -f "$MACRO"     ]] || { echo "ERROR: MACRO not found:     $MACRO";     exit 1; }
[[ -f "$INPUT"     ]] || { echo "ERROR: Input not found:     $INPUT";     exit 1; }
[[ -f "$RESP_FILE" ]] || { echo "ERROR: Resp. file not found: $RESP_FILE"; exit 1; }
[[ -f "$EFF_FILE"  ]] || { echo "ERROR: Eff. file not found: $EFF_FILE"; exit 1; }

mkdir -p "$OUT_DIR"

########################
# Run inside container
########################

# IMPORTANT: bind /gpfs01 so the container sees the same paths
apptainer exec -e -B /gpfs01 \
  "$SIF" \
  root -l -b <<EOF
gSystem->Load("libRooUnfold");
.x ${MACRO}+("${INPUT}","${RESP_FILE}","${EFF_FILE}","${OUT_DIR}",${NITER});
.q
EOF

echo "----------------------------------------"
echo "Done."
echo "----------------------------------------"
