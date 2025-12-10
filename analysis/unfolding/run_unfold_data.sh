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
#MACRO="${SCRIPT_DIR}/unfold_data_test.cxx"
MACRO="$(realpath ${SCRIPT_DIR}/unfold_data_test.cxx)"

########################
# Arguments
########################


#1st arg: method of unfolding (Bayes or SVD) - default: Bayes

METHOD="${1:-Bayes}"

# 2nd arg: data input (either absolute path or basename under ${BASE}/trees)

#if [[ $# -ge 2 ]]; then
#  if [[ "$2" = /* ]]; then
#    INPUT="$2"
#  else
#    INPUT="${BASE}/trees/$2"
#  fi
#else
#  INPUT="${BASE}/trees/data_merged.root"
#fi

INPUT="${2:-/gpfs/mnt/gpfs01/star/pwg/svomich/JetsTrees/trees/data_merged.root}"

# 3rd arg: RESPONSE ROOT FILE (single file with all tag directories)
# adjust the default name to whatever you actually wrote in the embedding macro
RESP_FILE="${3:-${SCRIPT_DIR}/out_SVD_Michalova_data/responses_embedding.root}"

# 4th arg: output directory for unfolded data spectra
OUT_DIR="${4:-${SCRIPT_DIR}/out_data_Bayes}"

# 5th arg: number of Bayes iterations
NITER="${5:-4}"

########################
# Checks
########################

echo "----------------------------------------"
echo "Running unfolding on REAL DATA"
echo Script dir : $SCRIPT_DIR
echo "Method      : $METHOD"
echo "SIF         : $SIF"
echo "Macro       : $MACRO"
echo "Input data  : $INPUT"
echo "Resp. file  : $RESP_FILE"
echo "Output dir  : $OUT_DIR"
echo "Iterations for Bayes : $NITER"
echo "----------------------------------------"

[[ -f "$SIF"      ]] || { echo "ERROR: SIF not found:      $SIF";      exit 1; }
[[ -f "$MACRO"    ]] || { echo "ERROR: MACRO not found:    $MACRO";    exit 1; }
[[ -f "$INPUT"    ]] || { echo "ERROR: Input not found:    $INPUT";    exit 1; }
[[ -f "$RESP_FILE" ]] || { echo "ERROR: Resp. file not found: $RESP_FILE"; exit 1; }

mkdir -p "$OUT_DIR"

########################
# Run inside container
########################

apptainer exec -e -B /gpfs/mnt/gpfs01,/gpfs01 \
  "$SIF" \
  bash -lc "root -l -b <<EOF
gSystem->Load(\"libRooUnfold\");
.x ${MACRO}+(\"${INPUT}\",\"${RESP_FILE}\",\"${OUT_DIR}\",${NITER},\"$METHOD\");
.q
EOF"

echo "----------------------------------------"
echo "Done."
echo "----------------------------------------"

