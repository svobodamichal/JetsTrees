#!/usr/bin/env bash
set -euo pipefail

########################
# Configuration
########################

# Base path as seen *inside the container*

BASE="/gpfs/mnt/gpfs01/star/pwg/svomich/JetsTrees"

SIF="${BASE}/analysis/unfolding/roounfold.sif"

# Directory where this script lives (macro is assumed to be here as well)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MACRO="${BASE}/analysis/unfolding/unfold_data.cxx"

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

# 2nd arg: directory with response files (unfold_response_*.root)
RESP_DIR="${2:-${BASE}/analysis/unfolding/out}"

# 3rd arg: output directory for unfolded data spectra
OUT_DIR="${3:-${BASE}/analysis/unfolding/out_data}"

# 4th arg: number of Bayes iterations
NITER="${4:-4}"

########################
# Checks
########################

echo "----------------------------------------"
echo "Running unfolding on REAL DATA"
echo "SIF        : $SIF"
echo "Macro      : $MACRO"
echo "Input data : $INPUT"
echo "Resp. dir  : $RESP_DIR"
echo "Output dir : $OUT_DIR"
echo "Iterations : $NITER"
echo "----------------------------------------"

[[ -f "$SIF"    ]] || { echo "ERROR: SIF not found:    $SIF";    exit 1; }
[[ -f "$MACRO"  ]] || { echo "ERROR: MACRO not found:  $MACRO";  exit 1; }
[[ -f "$INPUT"  ]] || { echo "ERROR: Input not found:  $INPUT";  exit 1; }
[[ -d "$RESP_DIR" ]] || { echo "ERROR: Resp. dir not found: $RESP_DIR"; exit 1; }

mkdir -p "$OUT_DIR"

########################
# Run inside container
########################

apptainer exec -e -B /gpfs/mnt/gpfs01 \
  "$SIF" \
  bash -lc "cd '$SCRIPT_DIR'; root -l -b <<EOF
gSystem->Load(\"libRooUnfold\");
.x ${MACRO}+(\"${INPUT}\",\"${RESP_DIR}\",\"${OUT_DIR}\",${NITER});
.q
EOF"

echo "----------------------------------------"
echo "Done."
echo "----------------------------------------"
