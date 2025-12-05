#!/usr/bin/env bash
set -euo pipefail

# --- edit these three if your layout changes ---
BASE="/gpfs/mnt/gpfs01/star/pwg/walz/Analysis_trees"
SIF="/gpfs/mnt/gpfs01/star/pwg/prozorov/singularity/roounfold.sif"
#SIF="/gpfs/mnt/gpfs01/star/pwg/walz/Analysis_trees/root_latest.sif"
MACRO="/gpfs/mnt/gpfs01/star/pwg/walz/Analysis_trees/analysis/unfolding/unfold_embedding_test.cxx"
# ----------------------------------------------

METHOD="${1:-Bayes}"
IN="${2:-$BASE/trees/merged_all/embedding_merged.root}"
OUT="${3:-$BASE/analysis/unfolding/out_test_Bayes_new}"

if [[ ! -f "$IN" ]]; then
  echo "Input not found: $IN"
  exit 1
fi




if [[ "$METHOD" != "Bayes" && "$METHOD" != "SVD" ]]; then
  echo "ERROR: Unknown method '$METHOD'. Use: Bayes or SVD"
  exit 1
fi

echo "Unfolding with method: $METHOD"

mkdir -p "$OUT"

#apptainer exec -e -B /gpfs/mnt/gpfs01 \
#  "$SIF" \
#  root -l -b -q "{gSystem->Load(\"libRooUnfold\"); .x $MACRO++(\"$IN\",\"$OUT\")}"
apptainer exec -e -B /gpfs/mnt/gpfs01 \
  "$SIF" \
  root -l -b -q \
  -e 'gSystem->Load("libRooUnfold");' \
  -e ".x ${MACRO}+(\"$IN\",\"$OUT\",\"$METHOD\")"

