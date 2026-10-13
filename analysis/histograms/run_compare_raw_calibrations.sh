#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
# Relative paths are resolved against this wrapper's directory.
resolve_path() {
  case "$1" in /*) printf '%s' "$1" ;; *) printf '%s/%s' "$SCRIPT_DIR" "$1" ;; esac
}
ORIG="$(resolve_path "${1:-hists_data_OrigCalib.root}")"
HS="$(resolve_path "${2:-hists_data_HSCalib.root}")"
OUT="$(resolve_path "${3:-raw_calibration_comparison}")"
REBIN="${4:-1}"
XMIN="${5:--20}"
XMAX="${6:-60}"
[[ $# -le 6 ]] || { echo "Usage: $0 [ORIG.root [HS.root [OUT_DIR [REBIN [XMIN [XMAX]]]]]]" >&2; exit 1; }
[[ "$REBIN" =~ ^[1-9][0-9]*$ ]] || { echo 'REBIN must be a positive integer' >&2; exit 1; }
for value in "$XMIN" "$XMAX"; do
  [[ "$value" =~ ^-?[0-9]+([.][0-9]+)?$ ]] || { echo 'X limits must be numbers' >&2; exit 1; }
done
command -v root >/dev/null || { echo 'ROOT is not available. Run in your ROOT environment/container.' >&2; exit 1; }
[[ -f "$ORIG" && -f "$HS" ]] || { echo "Missing input: $ORIG or $HS" >&2; exit 1; }
# Escape paths for C++ string literals, preserving spaces and shell metacharacters.
cpp_quote() { local s="$1"; s="${s//\\/\\\\}"; s="${s//\"/\\\"}"; printf '"%s"' "$s"; }
MACRO="$SCRIPT_DIR/compare_raw_calibrations.C"
printf 'Original: %s\nHS:       %s\nOutput:   %s\n' "$ORIG" "$HS" "$OUT"
# ACLiC compiles C++03 with the ROOT 5 toolchain, avoiding CINT execution.
root -l -b -q "$MACRO+($(cpp_quote "$ORIG"),$(cpp_quote "$HS"),$(cpp_quote "$OUT"),$REBIN,$XMIN,$XMAX)"
