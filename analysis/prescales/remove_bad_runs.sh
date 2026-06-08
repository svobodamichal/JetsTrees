#!/bin/bash

# Input files
bht2_matched="matched_output_BHT2VPDMB.txt"
vpdmb_matched="matched_output_VPDMB.txt"
bad_run_list="BadRunList_14.list"

# Output files (new versions of matched files without bad runs)
bht2_output="BHT2VPDMB-30_matched_cleaned.txt"
vpdmb_output="VPDMB-30_matched_cleaned.txt"

# Read the bad run numbers into an array
bad_runs=($(cat "$bad_run_list"))

# Function to remove bad run numbers from a file
remove_bad_runs() {
  input_file="$1"
  output_file="$2"

  # Copy the header (first row) to the output file
  head -n 1 "$input_file" > "$output_file"

  # Iterate over each line in the input file (skipping the header)
  tail -n +2 "$input_file" | while read -r line; do
    # Extract the run number (first column)
    run_number=$(echo "$line" | awk '{print $1}')

    # Check if the run number is in the list of bad runs
    if [[ ! " ${bad_runs[*]} " =~ " $run_number " ]]; then
      # If it's not a bad run, append the line to the output file
      echo "$line" >> "$output_file"
    fi
  done
}

# Remove bad runs from BHT2VPDMB-30_matched.txt
echo "Processing $bht2_matched..."
remove_bad_runs "$bht2_matched" "$bht2_output"

# Remove bad runs from VPDMB-30_matched.txt
echo "Processing $vpdmb_matched..."
remove_bad_runs "$vpdmb_matched" "$vpdmb_output"

echo "Done. Cleaned files:"
echo "$bht2_output"
echo "$vpdmb_output"
