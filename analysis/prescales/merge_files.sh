#!/bin/bash

# Input file prefixes
prefix1="BHT2VPDMB"
prefix2="VPDMB"

# Output file names for both sets
matched_output1="matched_output_BHT2VPDMB.txt"
unmatched_output1="unmatched_output_BHT2VPDMB.txt"
matched_output2="matched_output_VPDMB.txt"
unmatched_output2="unmatched_output_VPDMB.txt"

# Temporary files to store sorted data
sorted_file1="sorted_file1.txt"
sorted_file2="sorted_file2.txt"
temp_file="temp_merged.txt"

# Counters for both sets of files
matched_counter1=0
unmatched_file1_counter1=0
unmatched_file2_counter1=0

matched_counter2=0
unmatched_file1_counter2=0
unmatched_file2_counter2=0

# Function to process files
process_files() {
  local file1=$1
  local file2=$2
  local matched_output=$3
  local unmatched_output=$4
  local matched_counter=$5
  local unmatched_file1_counter=$6
  local unmatched_file2_counter=$7

  # Check if both files exist
  if [[ ! -f "$file1" || ! -f "$file2" ]]; then
    echo "One or both input files do not exist for $file1 and $file2."
    exit 1
  fi

  # Sort both files based on the appropriate columns
  echo "Sorting input files..."
  sort -k 1,1 "$file1" -o "$sorted_file1"      # Sort File 1 by first column (start time)
  sort -k 2,2 "$file2" -o "$sorted_file2"      # Sort File 2 by second column (start time)

  # Join files on start_time (second column in file2, first column in file1)
  # -a 1: Include all lines from file1, even if they don't match
  # -a 2: Include all lines from file2, even if they don't match
  echo "Joining sorted files..."
  join -1 1 -2 2 -a 1 -a 2 -o 2.1 1.2 2.5 2.6 2.7 "$sorted_file1" "$sorted_file2" > "$temp_file"

  # Initialize output files
  echo "Run_Number Number_of_Events Sampled_Luminosity Prescale Livetime" > "$matched_output"
  echo "Run_Number Number_of_Events Sampled_Luminosity Prescale Livetime" > "$unmatched_output"

  # Process merged data with better handling for missing columns
  echo "Processing merged data..."
  while IFS= read -r line; do
    # Split the line into fields using `awk`
    run_number=$(echo "$line" | awk '{print $1}')
    number_of_events=$(echo "$line" | awk '{print $2}')
    sampled_luminosity=$(echo "$line" | awk '{print $3}')
    prescale=$(echo "$line" | awk '{print $4}')
    livetime=$(echo "$line" | awk '{print $5}')

    # Check for missing fields
    if [[ -z "$run_number" || -z "$number_of_events" || -z "$sampled_luminosity" || -z "$prescale" || -z "$livetime" ]]; then
      # If any field is missing, write to the unmatched output
      echo "$line" >> "$unmatched_output"

      # Determine which file the missing data comes from
      if [[ -z "$number_of_events" ]]; then
        ((unmatched_file1_counter++))  # File 1 is missing data
      else
        ((unmatched_file2_counter++))  # File 2 is missing data
      fi
    else
      # If all fields are present, write to the matched output
      echo "$line" >> "$matched_output"
      ((matched_counter++))
    fi
  done < "$temp_file"

  # Output the counters
  echo "Matched lines: $matched_counter"
  echo "Unmatched lines from File 1: $unmatched_file1_counter"
  echo "Unmatched lines from File 2: $unmatched_file2_counter"
}

# Process the first set (BHT2VPDMB)
file1_set1="BHT2VPDMB-30_2_subtr.txt"
file2_set1="BHT2VPDMB-30_9.txt"
process_files "$file1_set1" "$file2_set1" "$matched_output1" "$unmatched_output1" \
              matched_counter1 unmatched_file1_counter1 unmatched_file2_counter1

# Process the second set (VPDMB)
file1_set2="VPDMB-30_2_subtr.txt"
file2_set2="VPDMB-30_9.txt"
process_files "$file1_set2" "$file2_set2" "$matched_output2" "$unmatched_output2" \
              matched_counter2 unmatched_file1_counter2 unmatched_file2_counter2

# Cleanup temporary files
echo "Cleaning up..."
rm "$temp_file" "$sorted_file1" "$sorted_file2"

echo "Done."
