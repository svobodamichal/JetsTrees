#!/bin/bash

# Input file
input_file1="BHT2VPDMB-30_2.txt"
output_file1="BHT2VPDMB-30_2_subtr.txt"

input_file2="VPDMB-30_2.txt"
output_file2="VPDMB-30_2_subtr.txt"

# Use awk to process the file
awk 'NR == 1 {print $1, int($2); prev=int($2); next} {print $1, $2 - prev; prev = $2}' "$input_file1" > "$output_file1"
awk 'NR == 1 {print $1, int($2); prev=int($2); next} {print $1, $2 - prev; prev = $2}' "$input_file2" > "$output_file2"

echo "Processed output saved to $output_file1 and $output_file2"
