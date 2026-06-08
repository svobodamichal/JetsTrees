#!/bin/bash

# Extract run numbers from the first column and sort them
cut -d' ' -f1 VPDMB-30_matched_cleaned.txt | sort -n > vpdmb_runs_sorted.txt
cut -d' ' -f1 BHT2VPDMB-30_matched_cleaned.txt | sort -n > bht2vpdmb_runs_sorted.txt

# Use comm to find differences between the two sorted files
# comm outputs:
#   - Column 1: Lines unique to vpdmb_runs_sorted.txt
#   - Column 2: Lines unique to bht2vpdmb_runs_sorted.txt
#   - Column 3: Lines common to both files

echo "Run numbers in VPDMB but not in BHT2VPDMB:" > missing_runs.txt
comm -23 vpdmb_runs_sorted.txt bht2vpdmb_runs_sorted.txt >> missing_runs.txt

echo "Run numbers in BHT2VPDMB but not in VPDMB:" >> missing_runs.txt
comm -13 vpdmb_runs_sorted.txt bht2vpdmb_runs_sorted.txt >> missing_runs.txt

echo "Comparison complete. Check 'missing_runs.txt' for results."

# Cleanup temporary files
rm vpdmb_runs_sorted.txt bht2vpdmb_runs_sorted.txt

