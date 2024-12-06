#!/usr/bin/env bash

# Ensure the script exits on any error and any unset variable
set -euo pipefail

# Ensure output directory exists
output_dir="$PWD/output"
mkdir -p "$output_dir"

# Find all ROOT files and list them
find condorout/OutDir* -name "OUTHIST*.root" > root_files.txt

# File to keep track of the global counter
counter_file="$output_dir/run_counter.txt"

# Initialize counter if the file does not exist
if [ ! -f "$counter_file" ]; then
  echo 0 > "$counter_file"
fi

# Read the current counter value
global_counter=$(cat "$counter_file")

# Increment the counter
global_counter=$((global_counter + 1))

# Save the updated counter value
echo "$global_counter" > "$counter_file"


# Check if root_files.txt is not empty
if [ -s root_files.txt ]; then
  # Extract the common part of the file names, ignoring directory names
  common_part=$(awk '
  {
    # Extract the filename from the full path
    n = split($0, arr, "/")
    filename = arr[n]
    
    # Initialize the common part or compare with the next filename
    if (NR == 1) {
      common = filename
    } else {
      i = 1
      while (substr(common, i, 1) == substr(filename, i, 1)) {
        i++
      }
      common = substr(common, 1, i - 1)
    }
  }
  END {
    print common
  }' root_files.txt)

  # Ensure the common part is not empty
  if [ -z "$common_part" ]; then
    echo "Failed to determine the common part of the ROOT file names."
    exit 1
  fi

  # Add the word "merged" to the common part
  common_part="${common_part}merged_V${global_counter}"

  # Clean up the common part to make it suitable for a filename
  common_part=$(echo "$common_part" | sed 's/[^a-zA-Z0-9]/_/g')

  # Set the output file name
  output_file="${output_dir}/${common_part}.root"

  # Merge the ROOT files
  hadd -f  -j 8 "$output_file" $(<root_files.txt)
else
  echo "No ROOT files found to merge."
  exit 1
fi

echo "Merged ROOT file: $output_file"
echo "Done."
