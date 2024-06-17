#!/usr/bin/env bash

# Ensure the script exits on any error and any unset variable
set -euo pipefail

# Ensure output directory exists
output_dir="$PWD/output"
mkdir -p "$output_dir"

# Find all ROOT files and merge them
find condorout/OutDir* -name "CALOHIST_DST_CALO_run2pp_new_2024p001-00043404-*.root" > root_files.txt

# Check if root_files.txt is not empty
if [ -s root_files.txt ]; then
  hadd -f "$output_dir/merged_file.root" $(<root_files.txt)
else
  echo "No ROOT files found to merge."
  exit 1
fi
