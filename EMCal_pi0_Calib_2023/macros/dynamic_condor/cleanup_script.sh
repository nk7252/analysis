#!/bin/bash

#important setup for path if you are running on condor, -u -n
export USER="nkumar -u -n"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}
baseDir=${HOME}/analysis/EMCal_pi0_calib_2023/macros/dynamic_condor

# Run the hadd command to merge ROOT files, removed ./ from both
hadd -f "${baseDir}/output/merged_file.root" ${baseDir}/condorout/OutDir*/*.root

echo "Cleanup complete!"