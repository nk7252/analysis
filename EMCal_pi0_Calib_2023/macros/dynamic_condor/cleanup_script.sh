#!/bin/bash

#important setup for path if you are running on condor, -u -n
export USER="nkumar"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}
baseDir=${HOME}/analysis/EMCal_pi0_calib_2023/macros/dynamic_condor

#always run this setup anyway
source /opt/sphenix/core/bin/sphenix_setup.sh -n new
export MYINSTALL=/sphenix/u/nkumar/Install
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL 

# Run the hadd command to merge ROOT files, removed ./ from both
hadd -f "${baseDir}/output/merged_file.root" ${baseDir}/condorout/OutDir*/*.root

echo "Cleanup complete!"