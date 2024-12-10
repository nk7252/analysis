#!/bin/bash

#./test_run.sh 100 dst_calo_cluster.list
#which runs it as 
# root -l -q "Fun4All_G4_Waveform.C(100, dst_calo_cluster.list)"


if [ "$#" -ne 2 ]; then
    echo "Usage: $0  <n_events> <cluster_list> "
    echo "Default Usage is: $0  100 dst_calo_cluster.list"
    exit 1
fi

n_events="$1"
cluster_list="$2"

#create absolute path for input files
input_file1="/sphenix/user/nkumar/analysis/EMCal_pi0_Calib_2023/macros/condor/${cluster_list}"

#important setup for path if you are running on condor
export USER="nkumar -u -n"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}
baseDir=/sphenix/user/nkumar/analysis/EMCal_pi0_Calib_2023/macros

#always run this setup anyway
source /opt/sphenix/core/bin/sphenix_setup.sh -n new
export MYINSTALL=/sphenix/u/nkumar/Install
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL 

#check if everything is set correctly
echo $MYINSTALL 
echo $baseDir

# cd in to basedir to run script. print environment for debugging
printenv
cd $baseDir

# Execute the custom ROOT script with arguments
root -l -q "Fun4All_G4_Waveform.C($n_events,\"$input_file1\")"

