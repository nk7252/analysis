#!/bin/bash

#./test_run.sh 100 /sphenix/user/nkumar/analysis/EMCal_pi0_Calib_2023/macros/listfiles/single/run24/pi0/dst_calo_cluster.list /sphenix/user/nkumar/analysis/EMCal_pi0_Calib_2023/macros/listfiles/single/run24/pi0/g4hits.list
#which runs it as 
# root -l -q "Fun4All_EMCal_sp.C(100, /sphenix/user/nkumar/analysis/EMCal_pi0_Calib_2023/macros/listfiles/single/run24/pi0/dst_calo_cluster.list, /sphenix/user/nkumar/analysis/EMCal_pi0_Calib_2023/macros/listfiles/single/run24/pi0/g4hits.list)"


if [ "$#" -ne 3 ]; then
    echo "Usage: $0  <n_events> <cluster_list> <truth_list> <mbd_list> <global_list> <truth_jet_list>"
    echo "Default Usage is: $0  100 dst_calo_cluster.list dst_truth.list dst_mbd_epd.list dst_global.list dst_truth_jet.list"
    exit 1
fi

n_events="$1"
cluster_list="$2"
truth_list="$3"

#create absolute path for input files
#input_file1="/sphenix/user/nkumar/analysis/EMCal_pi0_Calib_2023/macros/condorscripts/condor/${cluster_list}"
#input_file2="/sphenix/user/nkumar/analysis/EMCal_pi0_Calib_2023/macros/condorscripts/condor/${truth_list}"
#export listfile="/sphenix/user/nkumar/analysis/EMCal_pi0_Calib_2023/macros/listfiles/single/run24/pi0/dst_calo_cluster.list"
#export listfile2="/sphenix/user/nkumar/analysis/EMCal_pi0_Calib_2023/macros/listfiles/single/run24/pi0/g4hits.list"



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
root -l -q "Fun4All_G4_Waveform.C($n_events,\"$cluster_list\",\"$truth_list\")"
