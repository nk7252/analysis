#!/bin/bash

# Check for the correct number of arguments. for now only worry about the first
#the others are const std::string &fname = "inputdata_sp.txt", const std::string &fname_truth = "G4Hits.list"
#./NK_spmc.sh  100 inputdata_sp.txt g4hits.list
#which runs it as 
# root -l -q "Fun4All_EMCal_sp.C(100,inputdata_sp.txt, g4hits.list)"

#or try : ./NK_spmc.sh  100 dst_calo_waveform.list dst_truth.list
#root -l -q "Fun4All_EMCal_sp.C(100,dst_calo_waveform.list, dst_truth.list)"
#or try : ./NK_spmc.sh  100 dst_calo_waveform.list dst_truth.list dst_global.list

if [ "$#" -ne 4 ]; then
    echo "Usage: $0  <n_events> <inputdata_sp> <sp_truth_list> <input_global>"
    echo "Default Usage is: $0  10000 inputdata_sp.txt g4hits.list dst_global.txt"
    exit 1
fi

n_events="$1"
inputdata_sp="$2"
sp_truth_list="$3"
input_global="$4"

#create absolute path for input files
input_file1="/sphenix/user/nkumar/analysis/EMCal_pi0_Calib_2023/macros/condor/${inputdata_sp}"
input_file2="/sphenix/user/nkumar/analysis/EMCal_pi0_Calib_2023/macros/condor/${sp_truth_list}"
input_file3="/sphenix/user/nkumar/analysis/EMCal_pi0_Calib_2023/macros/condor/${input_global}"

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
root -l -q "Fun4All_EMCal_sp.C($n_events,\"$input_file1\",\"$input_file2\",\"$input_file3\")"

