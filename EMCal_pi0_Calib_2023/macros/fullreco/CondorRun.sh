#!/bin/bash

source /sphenix/u/nkumar/setup.sh

echo "------------------setting up environment--------------------"
export Cur_dir=$(pwd)
echo "running area:" ${Cur_dir}
echo "-------------------------running----------------------------"
cd ${Cur_dir}
ls
root '' > notes.log



root -b "Fun4All_G4_sPHENIX.C(0,\"inputdatahits.txt\")"
#root -b "Fun4All_EMCal_sp.C(0,\"inputdata.txt\",\"inputdatahits.txt\")"

echo "JOB COMPLETE!"