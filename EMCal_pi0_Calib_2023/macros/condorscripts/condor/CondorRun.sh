#!/bin/bash

source /sphenix/u/nkumar/setup.sh

echo "------------------setting up environment--------------------"
export Cur_dir=$(pwd)
echo "running area:" ${Cur_dir}
echo "-------------------------running----------------------------"
cd ${Cur_dir}
ls
root '' > notes.log



root -b "Fun4All_G4_Waveform.C(0,\"inputdata.txt\",\"inputdatahits.txt\",\"inputdatambd.txt\",\"inputdataglobal.txt\")"
#root -b "Fun4All_G4_Waveform.C(0,\"inputdata.txt\",\"inputdatahits.txt\")"

echo "JOB COMPLETE!"
