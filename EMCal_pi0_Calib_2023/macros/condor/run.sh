#!/usr/bin/bash

export TargetDir="$PWD"/condorout


if [ -d ${TargetDir} ]; then
  if [ -n "$(ls -A ${TargetDir}/OutDir*)" ]; then
    rm -rf ${TargetDir}/OutDir*
  fi
else
  mkdir ${TargetDir}
fi

  export listfile="dst_calo_waveform.list"
  #export listfile="dst_calo_cluster.list"
  #export listfile2="g4hits.list"
  export listfile2="dst_truth.list"
  export listfile3="dst_global.list"
  #export listfile3="dst_truth_g4hit.list"
  #export listfile3="g4hits.list"

  #only delete and regenerate if switching file source or number of events
  #rm $listfile
  #rm $listfile2
  #rm $listfile3

  #single particle pion
  #CreateFileList.pl -type 14  -run 13 -particle pi0 -pmin 200 -pmax 10000 DST_CALO_CLUSTER G4Hits 
  #pythia pp. run 11 is pythia pp min bias. run 15 is pythia pp 20 micro-s streaming
  #CreateFileList.pl DST_CALO_CLUSTER G4Hits -type 3 -run 11 -nopileup

  #CreateFileList.pl DST_CALO_WAVEFORM DST_TRUTH -type 3 -run 15 -nopileup -n 1000000000
  
  #
  #CreateFileList.pl -run 15 -type 3 -nop DST_CALO_CLUSTER DST_TRUTH -n 100000
  #G4Hits || ! -f $listfile3 
  #run 111 also works see https://wiki.sphenix.bnl.gov/index.php?title=MDC2_2022
# to test use a small set. like  -n 1000
  # Check if the list files were created successfully
  if [[ ! -f $listfile || ! -f $listfile2 || ! -f $listfile3 ]]; then
      echo "Error: One or more list files were not created successfully."
      exit 1
  fi

  echo "All list files created successfully."


  #DST_GLOBAL-nopileup-n 10000000 DST_CALO_CLUSTER
  #DST_CALO_WAVEFORM
  #number of jobs 
  j=5000


  # Count the number of lines in dst_calo_cluster.list
  #num_lines=$(wc -l < dst_calo_cluster.list)

  # Set j to half the number of lines (round down)
  #j=$((num_lines / 2))

  # Cap j at 1000 if it exceeds this value
  #if [ $j -gt 1000 ]; then
  #  j=1000
  #fi

  tot_files=$( cat ${listfile} | wc -l )
  echo "total files: $tot_files"
  rem=$(( $tot_files%$j ))
  files_per_job=$(( $tot_files/$j ))
  njob=$j
  if [ $rem -ne 0 ]; then
    files_per_job=$(( $files_per_job+1 ))
  fi
  rem2=$(( $tot_files%$files_per_job ))
  njob=$(( $tot_files/$files_per_job ))
  if [ $rem2 -ne 0 ]; then
    njob=$(( ($tot_files/$files_per_job)+1 ))
  fi
  echo "files per job: $files_per_job"
  echo "njob: $njob"


  for((q=0;q<$njob;q++));
  do

    mkdir ${TargetDir}/OutDir$q
    export WorkDir="${TargetDir}/OutDir$q"
    echo "WorkDir:" ${WorkDir}

    start_file=$(( $q*$files_per_job+1 ))
    end_file=$(( $start_file+$files_per_job-1 ))
    echo "start file: $start_file   end file: $end_file"

    sed -n $start_file\,${end_file}p ${listfile} > tmp.txt
    sed -n $start_file\,${end_file}p ${listfile2} > tmp2.txt
    sed -n $start_file\,${end_file}p ${listfile3} > tmp3.txt
    mv tmp.txt ${WorkDir}/inputdata.txt
    mv tmp2.txt ${WorkDir}/inputdatahits.txt
    mv tmp3.txt ${WorkDir}/inputdataglobal.txt
    
    pushd ${WorkDir}

    
    cp -v "$PWD"/../../CondorRun.sh CondorRunJob$li.sh
    cp "$PWD"/../../../Fun4All_EMCal_sp.C .

    chmod +x CondorRunJob$li.sh
        
    
    cat >>ff.sub<< EOF
+JobFlavour                   = "workday"
transfer_input_files          = ${WorkDir}/CondorRunJob$li.sh , ${WorkDir}/Fun4All_EMCal_sp.C , ${WorkDir}/inputdata.txt , ${WorkDir}/inputdatahits.txt , ${WorkDir}/inputdataglobal.txt
Executable                    = CondorRunJob$li.sh
request_memory                = 10GB
Universe                      = vanilla
Notification                  = Never
GetEnv                        = True
Priority                      = +12
Output                        = condor.out
Error                         = condor.err
Log                           = /tmp/condor$li.log
PeriodicHold                  = (NumJobStarts>=1 && JobStatus == 1)
concurrency_limits            = CONCURRENCY_LIMIT_DEFAULT:100


Queue
EOF

    condor_submit ff.sub
    popd
  
done

