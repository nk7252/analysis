#!/usr/bin/bash

export TargetDir="$PWD"/condorout


if [ -d ${TargetDir} ]; then
  if [ -n "$(ls -A ${TargetDir}/OutDir*)" ]; then
    rm -rf ${TargetDir}/OutDir*
  fi
else
  mkdir ${TargetDir}
fi

  export listfile="g4hits.list"

  #single particle pion
  #G4Hits || ! -f $listfile3 
  # Check if the list files were created successfully

  if [[ ! -f $listfile ]]; then
      echo "Error: One or more list files were not created successfully."
      exit 1
  fi

  echo "All list files created successfully."

  #number of jobs 
  j=1000

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

    mv tmp.txt ${WorkDir}/inputdatahits.txt
    
    pushd ${WorkDir}

    
    cp -v "$PWD"/../../CondorRun.sh CondorRunJob$li.sh
    cp "$PWD"/../../../Fun4All_G4_sPHENIX.C .
    cp "$PWD"/../../../G4Setup_sPHENIX.C .

    chmod +x CondorRunJob$li.sh
        
    
    cat >>ff.sub<< EOF
+JobFlavour                   = "workday"
transfer_input_files          = ${WorkDir}/CondorRunJob$li.sh, ${WorkDir}/inputdatahits.txt, ${WorkDir}/Fun4All_G4_sPHENIX.C, ${WorkDir}/G4Setup_sPHENIX.C
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

