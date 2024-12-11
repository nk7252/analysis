#!/usr/bin/bash

#important, setup for path if you are running on condor, -u -n removed
export USER="nkumar"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}
baseDir=${HOME}/analysis/EMCal_pi0_Calib_2023/macros
echo "Base Directory: $baseDir"
# Setting the target directory
export TargetDir=${baseDir}/dynamic_condor/condorout
echo "Target Directory: $TargetDir"
# Cleaning or creating the target directory
if [ -d ${TargetDir} ]; then
  echo "Target directory exists."
  if [ -n "$(ls -A ${TargetDir}/OutDir*)" ]; then
    echo "Cleaning target directory."
    rm -rf ${TargetDir}/OutDir*
  else
    echo "Target directory is empty."
  fi
else
  echo "Creating target directory."
  mkdir -p ${TargetDir}
fi

# Check if directory was created successfully
if [ ! -d ${TargetDir} ]; then
  echo "Failed to create target directory: ${TargetDir}"
  exit 1
fi


# Defining list files and cleaning old ones
export listfile="dst_calo_cluster.list"
export listfile2="g4hits.list"
rm -f $listfile
rm -f $listfile2

# Creating new file lists
CreateFileList.pl -type 14  -run 13 -particle pi0 -pmin 200 -pmax 10000 DST_CALO_CLUSTER G4Hits

# Calculating the number of jobs
j=500
# Count the total number of files in the listfile
tot_files=$(cat ${listfile} | wc -l)
echo "total files: $tot_files"
# Calculate the remainder when total files are divided by the number of jobs
rem=$((tot_files % j))
# Calculate the number of files per job by dividing the total files by the number of jobs
files_per_job=$((tot_files / j))
# Initially set the number of jobs to the desired value
njob=$j
# If there is a remainder, increment files per job by one to account for the remaining files
if [ $rem -ne 0 ]; then
  files_per_job=$((files_per_job + 1))
fi
# Recalculate the remainder when total files are divided by the updated files per job
rem2=$((tot_files % files_per_job))
# Calculate the number of jobs by dividing the total files by files per job
njob=$((tot_files / files_per_job))
# If there is a remainder, increment the number of jobs by one to account for the remaining files
if [ $rem2 -ne 0 ]; then
  njob=$((tot_files / files_per_job + 1))
fi
# Output the calculated number of files per job and the number of jobs
echo "files per job: $files_per_job"
echo "njob: $njob"

# Initializing the DAG file
dagfile="jobs.dag"
rm -f $dagfile

# Verify the source files exist before entering the loop
if [ ! -f ${baseDir}/dynamic_condor/CondorRun.sh ]; then
  echo "CondorRun.sh does not exist at ${baseDir}/dynamic_condor/CondorRun.sh"
  exit 1
fi

if [ ! -f ${baseDir}/Fun4All_EMCal_sp.C ]; then
  echo "Fun4All_EMCal_sp.C does not exist at ${baseDir}/Fun4All_EMCal_sp.C"
  exit 1
fi

# Loop to create subjobs
for ((q = 0; q < njob; q++)); do
  #create output dir for each job
  mkdir -p ${TargetDir}/OutDir$q
  export WorkDir="${TargetDir}/OutDir$q"
  # Debugging: Check if WorkDir was created successfully
  if [ ! -d ${WorkDir} ]; then
    echo "Failed to create working directory: ${WorkDir}"
    exit 1
  fi
  echo "WorkDir: ${WorkDir}"
  #calculate start and end files(of total list) for each job
  start_file=$((q * files_per_job + 1))
  end_file=$((start_file + files_per_job - 1))
  echo "start file: $start_file   end file: $end_file"
  #extract relevant filenames from master listfile and makes a new listfile for the individual job in the job workdir
  sed -n $start_file,${end_file}p ${listfile} > ${WorkDir}/inputdata.txt
  sed -n $start_file,${end_file}p ${listfile2} > ${WorkDir}/inputdatahits.txt
  # Debugging: Ensure inputdata.txt and inputdatahits.txt are created
  if [ ! -f ${WorkDir}/inputdata.txt ]; then
    echo "Failed to create inputdata.txt in ${WorkDir}"
    exit 1
  fi

  if [ ! -f ${WorkDir}/inputdatahits.txt ]; then
    echo "Failed to create inputdatahits.txt in ${WorkDir}"
    exit 1
  fi
  #copy job template and macro in to job workdir and rename job template
  cp -v ${baseDir}/dynamic_condor/CondorRun.sh ${WorkDir}/CondorRunJob$q.sh
  cp ${baseDir}/Fun4All_EMCal_sp.C ${WorkDir}
  # Debugging: Ensure files are copied successfully
  if [ ! -f ${WorkDir}/CondorRunJob$q.sh ]; then
    echo "Failed to copy CondorRun.sh to ${WorkDir}/CondorRunJob$q.sh"
    exit 1
  fi

  if [ ! -f ${WorkDir}/Fun4All_EMCal_sp.C ]; then
    echo "Failed to copy Fun4All_EMCal_sp.C to ${WorkDir}"
    exit 1
  fi

  chmod +x ${WorkDir}/CondorRunJob$q.sh
  #create submission file for each job
  cat > ${WorkDir}/ff.sub << EOF
+JobFlavour                   = "workday"
transfer_input_files          = ${WorkDir}/CondorRunJob$q.sh, ${WorkDir}/inputdata.txt, ${WorkDir}/Fun4All_EMCal_sp.C, ${WorkDir}/inputdatahits.txt
Executable                    = CondorRunJob$q.sh
request_memory                = 2GB
Universe                      = vanilla
Notification                  = Never
GetEnv                        = True
Priority                      = +12
Output                        = condor.out
Error                         = condor.err
Log                           = /tmp/condor$q.log
Queue
EOF
 #should_transfer_files         = YES
 #when_to_transfer_output       = ON_EXIT
  #transfer_output_files         = condor.out, condor.err, condor.log
  # each job submission file is added to the dag file. with a jobname (Job$q) and path to the submission file
  echo "JOB Job$q ${WorkDir}/ff.sub" >> $dagfile
done

# Add the final cleanup job
cat > cleanup.sub << EOF
universe        = vanilla
executable      = ${baseDir}/dynamic_condor/cleanup_script.sh
arguments       =
output          = cleanup_output.txt
error           = cleanup_error.txt
log             = cleanup_log.txt
queue
EOF
#transfer_output_files = ./output/merged_file.root
#should_transfer_files = YES
#when_to_transfer_output = ON_EXIT
#add a cleanup job to the end of the dagfile
echo "JOB CLEANUP cleanup.sub" >> $dagfile
#mark all jobs from 0->njob-1 as parents, and the cleanup as child. child jobs run only after ALL parent jobs are completed.
echo "PARENT $(seq -s ' ' -f "Job%.0f" 0 $((njob - 1))) CHILD CLEANUP" >> $dagfile

# Submit the DAG file
condor_submit_dag -f -debug 3 -maxjobs 50 $dagfile