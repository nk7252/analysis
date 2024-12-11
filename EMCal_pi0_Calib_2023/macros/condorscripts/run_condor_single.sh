#!/bin/bash

export macropath=$(pwd)

# Path to the dst files
export dstpath="/sphenix/user/shuhangli/FunWithML/TreeMaker/macro/condorphotonlist"
# Path to the output directory
export TargetDir="/sphenix/user/shuhangli/FunWithML/TreeMaker/macro/condorphotonlist"
# Total events for this condor. (1 line = 1000 events)
total_line=30000
#total_line=10031
# Events per job (1 line = 1000 events)
lines_per_job=1
# Skip lines (for the case of resubmission or continuation)
skip_lines=0

total_dst_lines=$(wc -l < $dstpath/dst_calo_cluster.list)
if [ $total_line -gt $total_dst_lines ]; then
    echo "Error: total_lines is greater than the number of lines in the dst files."
    echo "Total DST lines: $total_dst_lines"
    exit 1
fi
total_jobs=$(( (total_line + lines_per_job - 1) / lines_per_job ))

for ((i=0; i<total_jobs; i++)); do
    start_line=$(( i * lines_per_job + 1 + skip_lines))
    end_line=$(( start_line + lines_per_job - 1 ))

    mkdir -p "${TargetDir}/OutDir$i"
    export WorkDir="${TargetDir}/OutDir$i"
    pushd "${WorkDir}"

    cp -v "$macropath/CondorRunSim.sh" "CondorRunTC$i.sh"
    cp "$macropath/Fun4All_run_dst.C" "${WorkDir}/"
    chmod +x "CondorRunTC$i.sh"

    sed -n "${start_line},${end_line}p" $dstpath/dst_calo_cluster.list > dst_calo_cluster.list
    sed -n "${start_line},${end_line}p" $dstpath/g4hits.list > g4hits.list
    #sed -n "${start_line},${end_line}p" $dstpath/dst_mbd_epd.list > dst_mbd_epd.list


    cat >>ff.sub<< EOF
+JobFlavour                   = "espresso"
transfer_input_files          = ${WorkDir}/CondorRunTC$i.sh,${WorkDir}/Fun4All_run_dst.C,${WorkDir}/dst_calo_cluster.list,${WorkDir}/g4hits.list
Executable                    = CondorRunTC$i.sh
request_memory                = 10GB
Universe                      = vanilla
Notification                  = Never
GetEnv                        = True
Priority                      = +80
Output                        = test.out
Error                         = test.err
Log                           = /tmp/sli_$i.log
Notify_user                   = sl4859@columbia.edu

Queue
EOF

    condor_submit ff.sub
    popd
done
