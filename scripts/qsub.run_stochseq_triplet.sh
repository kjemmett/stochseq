#!/bin/bash
#$ -S /bin/bash

######################################
time1=$( date "+%s" )
echo `pwd`
echo BEGIN `date`
echo MACHINE `hostname`
######################################

if [ -n $SGE_TASK_ID ]
then
    run_id=$1.$SGE_TASK_ID
else
    run_id=$1
fi

run_path=$2
L=$3
p=$4
e=$5
N=$6
if [ -n $7 ]
then
    workers=$7
fi

cd $run_path/src

/nfs/apps/matlab/current/bin/matlab -nodisplay -r "run_lambda_stochseq_triplet('$run_id', '$run_path', '$L', '$p', '$e', '$N', 'par', 'false', 'workers', '$workers'); exit;"


########################################
echo END `date`
time2=$( date "+%s" )
echo [deltat] $(( $time2 - $time1 ))
########################################
