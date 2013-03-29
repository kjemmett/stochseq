#!/bin/bash
#$ -S /bin/bash

######################################
time1=$( date "+%s" )
echo `pwd`
echo BEGIN `date`
echo MACHINE `hostname`
######################################

if [ $SGE_TASK_ID ]
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

observation_model=$7

if [ $8 ]
then
    # controls smp execution
    par='true'
    workers=$8
else
    par='false'
fi

cd $run_path/src

if [ $observation_model == "singlet" ]
then
/nfs/apps/matlab/current/bin/matlab -nodisplay -r "run_lambda_stochseq('$run_id', '$run_path', '$L', '$p', '$e', '$N', 'par', '$par', 'workers', '$workers'); exit;"
elif [$observation_model == "triplet" ]
then
/nfs/apps/matlab/current/bin/matlab -nodisplay -r "run_lambda_stochseq_triplet('$run_id', '$run_path', '$L', '$p', '$e', '$N', 'par', '$par', 'workers', '$workers'); exit;"
fi


########################################
echo END `date`
time2=$( date "+%s" )
echo [deltat] $(( $time2 - $time1 ))
########################################
