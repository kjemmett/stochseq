#!/bin/bash
#$ -S /bin/bash

######################################
time1=$( date "+%s" )
echo `pwd`
echo BEGIN `date`
echo MACHINE `hostname`
######################################

run_path=$1
cd $run_path

# pull parameters from param_list
foo=`sed -n "$SGE_TASK_ID"p param_list`
IFS=' ' read -a bar <<< $foo
L=${bar[0]}
p=${bar[1]}
e=${bar[2]}
N=${bar[3]}
T=${bar[4]}

run_id=L$L-p$p-e$e-N$N
cd $run_path/src
/nfs/apps/matlab/2012a/bin/matlab -nodisplay -r "run_random_stochseq('$run_id', '$run_path', '$L', '$p', '$e', '$N', '$T', '$SGE_TASK_ID'); exit;"

########################################
echo END `date`
time2=$( date "+%s" )
echo [deltat] $(( $time2 - $time1 ))
########################################
