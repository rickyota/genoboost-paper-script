#!/bin/bash

#PJM -L rscgrp=gg57,node=1
#PJM -L elapse=4:00:00
#PJM -j
#PJM -g gg57
#PJM --bulk

# necessary?
##PJM --omp thread=64

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 28
#$ -q all.q
## run in centos8
#$ -l hostname=a02|a03

set -eux

kind="default"
program="boosting"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh

if [ -v SGE_TASK_ID ]; then
    job_id=$((SGE_TASK_ID - 1))
    threads=$NSLOTS
else
    job_id=$4
    threads=16
    #cvi="$2"
fi

if [ -v PJM_BULKNUM ]; then
    # load from env variables
    :
else
    #dataname="$1"
    #sex="$2" # male, female, both
    phes="$1"
    cvn="$2"
    model="$3"
fi

phei=$((job_id / cvn))
cvi=$((job_id % cvn))

phes_ar=(${phes})
phe=${phes_ar[$phei]}

ddata="./data/${dataname}/"
dresult="./result/nonaddpgs/${dataname}/${program}/${kind}/"
mkdir -p $dresult

dout="${dresult}/${model}/${phe}/score.cv${cvi}/"
mkdir -p $dout
dtr="${dresult}/${model}/${phe}/train.cv${cvi}/"

if [ $HPCL_HOST = "pg" ]; then
    fgenot="/grid2/ricky/data/${dataname}/ukb_imp"
else
    fgenot="${ddata}/ukb_imp"
fi
fphe="${ddata}ukb_imp.phe"
#fcv="${ddata}ukb_imp.cv"

flog="./log/nonaddpgs/${dataname}/${program}/${kind}/${model}/${phe}/score.cv${cvi}/score.$(date +'%y%m%d-%H%M%S').log"
mkdir -p "$(dirname $flog)"

bash ./script/lib/boosting/score.sh $dout $dtr $fgenot $fphe $threads |& tee $flog

# boosting.score and boosting.scorecov
gzip -f ${dout}/${program}.score
gzip -f ${dout}/${program}.scorecov

#for f in ${dout}/*.score; do
#    [ ! -f $f ] && continue
#    gzip -f $f
#done
#
#for f in ${dout}/*.scorecov; do
#    [ ! -f $f ] && continue
#    gzip -f $f
#done
