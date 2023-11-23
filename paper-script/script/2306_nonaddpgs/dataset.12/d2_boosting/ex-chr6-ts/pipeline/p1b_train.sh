#!/bin/bash

#PJM -L rscgrp=gg57,node=1
##PJM -L elapse=48:00:00
#PJM -L elapse=24:00:00
##PJM -X
#PJM -j
#PJM -g gg57
#PJM --bulk

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 48
#$ -q all.q
#$ -l hostname=a02|a03

set -eux

kind="ex-chr6-ts"
program="boosting"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh

#if [ -v SGE_TASK_ID ]; then
#    job_id=$((SGE_TASK_ID - 1))
#    #threads=$NSLOTS
#elif [ -v PJM_BULKNUM ]; then
#    job_id=$PJM_BULKNUM
#    ## FIXME
#    #threads=64
#    #threads=$PJM_MPI_PROC
#else
#    job_id=$8
#    #threads=16
#fi

if [ -v PJM_BULKNUM ]; then
    # load from env variables
    :
else
    #dataname="$1"
    #sex="$2" # male, female, both
    phes="$1"
    cvn="$2"
    #lrs="$5"
    model="$3"
    #itern="$7"
fi

#phei=$((job_id / cvn))
#cvi=$((job_id % cvn))
#phes_ar=(${phes})
#phe=${phes_ar[$phei]}

#ddata="./data/${dataname}/"
dresult="./result/nonaddpgs/${dataname}/${program}/${kind}/"
mkdir -p $dresult

for phe in $phes; do
    for cvi in $(seq 0 $((cvn - 1))); do
        dout="${dresult}/${model}/${phe}/train.cv${cvi}/"
        mkdir -p $dout

        dout_tr="./result/nonaddpgs/${dataname}/${program}/default/${model}/${phe}/train.cv${cvi}/"

        for f in ${dout_tr}/*.wgt; do
            [ ! -f $f ] && continue

            b=$(basename $f)
            fwgt="${dout}/${b}"
            # reset iteration

            if [[ $model == logitnomissing_* ]]; then
                cat $f | awk '(NR==1 || $2==CONST || $2==COV || $12!=6)' | cut -f2- | awk '{if(NR==1){print "iteration\t"$0}else{print NR-2"\t"$0}}' >$fwgt
            elif [[ $model == logitadd_* ]]; then
                cat $f | awk '(NR==1 || $2==CONST || $2==COV || $9!=6)' | cut -f2- | awk '{if(NR==1){print "iteration\t"$0}else{print NR-2"\t"$0}}' >$fwgt
            else
                exit 1
            fi

            # for score
            lr=${b/boosting_lr-/}
            lr=${lr/.wgt/}
            echo $lr
            mkdir -p ${dout}/raw/para_lr-${lr}/
            ln -sf ../../${b} ${dout}/raw/para_lr-${lr}/boosting.wgt

        done

    done
done
