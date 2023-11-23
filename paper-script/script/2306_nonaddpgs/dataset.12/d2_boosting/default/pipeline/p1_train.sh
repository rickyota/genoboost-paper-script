#!/bin/bash

#PJM -L rscgrp=gg57,node=1
##PJM -L elapse=48:00:00
#PJM -L elapse=24:00:00
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
#$ -l hostname=a02|a03
## run in centos8
##$ -l hostname=a02|a03|z01|z02|g01|g02

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
elif [ -v PJM_BULKNUM ]; then
    job_id=$PJM_BULKNUM
    ## FIXME
    threads=64
    #threads=$PJM_MPI_PROC
else
    job_id=$6
    threads=16
fi

if [ -v PJM_BULKNUM ]; then
    # load from env variables
    :
else
    #dataname="$1"
    sex="$1" # male, female, both
    phes="$2"
    cvn="$3"
    model="$4"
    lrs="$5"
    #itern="$6"
fi

phei=$((job_id / cvn))
cvi=$((job_id % cvn))

phes_ar=(${phes})
phe=${phes_ar[$phei]}

ddata="./data/${dataname}/"
dresult="./result/nonaddpgs/${dataname}/${program}/${kind}/"
mkdir -p $dresult

dout="${dresult}/${model}/${phe}/train.cv${cvi}/"
mkdir -p $dout

if [ $HPCL_HOST = "pg" ]; then
    fgenot="/grid2/ricky/data/${dataname}/ukb_imp"
else
    fgenot="${ddata}/ukb_imp"
fi
fphe="${ddata}ukb_imp.phe"
fcv="${ddata}ukb_imp.cv"

flog="./log/nonaddpgs/${dataname}/${program}/${kind}/${model}/${phe}/train.cv${cvi}/train.$(date +'%y%m%d-%H%M%S').log"
mkdir -p "$(dirname $flog)"

#bash ./script/lib/boosting/freemodel_230507.sh \
bash ./script/lib/boosting/freemodel.sh \
    $phe $cvi "$lrs" $model $sex $dout $fgenot $fphe $fcv $threads |&
    tee $flog

# $phe $cvi "$lrs" $itern $model $sex $dout $fgenot $fphe $fcv $threads |&
