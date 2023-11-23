#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 8
#$ -q all.q
#$ -l hostname=a02|a03

set -eux

program="snpnet"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh

if [ -v SGE_TASK_ID ]; then
    job_id=$((SGE_TASK_ID - 1))
    threads=$NSLOTS
else
    job_id=$4
    threads=16
fi

phes="$1"
cvn="$2"

cvi=$((job_id % cvn))
phei=$((job_id / cvn))
phes_ar=(${phes})
phe=${phes_ar[$phei]}

ddata="./data/${dataname}/"
dresult="./result/nonaddpgs/${dataname}/${program}/"
mkdir -p $dresult

if [ $HPCL_HOST = "pg" ]; then
    fgenot="/grid2/ricky/data/${dataname}/ukb_imp"
else
    fgenot="${ddata}/ukb_imp"
fi

dwgt="${dresult}/${phe}/train.cv${cvi}/"

dout="${dresult}/${phe}/score.cv${cvi}/raw/"
mkdir -p $dout

if [ $HPCL_HOST = "pg" ]; then
    fgenot="/grid2/ricky/data/${dataname}/ukb_imp"
else
    fgenot="${ddata}/ukb_imp"
fi
fphe="${ddata}/ukb_imp.phe"

bash ./script/lib/genetics/score.sh \
    $dout $dwgt $fgenot $fphe $phe $threads
