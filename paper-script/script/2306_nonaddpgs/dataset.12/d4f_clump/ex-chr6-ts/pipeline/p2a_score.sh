#!/bin/bash

#PJM -L rscgrp=gg57,node=1
#PJM -L elapse=3:00:00
#PJM -j
#PJM -g gg57
#PJM --bulk

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 16
#$ -q all.q
#$ -l hostname=a02|a03

set -eux

kind="ex-chr6-ts"
program="clump"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh

if [ -v SGE_TASK_ID ]; then
    job_id=$((SGE_TASK_ID - 1))
    threads=$NSLOTS
elif [ -v PJM_BULKNUM ]; then
    job_id=$PJM_BULKNUM
    ## FIXME
    threads=64
else
    job_id=$3
    threads=16
fi

if [ -v PJM_BULKNUM ]; then
    # load from env variables
    :
else
    phes="$1"
    cvn="$2"
fi

cvi=$((job_id % cvn))
phei=$((job_id / cvn))
phes_ar=(${phes})
phe=${phes_ar[$phei]}

flog="./log/nonaddpgs/${dataname}/${program}/${kind}/${phe}/score.cv${cvi}/score.$(date +'%y%m%d-%H%M%S').log"
mkdir -p "$(dirname $flog)"

ddata="./data/${dataname}/"
dresult="./result/nonaddpgs/${dataname}/${program}/${kind}/"
mkdir -p $dresult

dwgt="${dresult}/${phe}/train.cv${cvi}/"

dout="${dresult}/${phe}/score.cv${cvi}/"
mkdir -p $dout

if [ $HPCL_HOST = "pg" ]; then
    fgenot="/grid2/ricky/data/${dataname}/ukb_imp"
else
    fgenot="${ddata}/ukb_imp"
fi

bash ./script/lib/genetics/score_no-concat_use-snv-pos.sh \
    $dout $dwgt $fgenot $threads "n" |& tee $flog

bash ./script/lib/genetics/score_concat_use-snv-pos.sh \
    $dout $dwgt $fgenot $threads "n" |& tee $flog

for f in ${dout}/*.score; do
    [ ! -f $f ] && continue
    gzip -f $f
done
