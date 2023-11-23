#!/bin/bash

#PJM -L rscgrp=gg57,node=1
#PJM -L elapse=4:00:00
# too short
##PJM -L elapse=4:00:00
#PJM -j
#PJM -g gg57
#PJM --bulk

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 16
##$ -pe smp 28
#$ -q all.q
#$ -l hostname=a02|a03

set -eux

kind="simu"
program="ldpred"
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
    #dataname="$1"
    #phes="$1"
    #cvn="$2"
    simu_type="$1"
    phes="$2"
fi

cvi=0
phei=${job_id}
phes_ar=(${phes})
phe=${phes_ar[$phei]}

ddata="./data/${dataname}/"
dresult="./result/nonaddpgs/${dataname}/${program}/${kind}/${simu_type}/"

dtr="${dresult}/${phe}/train.cv${cvi}/"

# raw/ for flip score
dout="${dresult}/${phe}/score.cv${cvi}/raw/"
mkdir -p $dout

if [ $HPCL_HOST = "pg" ]; then
    fgenot="/grid2/ricky/data/${dataname}/ukb_imp"
else
    fgenot="${ddata}/ukb_imp"
fi

flog="./log/nonaddpgs/${dataname}/${program}/${kind}/${simu_type}/${phe}/score.cv${cvi}/score.$(date +'%y%m%d-%H%M%S').log"
mkdir -p "$(dirname $flog)"

bash ./script/lib/genetics/score.sh $dout $dtr $fgenot $threads

for f in ${dout}/*.score; do
    [ ! -f $f ] && continue
    gzip -f $f
done
