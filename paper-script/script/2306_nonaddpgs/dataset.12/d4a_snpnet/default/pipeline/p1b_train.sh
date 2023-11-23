#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 28
#$ -q all.q
#$ -l hostname=a02|a03|z01|z02|g01|g02

## Extremely slow at a02, a03  or simultaneously running at the same node

set -eux

kind="default"
program="snpnet"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

if [ -v SGE_TASK_ID ]; then
    job_id=$((SGE_TASK_ID - 1))
    threads=$NSLOTS
else
    job_id=$4
    threads=32
fi

sex="$1" # male, female, both
phes="$2"
cvn="$3"

phei=$((job_id / cvn))
cvi=$((job_id % cvn))

phes_ar=(${phes})
phe=${phes_ar[$phei]}

flog="./log/nonaddpgs/${dataname}/${program}/${kind}/${phe}/train.cv${cvi}/train.$(date +'%y%m%d-%H%M%S').log"
mkdir -p "$(dirname $flog)"

ddata="./data/${dataname}/"
dtmp="${ddata}/tmp/${program}/"
dresult="./result/nonaddpgs/${dataname}/${program}/${kind}/"
mkdir -p $dresult

dout="${dresult}/${phe}/train.cv${cvi}/raw/"
mkdir -p $dout

#if [ $HPCL_HOST = "pg" ]; then
#    :
#    # ny
#    #fgenot="/grid2/ricky/data/${dataname}/tmp/${program}/ukb_imp"
#else
#    fgenot="${dtmp}/ukb_imp"
#fi

fgenot="${dtmp}ukb_imp"
fphe="${dtmp}ukb_imp.phe"

if [ $sex == "both" ]; then
    fsamples="${dtmp}/ukb_imp.cv${cvi}.samples"
elif [ $sex == "female" ]; then
    fsamples="${dtmp}/ukb_imp.cv${cvi}.female.samples"
fi

bash ./script/lib/${program}/train.sh \
    $phe $dout $fgenot $fphe $fsamples $sex $threads |& tee $flog
