#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 28
#$ -q all.q
#$ -l hostname=a02|a03|z01|z02|g01|g02

set -eux

program="sbayesr"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

if [ -v SGE_TASK_ID ]; then
    job_id=$((SGE_TASK_ID - 1))
    #threads=$NSLOTS
else
    job_id=$3
    #threads=16
fi

phes="$1"
cvn="$2"

phei=$((job_id / cvn))
cvi=$((job_id % cvn))

phes_ar=(${phes})
phe=${phes_ar[$phei]}

ddata="./data/${dataname}/"
dtmp="${ddata}/tmp/${program}/"
dresult="./result/nonaddpgs/${dataname}/ldpred/"
mkdir -p $dresult

dout="${dresult}/${phe}/train.cv${cvi}/raw/"
mkdir -p $dout

fma="${dtmp}${phe}.cv${cvi}.ma"

fldm="./data/sbayesr/ldm_sparse/ukbEURu_imp_v3_HM3_n50k.chisq10.ldm.sparse"

bash ./script/lib/${program}/train.sh $dout $fma $fldm
