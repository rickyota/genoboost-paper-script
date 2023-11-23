#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 4
#$ -q all.q
#$ -l hostname=a02|a03|z01|z02|g01|g02

set -eux

program="covonly"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh

if [ -v SGE_TASK_ID ]; then
    job_id=$((SGE_TASK_ID - 1))
    #threads=$NSLOTS
else
    job_id=$4
    #job_id=$5
    #threads=16
fi

sex="$1" # male, female, both
phes="$2"
cvn="$3"

cvi=$((job_id % cvn))
phei=$((job_id / cvn))
phes_ar=(${phes})
phe=${phes_ar[$phei]}

ddata="./data/${dataname}/"
dresult="./result/nonaddpgs/${dataname}/"
mkdir -p $dresult

fcv="${ddata}ukb_imp.cv"
fphe="${ddata}ukb_imp.phe"

bash ./script/lib/${program}/score.sh \
    $phe $cvi $sex $dresult $fphe $fcv

# TMP gzip
#f=${dresult}/${program}/${phe}/score.cv${cvi}/${program}.scorecov_regontr
#gzip -f $f
