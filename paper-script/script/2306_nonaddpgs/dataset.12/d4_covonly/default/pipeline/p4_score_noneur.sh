#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 4
#$ -q all.q
#$ -l hostname=a02|a03|z01|z02|g01|g02

set -eux

kind="default"
mode="noneur"
program="covonly"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

if [ -v SGE_TASK_ID ]; then
    job_id=$((SGE_TASK_ID - 1))
    threads=$NSLOTS
else
    job_id=$4
    threads=4
fi

source ./script/lib/system/python_threads.sh $threads

#dataname="$1"
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

fcv="${ddata}ukb_imp_noneur.cv"
fphe="${ddata}ukb_imp_noneur.phe"

bash ./script/lib/${program}/score-withcov_mode.sh \
    $phe $cvi $sex $dresult $fphe $fcv $kind $mode

# TMP gzip
#f=${dresult}/${program}/${phe}/score.cv${cvi}/${program}.scorecov_regontr
#gzip -f $f
