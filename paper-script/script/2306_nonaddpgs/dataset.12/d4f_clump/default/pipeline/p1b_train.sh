#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 16
#$ -q all.q
#$ -l hostname=a02|a03|z01|z02|g01|g02

set -eux

program="clump"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

if [ -v SGE_TASK_ID ]; then
    job_id=$((SGE_TASK_ID - 1))
    threads=$NSLOTS
else
    job_id=$6
    threads=16
fi

sex="$1" # male, female, both
phes="$2"
cvn="$3"
clump_p_max="$4"
clump_r2s="$5"

phei=$((job_id / cvn))
cvi=$((job_id % cvn))

phes_ar=(${phes})
phe=${phes_ar[$phei]}

ddata="./data/${dataname}/"
dresult="./result/nonaddpgs/${dataname}/${program}/"
mkdir -p $dresult

dout="${dresult}/${phe}/train.cv${cvi}/raw/"
mkdir -p $dout

fgenot="${ddata}ukb_imp"
fphe="${ddata}ukb_imp.plink1.phe"
fcv="${ddata}ukb_imp.cv"

fassoc="./result/nonaddpgs/${dataname}/assoc/ukb_imp.cv${cvi}.${phe}.glm.logistic.hybrid"
fcv="${ddata}ukb_imp.cv"

bash ./script/lib/${program}/train.sh \
    $phe $cvi $dout $fgenot $fphe $fassoc $fcv $sex $threads "$clump_p_max" "$clump_r2s"
