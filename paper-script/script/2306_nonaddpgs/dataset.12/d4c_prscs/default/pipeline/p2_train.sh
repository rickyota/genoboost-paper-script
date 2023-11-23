#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 28
#$ -q all.q
#$ -l hostname=a02|a03|z01|z02|g01|g02

set -eux

program="prscs"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh

if [ -v SGE_TASK_ID ]; then
    job_id=$((SGE_TASK_ID - 1))
    #threads=$NSLOTS
else
    job_id=$5
    #threads=16
fi

dataname="$1"
sex="$2" # male, female, both
phes="$3"
cvn="$4"

phei=$((job_id / cvn))
cvi=$((job_id % cvn))

phes_ar=(${phes})
phe=${phes_ar[$phei]}

seed=59
ddata="./data/${dataname}/"
dresult="./result/nonaddpgs/${dataname}/${program}/"
mkdir -p $dresult

dout="${dresult}/${phe}/train.cv${cvi}/raw/"
mkdir -p $dout

fgenot_ref="./data/1kg/eur.qc/1kg.eur.qc"
fcv="${ddata}ukb_imp.cv"

dtmp="${dout}/tmp"
fss="${dtmp}/ukb_imp.cv${cvi}.${phe}.glm.logistic.hybrid.ss"

if [ $sex == "female" ]; then
    samplen=$(awk -v cv=${cvi} '($2==0 && $(3+cv) == "tr"){ print $1 }' $fcv | awk 'END{print NR}')
elif [ $sex == "both" ]; then
    samplen=$(awk -v cv=${cvi} '($(3+cv) == "tr"){ print $1 }' $fcv | awk 'END{print NR}')
fi

echo $samplen

bash ./script/lib/${program}/train.sh $dout $fss $fgenot_ref $samplen $seed
