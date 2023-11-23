#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 28
#$ -q all.q
#$ -l hostname=a02|a03|z01|z02|g01|g02

## Extremely slow at a02, a03  or simultaneously running at the same node

set -eux

kind="monitor"
program="snpboost"
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

fgenot="${dtmp}ukb_imp"
fphe="${dtmp}ukb_imp.cv.phe"
#fphe="${dtmp}ukb_imp.phe"

# use standardized cov
# -> no effect
#fphe="${dtmp}ukb_imp_cov_standard.phe"

if [ $sex == "both" ]; then
    fsamples="${dtmp}/ukb_imp.cv${cvi}.monitor.samples"
    #fsamples_va="${dtmp}/ukb_imp.cv${cvi}.va.samples"
elif [ $sex == "female" ]; then
    fsamples="${dtmp}/ukb_imp.cv${cvi}.female.monitor.samples"
    #fsamples_va="${dtmp}/ukb_imp.cv${cvi}.female.va.samples"
fi

bash ./script/lib/${program}/train_monitor.sh \
    $phe $dout $fgenot $fphe $fsamples $sex $cvi $threads |& tee $flog
