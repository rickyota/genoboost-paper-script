#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 28
#$ -q all.q
#$ -l hostname=a02|a03|z01|z02|g01|g02
##$ -l hostname=a02|a03

set -eux

kind="ex-chr6-ts"
program="ldpred"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

#if [ -v SGE_TASK_ID ]; then
#    job_id=$((SGE_TASK_ID - 1))
#    threads=$NSLOTS
#else
#    job_id=$5
#    threads=16
#fi

#dataname="$1"
#fgenot_ref="$2"
#sex="$2" # male, female, both
phes="$1"
cvn="$2"

#cvi=$((job_id % cvn))
#phei=$((job_id / cvn))
#phes_ar=(${phes})
#phe=${phes_ar[$phei]}

#flog="./log/nonaddpgs/${dataname}/${program}/${kind}/${phe}/train.cv${cvi}/train.$(date +'%y%m%d-%H%M%S').log"
#mkdir -p "$(dirname $flog)"

#ddata="./data/${dataname}/"
dresult="./result/nonaddpgs/${dataname}/${program}/${kind}/"
mkdir -p $dresult

rhos="1.0 0.1 0.01 0.001"

for phe in $phes; do
    for cvi in $(seq 0 $((cvn - 1))); do
        dout="${dresult}/${phe}/train.cv${cvi}/"
        mkdir -p $dout

        dout_tr="./result/nonaddpgs/${dataname}/${program}/default/${phe}/train.cv${cvi}/"

        for rho in $rhos; do
            fwgt_tr="${dout_tr}/${program}_rho-${rho}.wgt"
            fwgt="${dout}/${program}_rho-${rho}.wgt"
            cat $fwgt_tr | awk '(NR==1 || $1!=6)' >$fwgt
        done
    done
done
