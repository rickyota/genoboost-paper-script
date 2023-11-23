#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 28
#$ -q all.q
#$ -l hostname=a02|a03|z01|z02|g01|g02

set -eux

kind="ex-chr6-ts"
program="snpnet"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

#if [ -v SGE_TASK_ID ]; then
#    job_id=$((SGE_TASK_ID - 1))
#    threads=$NSLOTS
#else
#    job_id=$5
#    threads=32
#fi

#dataname="$1"
#sex="$2" # male, female, both
phes="$1"
cvn="$2"

#phei=$((job_id / cvn))
#cvi=$((job_id % cvn))
#phes_ar=(${phes})
#phe=${phes_ar[$phei]}

#ddata="./data/${dataname}/"
#dtmp="${ddata}/tmp/${program}/"
dresult="./result/nonaddpgs/${dataname}/${program}/${kind}"
mkdir -p $dresult

for phe in $phes; do
    for cvi in $(seq 0 $((cvn - 1))); do
        dout="${dresult}/${phe}/train.cv${cvi}/"
        mkdir -p $dout

        dout_tr="./result/nonaddpgs/${dataname}/${program}/default/${phe}/train.cv${cvi}/"

        for f in ${dout_tr}/*.wgt; do
            [ ! -f $f ] && continue

            b=$(basename $f)
            fwgt="${dout}/${b}"
            cat $f | awk '(NR==1 || $2!=6)' >$fwgt

            #exit 0
        done
    done
done
