#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 4
#$ -q all.q
#$ -l hostname=a02|a03|z01|z02|g01|g02

set -eux

kind="monitor"
program="snpboost"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

threads=8

source ./script/lib/system/python_threads.sh $threads


sex="$1"
phes="$2"
cvn="$3"

#phei=$((job_id / cvn))
#cvi=$((job_id % cvn))
#phes_ar=(${phes})
#phe=${phes_ar[$phei]}

dresult="./result/nonaddpgs/${dataname}/${program}/${kind}/"
ddata="./data/${dataname}/"
dtmp="${ddata}/tmp/${program}/"
fgenot="${dtmp}ukb_imp"

source ./script/lib/system/conda_activate.sh genetics

for phe in $phes; do
    for cvi in $(seq 0 $((cvn - 1))); do
        dwgt="${dresult}/${phe}/train.cv${cvi}/raw/"
        dout="${dresult}/${phe}/train.cv${cvi}/"

        #python -m projects_py.genetics.main.snpnet_format_wgt \
        python -m projects_py.pgs.snpboost_format_wgt \
            --dout $dout \
            --dwgt $dwgt \
            --fgenot $fgenot \
            --sex $sex
        #exit 0

    done
done
