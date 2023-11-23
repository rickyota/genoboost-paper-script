#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 4
#$ -q all.q
#$ -l hostname=a02|a03|z01|z02|g01|g02

set -eux

kind="nonadd"
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
#    threads=8
#fi

threads=8

source ./script/lib/system/python_threads.sh $threads

#dataname="$1"
sex="$1"
phes="$2"
cvn="$3"

#phei=$((job_id / cvn))
#cvi=$((job_id % cvn))

#phes_ar=(${phes})
#phe=${phes_ar[$phei]}

# TODO
#flog="./log/nonaddpgs/${dataname}/${program}/${kind}/${phe}/train.cv${cvi}/p1c_post.$(date +'%y%m%d-%H%M%S').log"
#mkdir -p "$(dirname $flog)"

dresult="./result/nonaddpgs/${dataname}/${program}/${kind}/"
ddata="./data/${dataname}/"
dtmp="${ddata}/tmp/${program}/${kind}/"
fgenot="${dtmp}ukb_imp"

source ./script/lib/system/conda_activate.sh genetics

for phe in $phes; do
    for cvi in $(seq 0 $((cvn - 1))); do
        dwgt="${dresult}/${phe}/train.cv${cvi}/raw/"
        dout="${dresult}/${phe}/train.cv${cvi}/"

        #python -m projects_py.genetics.main.snpnet_format_wgt \
        python -m projects_py.pgs.snpnet_format_wgt \
            --dout $dout \
            --dwgt $dwgt \
            --fgenot $fgenot \
            --sex $sex
    done
done

source ./script/lib/system/conda_deactivate.sh
