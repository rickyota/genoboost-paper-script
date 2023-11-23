#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 2
#$ -q all.q
#$ -l hostname=a02|a03|z01|z02|g01|g02

set -eux

kind="ex-chr6-ts"
program="snpboost"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh
source ./script/lib/system/conda_activate.sh genetics

if [ -v SGE_TASK_ID ]; then
    job_id=$((SGE_TASK_ID - 1))
    threads=$NSLOTS
else
    job_id=$4
    threads=16
fi

source ./script/lib/system/python_threads.sh $threads

sex="$1"
phes="$2"
cvn="$3"

phei=$job_id
phes_ar=(${phes})
phe=${phes_ar[$phei]}

flog="./log/nonaddpgs/${dataname}/${program}/${kind}/${phe}/allstat.$(date +'%y%m%d-%H%M%S').log"
mkdir -p "$(dirname $flog)"

ddata="./data/${dataname}/"
dresult="./result/nonaddpgs/${dataname}/"
mkdir -p $dresult

fcv="${ddata}/ukb_imp.cv"
fphe="${ddata}/ukb_imp.phe"

python -m projects_py.boosting.main.allstat \
    --dresult $dresult \
    --method $program \
    --mid-path $kind \
    --fcv $fcv --cvn $cvn \
    --fpheno $fphe --phe $phe \
    --sex $sex |& tee $flog
