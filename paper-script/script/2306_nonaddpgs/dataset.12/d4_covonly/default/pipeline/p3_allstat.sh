#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -j y

#$ -pe smp 1
#$ -q all.q
#$ -l hostname=z02
##$ -l hostname=a02|a03|z01|z02|g01|g02

set -eux

kind='default'
program="covonly"

echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh
source ./script/lib/system/conda_activate.sh genetics

if [ -v SGE_TASK_ID ]; then
    job_id=$((SGE_TASK_ID - 1))
    threads=$NSLOTS
else
    job_id=$5
    threads=16
fi

source ./script/lib/system/python_threads.sh $threads

dataname="$1"
sex="$2" # calculate all scores including male for bc
phes="$3"
cvn="$4"

regon='tr'

phei=$job_id
phes_ar=(${phes})
phe=${phes_ar[$phei]}

ddata="./data/${dataname}/"
dresult="./result/nonaddpgs/${dataname}/"
mkdir -p $dresult

fcv="${ddata}/ukb_imp.cv"
fphe="${ddata}/ukb_imp.phe"

python -m projects_py.boosting.main.allstat \
    --dresult $dresult \
    --method $program \
    --regon $regon \
    --fcv $fcv --cvn $cvn \
    --fpheno $fphe --phe $phe \
    --sex $sex
