#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 4
#$ -q all.q
#$ -l hostname=a02|a03|z01|z02|g01|g02

set -eux

kind="noneur_covva"
program="covonly"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/conda_activate.sh genetics

if [ -v SGE_TASK_ID ]; then
    job_id=$((SGE_TASK_ID - 1))
    #threads=$NSLOTS
else
    job_id=$4
    #job_id=$5
    #threads=16
fi

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

fcv="${ddata}ukb_imp_noneur_vats.cv"
fphe="${ddata}ukb_imp_noneur.phe"

if [ $sex == "female" ]; then
    cov_names="age,PC_1-PC_10"
elif [ $sex == "both" ]; then
    cov_names="sex,age,PC_1-PC_10"
fi

# regon tr
python -m projects_py.boosting.main.withcov \
    --score \
    --regon va \
    --method $program \
    --mid-path $kind \
    --dresult $dresult \
    --fcv $fcv --cvi $cvi \
    --fpheno $fphe --phe $phe \
    --cov-name $cov_names \
    --sex $sex

# TMP gzip
#f=${dresult}/${program}/${phe}/score.cv${cvi}/${program}.scorecov_regontr
#gzip -f $f
