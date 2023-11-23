#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 4
#$ -q all.q
#$ -l hostname=a02|a03|z01|z02|g01|g02

set -eux

kind="monitor"
mode="noneur"
program="snpboost"
dataname="dataset.12"

echo "ny"
exit 1

hostname
pwd
echo "Start $(readlink -f $0)"

if [ -v SGE_TASK_ID ]; then
    job_id=$((SGE_TASK_ID - 1))
    #threads=$NSLOTS
else
    job_id=$3
    #threads=16
fi


phes="$1"
cvn="$2"

cvi=$((job_id % cvn))
phei=$((job_id / cvn))
phes_ar=(${phes})
phe=${phes_ar[$phei]}

flog="./log/nonaddpgs/${dataname}/${program}/${kind}/${phe}/score.cv${cvi}/score-withcov.$(date +'%y%m%d-%H%M%S').log"
mkdir -p "$(dirname $flog)"

ddata="./data/${dataname}/"
dresult="./result/nonaddpgs/${dataname}/"

fcv="${ddata}ukb_imp_noneur.cv"
#fcv="${ddata}ukb_imp.cv"
fphe="${ddata}ukb_imp_noneur_pcname.phe"
#fphe="${ddata}ukb_imp.phe"

cov_names="sex,age,PC_1-PC_10"

python -m projects_py.boosting.main.withcov \
    --score \
    --method $program \
    --mid-path $kind \
    --dresult $dresult \
    --fcv $fcv --cvi $cvi \
    --fpheno $fphe --phe $phe \
    --cov-names $cov_names \
    --mode-score $mode |& tee $flog

#bash ./script/lib/${program}/score-withcov_mode.sh \
#    $phe $cvi $dresult $fphe $fcv $kind $mode |& tee $flog
