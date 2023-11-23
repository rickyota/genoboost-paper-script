#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 16
#$ -q all.q
#$ -l hostname=a02|a03

set -eux

kind="monitor"
mode="noneur"
program="snpboost"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh

if [ -v SGE_TASK_ID ]; then
    job_id=$((SGE_TASK_ID - 1))
    threads=$NSLOTS
else
    job_id=$3
    threads=16
fi


#sex="$2" # calculate all scores including male for bc
phes="$1"
cvn="$2"

cvi=$((job_id % cvn))
phei=$((job_id / cvn))
phes_ar=(${phes})
phe=${phes_ar[$phei]}

flog="./log/nonaddpgs/${dataname}/${program}/${kind}/${phe}/score_${mode}.cv${cvi}/score.$(date +'%y%m%d-%H%M%S').log"
mkdir -p "$(dirname $flog)"

ddata="./data/${dataname}/"
dresult="./result/nonaddpgs/${dataname}/${program}/${kind}/"

dwgt="${dresult}/${phe}/train.cv${cvi}/"

dout="${dresult}/${phe}/score_${mode}.cv${cvi}/"
mkdir -p $dout

#fphe="${ddata}/ukb_imp.phe"

if [ $HPCL_HOST = "pg" ]; then
    fgenot="/grid2/ricky/data/${dataname}/ukb_imp_noneur"
else
    fgenot="${ddata}/ukb_imp_noneur"
fi

bash ./script/lib/genetics/score.sh \
    $dout $dwgt $fgenot $threads |& tee $flog

#bash ./script/lib/genetics/score_concat.sh \
#bash ./script/lib/genetics/score_concat_use-snv-pos.sh \
# $dout $dwgt $fgenot $threads "n" |& tee $flog

# TODO: gzip in python or score_concat.sh
gzip -f ${dout}/${program}.score
#gzip -f ${dout}/${program}_n.score
