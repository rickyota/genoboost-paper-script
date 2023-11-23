#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 16
#$ -q all.q
#$ -l hostname=a02|a03

set -eux

kind="default"
mode="noneur"
program="snpnet"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh

if [ -v SGE_TASK_ID ]; then
    job_id=$((SGE_TASK_ID - 1))
    threads=$NSLOTS
else
    job_id=$4
    threads=16
fi

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

if [ $HPCL_HOST = "pg" ]; then
    fgenot="/grid2/ricky/data/${dataname}/ukb_imp_noneur"
else
    fgenot="${ddata}/ukb_imp_noneur"
fi

#bash ./script/lib/genetics/score_concat.sh \
bash ./script/lib/genetics/score_concat_use-snv-pos.sh \
    $dout $dwgt $fgenot $threads "n" |& tee $flog

# TODO: gzip in python or score_concat.sh
gzip -f ${dout}/${program}_n.score

## test
#./projects_rust/target/release/genetics_res \
#    score --resume \
#    --verbose \
#    --threads $threads \
#    --dir-score . \
#    --file-genot ${fgenot} \
#    --genot-format "plink2-vzs" \
#    --file-phe ${fphe} --phe ${phe} \
#    --file-wgt ./snpnet_n-1_tmp.wgt
#
#    #--file-wgt ${dwgt}/snpnet_n-1.wgt
