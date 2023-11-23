#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 4
##$ -pe smp 8
##$ -pe smp 28
#$ -q all.q
#$ -l hostname=a02|a03

set -eux

kind="monitor"
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

flog="./log/nonaddpgs/${dataname}/${program}/${kind}/${phe}/score.cv${cvi}/score.$(date +'%y%m%d-%H%M%S').log"
mkdir -p "$(dirname $flog)"

ddata="./data/${dataname}/"
dtmp="${ddata}/tmp/${program}/"
#dtmp="${ddata}/tmp/${program}/${kind}/"
dresult="./result/nonaddpgs/${dataname}/${program}/${kind}/"

dwgt="${dresult}/${phe}/train.cv${cvi}/"

#dout="${dresult}/${phe}/score.cv${cvi}/raw/"
dout="${dresult}/${phe}/score.cv${cvi}/"
mkdir -p $dout

#fphe="${ddata}/ukb_imp.phe"

# or ${ddata}/ukb_imp
fgenot="${dtmp}ukb_imp"

bash ./script/lib/genetics/score.sh \
    $dout $dwgt $fgenot $threads |& tee $flog

gzip -f ${dout}/${program}.score

## test
#./projects_rust/target/release/genetics_res \
#    score  \
#    --verbose \
#    --threads $threads \
#    --dir-score . \
#    --file-genot ${fgenot} \
#    --genot-format "plink2-vzs" \
#    --file-phe ${fphe} --phe ${phe} \
#    --file-wgt ./snpnet_n-1_tmp.wgt
#
#    #--file-wgt ${dwgt}/snpnet_n-1.wgt
# --resume
