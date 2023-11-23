#!/bin/bash

# align 1000g to ukb

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -sync y
##$ -t 1:22
#$ -pe smp 16
#$ -q 2x.q
#$ -l hostname=z01|z02|g01|g02|a02|a03

# for -t 1:22
# 2x.q means that these array occupy half of a02 and distribute job
# according to -pe 32, which means 2(=128/2/32) jobs run simultaneously.

set -eux

dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh
source ./script/lib/system/module.sh load plink2

if [ -v SGE_TASK_ID ]; then
    cvi=$((SGE_TASK_ID - 1))
    threads=$NSLOTS
else
    cvi=$2
    threads=16
fi

dataname="$1"

ddata="./data/${dataname}/"

fcv="${ddata}/ukb_imp.cv"

# TODO: female only
plink2 --freq \
    cols='chrom,pos,ref,alt,altfreq,nobs' \
    --out ${ddata}ukb_imp.cv${cvi} \
    --threads $threads \
    --pfile ${ddata}ukb_imp vzs --keep <(awk -v cv=${cvi} '($(3+cv) == "tr"){ print $1 }' $fcv)
