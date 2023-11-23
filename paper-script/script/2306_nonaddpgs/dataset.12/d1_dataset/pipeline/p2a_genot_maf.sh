#!/bin/bash

# first narrow the candidates and speed up the following qc
# extracted by
# - maf
# - hwe
# -geno
# this should not affect the following results

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -sync y
#$ -t 1:22
#$ -pe smp 16
#$ -q 2x.q
#$ -l hostname=z01|z02|g01|g02|a02|a03

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
    chrom=$SGE_TASK_ID
    threads=$NSLOTS
else
    chrom=$5
    threads=16
fi

dgenot="$1"
maf_thre="$2"
hwe_thre="$3"
geno_thre="$4"

ddata="./data/${dataname}/"
ddata_tmp="${ddata}/tmp/"
mkdir -p $ddata_tmp

plink2 --make-just-pvar zs \
    --out ${ddata_tmp}ukb_imp_chr${chrom}.maf \
    --threads ${threads} \
    --pfile ${dgenot}ukb_imp_chr${chrom} \
    --keep ${ddata}/qc.sample \
    --maf $maf_thre \
    --hwe $hwe_thre \
    --geno $geno_thre
