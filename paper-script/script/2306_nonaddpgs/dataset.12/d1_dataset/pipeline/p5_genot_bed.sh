#!/bin/bash

# create pgen file
# it seems unavailable to simultaneously extract and merge

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -sync y
#$ -pe smp 16
##$ -q 2x.q
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

if [ -v SGE_O_WORKDIR ]; then
    #chrom=$SGE_TASK_ID
    threads=$NSLOTS
else
    #chrom=$6
    threads=16
fi

ddata="./data/${dataname}/"
ddata_tmp="${ddata}/tmp/"
mkdir -p $ddata_tmp

fgenot="${ddata}/ukb_imp"

# set major allele as reference
plink2 --make-bed \
    --keep-allele-order \
    --out $fgenot \
    --pfile $fgenot vzs \
    --threads ${threads}

fphe="${ddata}/ukb_imp.phe"
fphe_plink1="${ddata}/ukb_imp.plink1.phe"

# create .phe with FID
cat $fphe | awk 'BEGIN {OFS = "\t"} {print $1,$0}' | sed -e "1 s/IID/FID/" >$fphe_plink1
