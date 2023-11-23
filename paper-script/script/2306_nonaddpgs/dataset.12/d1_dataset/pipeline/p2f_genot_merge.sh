#!/bin/bash

# create pgen file
# it seems unavailable to simultaneously extract and merge

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -sync y
#$ -pe smp 110
##$ -q 2x.q
#$ -l hostname=z01|z02|g01|g02|a02|a03

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

dgenot="$1"
maf_thre="$2"
hwe_thre="$3"
geno_thre="$4"

ddata="./data/${dataname}/"
ddata_tmp="${ddata}/tmp/"

fgenot="${ddata}/ukb_imp"

# set major allele as reference
plink2 --make-pgen vzs \
    --out $fgenot \
    --pmerge-list <(printf "${ddata_tmp}ukb_imp_chr%s\n" {1..22}) pfile-vzs \
    --pmerge-output-vzs \
    --indiv-sort none \
    --maj-ref \
    --threads ${threads}

