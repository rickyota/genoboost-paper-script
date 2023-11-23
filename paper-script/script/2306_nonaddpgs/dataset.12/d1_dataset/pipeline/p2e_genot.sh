#!/bin/bash

# create pgen file

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -t 1:22
#$ -pe smp 64
#$ -q 2x.q
#$ -l hostname=z01|z02|g01|g02|a02|a03

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
    chrom=$6
    threads=32
fi

dataname="$1"
dgenot="$2"
maf_thre="$3"
hwe_thre="$4"
geno_thre="$5"

ddata="./data/${dataname}/"
ddata_tmp="${ddata}/tmp/"
mkdir -p $ddata_tmp

fgenot="${ddata_tmp}/ukb_imp_chr${chrom}"

snvs_include="${ddata}ukb_imp.valid_info_hm3.snvs"

#cat $snvs_include | cut -f3 | awk 'BEGIN{OFS="\t"}{split($0,a,"_"); print a[0]}' | head
#cat $snvs_include | cut -f3 | awk 'BEGIN{OFS="\t"}{split($1,a,"_"); print $1,a[1]}' | head

## untested
# set major allele as reference
# remove  _{REF}_{ALT} from rs
plink2 --make-pgen vzs \
    --out $fgenot \
    --maj-ref \
    --threads ${threads} \
    --pfile ${dgenot}ukb_imp_chr${chrom} vzs \
    --keep ${ddata}/qc.sample \
    --extract <(cat $snvs_include | cut -f3 | awk 'BEGIN{OFS="\t"}{split($1,a,"_"); print a[1]}') \
    --maf $maf_thre \
    --hwe $hwe_thre \
    --geno $geno_thre \
    --update-name <(cat $snvs_include | cut -f3 | awk 'BEGIN{OFS="\t"}{split($1,a,"_"); print $1,a[1]}')

