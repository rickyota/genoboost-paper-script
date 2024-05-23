#!/bin/bash

# create chrom pgen

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -sync y
#$ -pe smp 30
##$ -pe smp 110
#$ -l hostname=z01|z02|g01|g02
##$ -l hostname=z01|z02|g01|g02|a02|a03

set -eux

dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh
source ./script/lib/system/module.sh load plink2/20230426

if [ -v SGE_TASK_ID ]; then
	# this works even when -t 1:22 is not set
	# SGE_TASK_ID = "undefined"
	threads=$NSLOTS
else
	threads=32
fi

ddata="./data/${dataname}/"
mkdir -p $ddata

source ./script/lib/system/module.sh load csvkit

ddata_tmp="${ddata}/tmp/"
fgenot="${ddata}/ukb_imp_tras-ancestry"

plink2 --make-pgen vzs \
	--out $fgenot \
	--pmerge-list <(printf "${ddata_tmp}ukb_imp_trans-ancestry_chr%s\n" {1..22}) pfile-vzs \
	--pmerge-output-vzs \
	--indiv-sort none \
	--keep-allele-order \
	--threads ${threads}

source ./script/lib/system/module.sh unload plink2/20230426
