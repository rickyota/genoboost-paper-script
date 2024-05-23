#!/bin/bash

# create chrom pgen

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 30
#$ -t 1:22
##$ -q 2x.q
#$ -l hostname=z01|z02|g01|g02
##$ -l hostname=z01|z02|g01|g02|a02|a03

set -eux

dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh
source ./script/lib/system/conda_activate.sh genetics

if [ -v SGE_TASK_ID ]; then
	chrom=$SGE_TASK_ID
	threads=$NSLOTS
else
	chrom=$2
	threads=32
fi

dgenot="$1"

ddata="./data/${dataname}/"
mkdir -p $ddata

#for anc in afr sas eas; do
#	python -m projects_py.genetics.dataset.qc_sample \
#		--qc_sample \
#		--fout ${ddata}/qc.trans-ancestry.${anc}.sample \
#		--fgenot ${dgenot}/ukb_imp_chr% \
#		--genot-format plink2 \
#		--fukb $fukb_whole \
#		--fwithdrawn $fwithdrawn \
#		--ethnic-transancestry $anc
#done

# create .trans-ancestry.sample

#source ./script/lib/system/conda_deactivate.sh

#rm -rf ${ddata}/qc.trans-ancestry.sample
#for anc in afr sas eas; do
#cat ${ddata}/qc.trans-ancestry.${anc}.sample >>${ddata}/qc.trans-ancestry.sample
#done

source ./script/lib/system/module.sh load plink2/20230426
source ./script/lib/system/module.sh load csvkit

ddata_tmp="${ddata}/tmp/"
mkdir -p $ddata_tmp
fgenot="${ddata_tmp}/ukb_imp_trans-ancestry_chr${chrom}"

# TOFIX: somehow, rs1377872_G_A is duplicated in snvs_include.
snvs_include="${ddata}ukb_imp.valid_info_hm3.snvs"

zstdcat ${ddata}/ukb_imp.pvar.zst |
	plink2 --make-pgen vzs \
		--out $fgenot \
		--ref-allele force /dev/stdin 4 3 \
		--threads ${threads} \
		--pfile ${dgenot}ukb_imp_chr${chrom} vzs \
		--keep ${ddata}/qc.trans-ancestry.sample \
		--extract <(cat $snvs_include | cut -f3 | awk 'BEGIN{OFS="\t"}{split($1,a,"_"); print a[1]}') \
		--update-name <(cat $snvs_include | uniq | cut -f3 | awk 'BEGIN{OFS="\t"}{split($1,a,"_"); print $1,a[1]}')

source ./script/lib/system/module.sh unload csvkit
source ./script/lib/system/module.sh unload plink2/20230426
