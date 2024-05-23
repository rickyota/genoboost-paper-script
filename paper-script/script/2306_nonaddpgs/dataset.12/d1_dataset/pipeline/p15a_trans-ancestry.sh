#!/bin/bash

# create chrom pgen

#$ -S /bin/bash
#$ -cwd
#$ -j y
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
source ./script/lib/system/conda_activate.sh genetics

if [ -v SGE_TASK_ID ]; then
	#chrom=$SGE_TASK_ID
	threads=$NSLOTS
else
	#chrom=$4
	threads=32
fi

dgenot="$1"
fukb_whole="$2"
fwithdrawn="$3"

ddata="./data/${dataname}/"
mkdir -p $ddata

for anc in afr sas eas; do
	python -m projects_py.genetics.dataset.qc_sample \
		--qc_sample \
		--fout ${ddata}/qc.trans-ancestry.${anc}.sample \
		--fgenot ${dgenot}/ukb_imp_chr% \
		--genot-format plink2 \
		--fukb $fukb_whole \
		--fwithdrawn $fwithdrawn \
		--ethnic-transancestry $anc
done

# create .trans-ancestry.sample

source ./script/lib/system/conda_deactivate.sh

rm -rf ${ddata}/qc.trans-ancestry.sample
for anc in afr sas eas; do
	cat ${ddata}/qc.trans-ancestry.${anc}.sample >>${ddata}/qc.trans-ancestry.sample
done
