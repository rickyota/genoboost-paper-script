#!/bin/bash

# assumed to be called from ../dataset.12.sh

set -eux

dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh
source ./script/lib/system/conda_activate.sh genetics

dgenot="$1"
fukb_whole="$2"
fwithdrawn="$3"

ddata="./data/${dataname}/"
mkdir -p $ddata

python -m projects_py.genetics.dataset.qc_sample \
	--valid \
	--fout ${ddata}/qc.sample \
	--fgenot ${dgenot}/ukb_imp_chr% \
	--genot-format plink2 \
	--fukb $fukb_whole \
	--fwithdrawn $fwithdrawn

