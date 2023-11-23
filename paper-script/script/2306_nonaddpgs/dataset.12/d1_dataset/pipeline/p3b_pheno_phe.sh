#!/bin/bash

set -eux

dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh
source ./script/lib/system/conda_activate.sh genetics

dcase="$1"
phes="$2"

ddata="./data/${dataname}/"
ddata_tmp="${ddata}/tmp/"
mkdir -p $ddata_tmp

# cov
#python -m projects_py.genetics.dataset.phenotype \
#	--case \
#	--fout ${ddata_tmp}/case.phenotypes \
#	--fgenot ${ddata}/ukb_imp \
#	--genot-format plink2 \
#	--dukb_case $dcase \
#	--phes $phes

source ./script/lib/system/conda_deactivate.sh

source ./script/lib/system/module.sh load csvkit

csvjoin -c1,1 -t -I -v ${ddata_tmp}/sex_age_pc.cov ${ddata_tmp}/case.phenotypes |
	csvformat -T >${ddata}/ukb_imp.phe

source ./script/lib/system/module.sh unload csvkit
