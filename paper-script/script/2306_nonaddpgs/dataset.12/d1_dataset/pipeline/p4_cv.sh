#!/bin/bash

set -eux

dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh
source ./script/lib/system/conda_activate.sh genetics

dataname="$1"
prop_test="$2"
cvn="$3"
seed="$4"

ddata="./data/${dataname}/"
ddata_tmp="${ddata}/tmp/"
mkdir -p $ddata_tmp

python -m projects_py.genetics.dataset.cross_vali \
	--cv \
	--fout ${ddata}/ukb_imp.cv \
	--fgenot ${ddata}/ukb_imp \
	--genot-format plink2 \
	--prop_test $prop_test \
	--cvn $cvn \
	--seed $seed

source ./script/lib/system/conda_deactivate.sh
