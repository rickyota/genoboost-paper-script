#!/bin/bash

set -eux

dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh
source ./script/lib/system/conda_activate.sh genetics

fukb_whole="$1"
year_now="$2"
month_now="$3"

ddata="./data/${dataname}/"
ddata_tmp="${ddata}/tmp/"
mkdir -p $ddata_tmp

# cov
python -m projects_py.genetics.dataset.phenotype \
	--cov \
	--fout ${ddata_tmp}/sex_age_pc_noneur.cov \
	--fgenot ${ddata}/ukb_imp_noneur \
	--genot-format plink2 \
	--fukb $fukb_whole \
	--cov_column sex,age,PC_1-PC_10 \
	--now_date $year_now $month_now

source ./script/lib/system/conda_deactivate.sh
