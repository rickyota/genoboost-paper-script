#!/bin/bash

set -eux

dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh
source ./script/lib/system/conda_activate.sh genetics

ddata="./data/${dataname}/"

python -m projects_py.genetics.dataset.sida_rs \
    --fout ${ddata}/ukb_imp.sida_rs.tsv \
    --fgenot ${ddata}/ukb_imp \
    --genot-format plink2vzs

source ./script/lib/system/conda_deactivate.sh
