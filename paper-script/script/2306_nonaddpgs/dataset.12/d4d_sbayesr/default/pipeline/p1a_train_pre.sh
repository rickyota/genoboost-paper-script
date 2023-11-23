#!/bin/bash

set -eux

program="sbayesr"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh
source ./script/lib/system/conda_activate.sh $program


phes="$1"
cvn="$2"
#fassoc="$2"
#fgenot_ref="$2"

flog="./log/nonaddpgs/${program}/$(basename $0).$(date +'%y%m%d-%H%M%S').log"
mkdir -p "$(dirname $flog)"

ddata="./data/${dataname}/"
dresult="./result/nonaddpgs/${dataname}/${program}/"

# TODO: freq only for female

# create .ma
for cvi in $(seq 0 $((cvn - 1))); do
    ffreq="${ddata}ukb_imp.cv${cvi}.afreq"
    for phe in $phes; do
        fss="./result/nonaddpgs/${dataname}/assoc/ukb_imp.cv${cvi}.${phe}.glm.logistic.hybrid.ss"

        dout="${dresult}/${phe}/train.cv${cvi}/raw/"
        dtmp="${dout}/tmp"
        mkdir -p $dtmp

        fma="${dtmp}${phe}.cv${cvi}.ma"
        python ./lib/sbayesr/src/create_ma.py $fma ${fss} ${ffreq} |& tee ${flog}
    done
done
