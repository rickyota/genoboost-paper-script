#!/bin/bash

set -eux

kind="default"
program="prscs"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh
source ./script/lib/system/conda_activate.sh $program

phes="$1"
cvn="$2"

flog="./log/nonaddpgs/${program}/train_pre.$(date +'%y%m%d-%H%M%S').log"
mkdir -p "$(dirname $flog)"

dresult="./result/nonaddpgs/${dataname}/${program}/"

for cvi in $(seq 0 $((cvn - 1))); do
    for phe in $phes; do
        fss_ori="./result/nonaddpgs/${dataname}/assoc/ukb_imp.cv${cvi}.${phe}.glm.logistic.hybrid.ss"

        dout="${dresult}/${phe}/train.cv${cvi}/raw/"
        dtmp="${dout}/tmp"
        mkdir -p $dtmp
        fss="${dtmp}/ukb_imp.cv${cvi}.${phe}.glm.logistic.hybrid.ss"

        echo -e "SNP\tA1\tA2\tOR\tP" >$fss
        cat $fss_ori | tail -n +2 | awk '{ print $3"\t"$6"\t"$17"\t"$10"\t"$15 }' >>$fss

    done
done
