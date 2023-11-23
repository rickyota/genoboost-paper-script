#!/bin/bash

set -eux

program="clump"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh
source ./script/lib/system/conda_activate.sh genetics

phes="$1"
cvn="$2"
clump_p_max="$3"
clump_r2s="$4"

ddata="./data/${dataname}/"
dresult="./result/nonaddpgs/${dataname}/${program}/"
mkdir -p $dresult

for cvi in $(seq 0 $((cvn - 1))); do
    for phe in $phes; do
        dout="${dresult}/${phe}/train.cv${cvi}/"

        for clump_r2 in $clump_r2s; do
            echo "cvi, phe, clump_p,clump_r2:$cvi , $phe , $clump_r2"

            fclump="${dout}/clump.p${clump_p_max}.r${clump_r2}.clumped"
            fwgt="${dout}/clump.p${clump_p_max}.r${clump_r2}.wgt"
            fss="./result/nonaddpgs/${dataname}/assoc/ukb_imp.cv${cvi}.${phe}.glm.logistic.hybrid.ss"

            #TODO: female only freq
            python -m projects_py.pgs.clump \
                --convert-wgt \
                --fout $fwgt \
                --fclump $fclump \
                --ffreq ${ddata}ukb_imp.cv${cvi}.afreq \
                --fss $fss

        done

    done
done
