#!/bin/bash

set -eux

dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh

#dataname="$1"
phes="$1"
gmodels="$2"
cvn="$3"

#ddata="./data/${dataname}/"
dresult="./result/nonaddpgs/${dataname}/assoc/"
dresult_tmp="${dresult}/tmp/"

fbatch="${dresult}/snvs.batch"

batchn=$(cat $fbatch | tail -n1 | cut -f6)

check_all_batch_exist() {
    phe=$1
    gmodel=$2
    cvi=$3
    batchn=$4
    for batchi in $(seq 0 ${batchn}); do
        if [ $gmodel == "add" ]; then
            f=${dresult_tmp}ukb_imp.cv${cvi}.batch${batchi}.${phe}.glm.logistic.hybrid.zst
        else
            f=${dresult_tmp}ukb_imp.cv${cvi}.batch${batchi}.${gmodel}.${phe}.glm.logistic.hybrid.zst
        fi
        #f=${dresult_tmp}ukb_imp.cv${cvi}.batch${batchi}.${phe}.glm.logistic.hybrid.zst
        # echo "1" if any does not exit
        [ ! -f "$f" ] && echo "1" && exit 0
    done
    echo "0"
}

for cvi in $(seq 0 $((cvn - 1))); do
    for phe in $phes; do
        echo "phe $phe"
        #for gmodel in "add"; do
        for gmodel in $gmodels; do
            echo "gmodel $gmodel"

            # check if all files exist
            is_exist=$(check_all_batch_exist $phe $gmodel $cvi $batchn)
            if [[ $is_exist -eq 1 ]]; then
                continue
            fi

            #fassoc="${dresult}ukb_imp.cv${cvi}.${phe}.glm.logistic.hybrid"

            batchi=0
            if [ $gmodel == "add" ]; then
                f=${dresult_tmp}ukb_imp.cv${cvi}.batch${batchi}.${phe}.glm.logistic.hybrid.zst
                fassoc="${dresult}ukb_imp.cv${cvi}.${phe}.glm.logistic.hybrid"
            else
                f=${dresult_tmp}ukb_imp.cv${cvi}.batch${batchi}.${gmodel}.${phe}.glm.logistic.hybrid.zst
                fassoc="${dresult}ukb_imp.cv${cvi}.${phe}.${gmodel}.glm.logistic.hybrid"
            fi

            #zstdcat ${dresult_tmp}ukb_imp.cv${cvi}.batch0.${phe}.glm.logistic.hybrid.zst | head -n 1 >$fassoc
            zstdcat $f | head -n 1 >$fassoc
            for batchi in $(seq 0 ${batchn}); do
                if [ $gmodel == "add" ]; then
                    f=${dresult_tmp}ukb_imp.cv${cvi}.batch${batchi}.${phe}.glm.logistic.hybrid.zst
                else
                    f=${dresult_tmp}ukb_imp.cv${cvi}.batch${batchi}.${gmodel}.${phe}.glm.logistic.hybrid.zst
                fi
                zstdcat $f | tail -n +2 >>$fassoc
                #zstdcat ${dresult_tmp}ukb_imp.cv${cvi}.batch${batchi}.${phe}.glm.logistic.hybrid.zst | tail -n +2 >>$fassoc
            done

            # TODO: gzip -f $fassoc

            # TODO: rm tmp/

        done
    done
done
