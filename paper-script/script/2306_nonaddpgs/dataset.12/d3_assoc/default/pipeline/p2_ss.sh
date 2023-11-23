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

source ./script/lib/system/conda_activate.sh genetics

# add A2
# remove '#' from header
for cvi in $(seq 0 $((cvn - 1))); do
    for phe in $phes; do
        echo "phe $phe"
        #for gmodel in "add"; do
        for gmodel in $gmodels; do

            echo "gmodel $gmodel"
            if [ $gmodel == "add" ]; then
                fassoc="${dresult}ukb_imp.cv${cvi}.${phe}.glm.logistic.hybrid"
                fss="${dresult}ukb_imp.cv${cvi}.${phe}.glm.logistic.hybrid.ss"
            else
                fassoc="${dresult}ukb_imp.cv${cvi}.${phe}.${gmodel}.glm.logistic.hybrid"
                fss="${dresult}ukb_imp.cv${cvi}.${phe}.${gmodel}.glm.logistic.hybrid.ss"
            fi
            #fassoc="${dresult}ukb_imp.cv${cvi}.${phe}.glm.logistic.hybrid"
            #fss="${dresult}ukb_imp.cv${cvi}.${phe}.glm.logistic.hybrid.ss"
            # add p from logp
            [ ! -f "$fassoc" ] && continue

            python -m projects_py.genetics.main.format_ss \
                --add_a2 \
                --fout $fss \
                --fassoc $fassoc

        done
    done
done

# old from log10(p) -> p
#for cvi in $(seq 0 $((cvn - 1))); do
#    for phe in $phes; do
#        echo "phe $phe"
#        for gmodel in "add"; do
#            #for gmodel in $gmodels; do
#            echo "gmodel $gmodel"
#            fassoc="${dresult}ukb_imp.${phe}.glm.logistic.hybrid"
#            fss="${dresult}ukb_imp.${phe}.glm.logistic.hybrid.ss"
#            # add p from logp
#
#            python -m projects_py.genetics.assoc.assoc \
#                --add_col_p \
#                --fout $fss \
#                --fassoc $fassoc
#
#        done
#    done
#done

source ./script/lib/system/conda_deactivate.sh
