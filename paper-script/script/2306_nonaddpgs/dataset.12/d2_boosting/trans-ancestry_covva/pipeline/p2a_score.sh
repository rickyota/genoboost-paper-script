#!/bin/bash

# TODO: how to merge to lib/pipeline?

set -eux

kind="trans-ancestry_covva"
kind_ori="default"
program="boosting"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

phes="$1"
cvn="$2"
model="$3"

dresult="./result/nonaddpgs/${dataname}/${program}/${kind}/"
dresult_ori="./result/nonaddpgs/${dataname}/${program}/${kind_ori}/"
mkdir -p $dresult

for phe in $phes; do
    for cvi in $(seq 0 $((cvn - 1))); do
        for anc in afr sas eas; do
            dout="${dresult}/${model}/${anc}/${phe}/score.cv${cvi}/"
            mkdir -p $dout

            dout_ori="${dresult_ori}/${model}/${phe}/score_trans-ancestry.cv${cvi}/"

            for fscore_ori in ${dout_ori}/${program}*.score.gz; do
                [ ! -f $fscore_ori ] && continue
                fscore_name=$(basename $fscore_ori)

                fscore_rel="../../../../../${kind_ori}/${model}/${phe}/score_trans-ancestry.cv${cvi}/${fscore_name}"
                #fscore="${dout}/${fscore_name}"

                ln -sf $fscore_rel $dout
                #ln -sf $fscore_rel $fscore
            done
        done
    done
done
