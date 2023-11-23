#!/bin/bash

set -eux

kind="ex-chr6-ts"
program="lassosum"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

phes="$1"
cvn="$2"

dresult="./result/nonaddpgs/${dataname}/${program}/${kind}/"
kind_ori="default"

for phe in $phes; do
    for cvi in $(seq 0 $((cvn - 1))); do
        dout_ori="./result/nonaddpgs/${dataname}/${program}/${kind_ori}//${phe}/train.cv${cvi}/"
        dout="${dresult}/${phe}/train.cv${cvi}/"
        for fwgtcov in ${dout_ori}/*.wgtcov_regonva; do
            [ ! -f $fwgtcov ] && continue

            b=$(basename $fwgtcov)
            ln -sf ../../../${kind_ori}/${phe}/train.cv${cvi}/${b} ${dout}

        done
    done
done
