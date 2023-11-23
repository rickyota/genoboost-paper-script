#!/bin/bash

set -eux

program="ldpred"
kind="ex-chr6-gwasmethods"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"


phes="$1"
cvn="$2"

dresult="./result/nonaddpgs/${dataname}/${program}/${kind}/"

rhos_all="1.0 0.1 0.01 0.001"

for phe in $phes; do
    for cvi in $(seq 0 $((cvn - 1))); do

        dout="${dresult}/${phe}/score.cv${cvi}/"
        draw="${dout}/raw/"
        dtmp="${dout}/tmp/"
        mkdir -p $dtmp

        for rho in $rhos_all; do
            cp ${draw}/${program}_rho-${rho}.wgt.score ${dtmp}/${program}_rho-${rho}.score
            gzip -f ${dtmp}/${program}_rho-${rho}.score
        done

    done
done
