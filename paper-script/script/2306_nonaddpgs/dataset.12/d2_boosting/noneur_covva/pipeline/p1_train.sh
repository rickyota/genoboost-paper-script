#!/bin/bash

set -eux

kind="noneur_covva"
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
        dout="${dresult}/${model}/${phe}/train.cv${cvi}/"
        mkdir -p $dout

        dout_tr="${dresult_ori}/${model}/${phe}/train.cv${cvi}/"
        for fwgt_tr_abs in ${dout_tr}/*.wgt; do
            f=$(basename $fwgt_tr_abs)

            fwgt_tr="../../../../${kind_ori}/${model}/${phe}/train.cv${cvi}/${f}"
            fwgt="${dout}/${f}"

            ln -sf $fwgt_tr $fwgt
            #echo $fwgt $fwgt_tr
        done
    done
done

#rhos="1.0 0.1 0.01 0.001"
#
#for phe in $phes; do
#    for cvi in $(seq 0 $((cvn - 1))); do
#        dout="${dresult}/${phe}/train.cv${cvi}/"
#        mkdir -p $dout
#
#        #dout_tr="${dresult_ori}/${phe}/train.cv${cvi}/"
#
#        for rho in $rhos; do
#            fwgt_tr="../../../${kind_ori}/${phe}/train.cv${cvi}/${program}_rho-${rho}.wgt"
#            #fwgt_tr=$(readlink -f $fwgt_tr)
#            fwgt="${dout}/${program}_rho-${rho}.wgt"
#            ln -sf $fwgt_tr $fwgt
#        done
#    done
#done
