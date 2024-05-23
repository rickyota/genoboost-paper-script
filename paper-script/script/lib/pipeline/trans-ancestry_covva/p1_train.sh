#!/bin/bash

# TODO: allow mid_path="${kind} ${model}"
# for boosting
# by this, the number of ../../../ can be calculated by detecting space

set -eux

kind="trans-ancestry_covva"
#kind_ori="default"
#program="ldpred"
#dataname="dataset.12"

if [ $# -ne 6 ]; then
    echo "Wrong argument number."
    exit 1
fi

kind_ori=$1
program=$2
dataname=$3

hostname
pwd
echo "Start $(readlink -f $0)"

job_id=$6

phes="$4"
cvn="$5"

cvi=$((job_id % cvn))
phei=$((job_id / cvn))
phes_ar=(${phes})
phe=${phes_ar[$phei]}

dresult="./result/nonaddpgs/${dataname}/${program}/${kind}/"
dresult_ori="./result/nonaddpgs/${dataname}/${program}/${kind_ori}/"
mkdir -p $dresult

for anc in afr sas eas; do
    dout="${dresult}/${anc}/${phe}/train.cv${cvi}/"
    mkdir -p $dout

    dout_tr="${dresult_ori}/${phe}/train.cv${cvi}/"

    for fwgt_tr_abs in ${dout_tr}/*.wgt; do
        [ ! -f $fwgt_tr_abs ] && continue
        f=$(basename $fwgt_tr_abs)

        fwgt_tr="../../../../${kind_ori}/${phe}/train.cv${cvi}/${f}"

        #fwgt="${dout}/${f}"
        ln -sf $fwgt_tr $dout
        #ln -sf $fwgt_tr $fwgt
        #echo $fwgt
    done
done

#for phe in $phes; do
#    for cvi in $(seq 0 $((cvn - 1))); do
#        dout="${dresult}/${phe}/train.cv${cvi}/"
#        mkdir -p $dout
#
#        dout_tr="${dresult_ori}/${phe}/train.cv${cvi}/"
#        for fwgt_tr_abs in ${dout_tr}/*.wgt; do
#            f=$(basename $fwgt_tr_abs)
#            fwgt_tr="../../../${kind_ori}/${phe}/train.cv${cvi}/${f}"
#            fwgt="${dout}/${f}"
#            ln -sf $fwgt_tr $fwgt
#            #echo $fwgt
#        done
#    done
#done
