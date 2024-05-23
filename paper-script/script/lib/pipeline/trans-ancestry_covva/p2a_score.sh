#!/bin/bash

set -eux

kind="trans-ancestry_covva"
#kind_ori="default"
#program="clump"
#dataname="dataset.12"

if [ $# -ne 6 ]; then
    echo "Wrong argument number."
    exit 1
fi

kind_ori=$1
program=$2
dataname=$3

#program=$1
#dataname=$2

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
    dout="${dresult}/${anc}/${phe}/score.cv${cvi}/"
    #dout="${dresult}/${phe}/${anc}/score.cv${cvi}/"
    #dout="${dresult}/${phe}/score_${anc}.cv${cvi}/"
    #dout="${dresult}/${phe}/score.cv${cvi}/"
    mkdir -p $dout

    #dout_ori="${dresult_ori}/${phe}/score_${anc}.cv${cvi}/"
    dout_ori="${dresult_ori}/${phe}/score_trans-ancestry.cv${cvi}/"

    for fscore_ori in ${dout_ori}/${program}*.score.gz; do
        [ ! -f $fscore_ori ] && continue
        fscore_name=$(basename $fscore_ori)

        fscore_rel="../../../../${kind_ori}/${phe}/score_trans-ancestry.cv${cvi}/${fscore_name}"
        #fscore_rel="../../../${kind_ori}/${phe}/score_trans-ancestry.cv${cvi}/${fscore_name}"
        #fscore="${dout}/${fscore_name}"

        ln -sf $fscore_rel $dout
        #ln -sf $fscore_rel $fscore
    done
done

#for phe in $phes; do
#    for cvi in $(seq 0 $((cvn - 1))); do
#        for anc in afr sas eas; do
#            dout="${dresult}/${anc}/${phe}/score.cv${cvi}/"
#            #dout="${dresult}/${phe}/${anc}/score.cv${cvi}/"
#            #dout="${dresult}/${phe}/score_${anc}.cv${cvi}/"
#            #dout="${dresult}/${phe}/score.cv${cvi}/"
#            mkdir -p $dout
#
#            #dout_ori="${dresult_ori}/${phe}/score_${anc}.cv${cvi}/"
#            dout_ori="${dresult_ori}/${phe}/score_trans-ancestry.cv${cvi}/"
#
#            for fscore_ori in ${dout_ori}/${program}*.score.gz; do
#                [ ! -f $fscore_ori ] && continue
#                fscore_name=$(basename $fscore_ori)
#
#                fscore_rel="../../../../${kind_ori}/${phe}/score_trans-ancestry.cv${cvi}/${fscore_name}"
#                #fscore_rel="../../../${kind_ori}/${phe}/score_trans-ancestry.cv${cvi}/${fscore_name}"
#                #fscore="${dout}/${fscore_name}"
#
#                ln -sf $fscore_rel $dout
#                #ln -sf $fscore_rel $fscore
#            done
#        done
#    done
#done
