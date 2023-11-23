#!/bin/bash

# score flip is also here

set -eux

kind="default"
mode="noneur"
program="ldpred"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/conda_activate.sh genetics

#dataname="$1"
phes="$1"
cvn="$2"

ddata="./data/${dataname}/"
dresult="./result/nonaddpgs/${dataname}/${program}/${kind}/"

fphe="${ddata}/ukb_imp_noneur.phe"

rhos_all="1.0 0.1 0.01 0.001"

# TMP COMMENT
#for phe in $phes; do
#    for cvi in $(seq 0 $((cvn - 1))); do
#
#        dout="${dresult}/${phe}/score_${mode}.cv${cvi}/"
#        draw="${dout}/raw/"
#        #dtmp="${dout}/tmp/"
#        #mkdir -p $dtmp
#
#        for rho in $rhos_all; do
#            mv ${draw}/${program}_rho-${rho}.wgt.score ${draw}/${program}_rho-${rho}.score
#            #cp ${draw}/${program}_rho-${rho}.wgt.score ${dtmp}/${program}_rho-${rho}.score
#            gzip -f ${draw}/${program}_rho-${rho}.score
#        done
#
#    done
#done

# score filp
for phe in $phes; do
    for cvi in $(seq 0 $((cvn - 1))); do

        dout="${dresult}/${phe}/score_${mode}.cv${cvi}/"
        #dout="${dresult}/${phe}/score.cv${cvi}/"
        dtmp="${dout}/raw/"
        #dtmp="${dout}/tmp/"
        #dtmp_flip="${dout}/tmp_flipped/"
        #mkdir -p $dtmp_flip

        python -m projects_py.pgs.ldpred_score_flip \
            --dout ${dout} \
            --dscore ${dtmp} \
            --fpheno $fphe --phe $phe

        #dout="${dresult}/${kind}/${phe}/score.cv${cvi}/"
        #lr=${para/para_lr-/}
        #mv ${dout_tmp}/concat/${para}/score.nsnv/boosting_n.score ${dout}/boosting_lr-${lr}_n.score
        #mv ${dout_tmp}/concat/${para}//score.withcov.nsnv/boosting_n.score ${dout}/boosting_lr-${lr}_n.scorecov

        #gzip -f ${dout}/boosting_lr-${lr}_n.score
        #gzip -f ${dout}/boosting_lr-${lr}_n.scorecov

    done
done
