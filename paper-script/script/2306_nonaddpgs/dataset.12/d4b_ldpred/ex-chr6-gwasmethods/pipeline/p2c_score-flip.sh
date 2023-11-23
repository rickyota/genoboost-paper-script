#!/bin/bash

#PJM -L rscgrp=gg57,node=1
#PJM -L elapse=2:00:00
#PJM -j
#PJM -g gg57
#PJM --bulk

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 28
#$ -q all.q
#$ -l hostname=a02|a03

set -eux

program="ldpred"
kind="ex-chr6-gwasmethods"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh
source ./script/lib/system/conda_activate.sh genetics


phes="$1"
cvn="$2"

ddata="./data/${dataname}/"
dresult="./result/nonaddpgs/${dataname}/${program}/${kind}/"

fphe="${ddata}/ukb_imp.phe"

for phe in $phes; do
    for cvi in $(seq 0 $((cvn - 1))); do

        dout="${dresult}/${phe}/score.cv${cvi}/"
        dtmp="${dout}/tmp/"
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
