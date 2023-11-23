#!/bin/bash

# If this is the second time to run, delete tmp/ first.

set -eux

kind="ex-chr6-ts"
program="boosting"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh
source ./script/lib/system/conda_activate.sh genetics

#dataname="$1"
model="$1"
phes="$2"
cvn="$3"

dresult="./result/nonaddpgs/${dataname}/${program}/${kind}/"

for phe in $phes; do
    for cvi in $(seq 0 $((cvn - 1))); do
        dout_raw="${dresult}/${model}/${phe}/score.cv${cvi}/raw"

        # concat
        for d in ${dout_raw}/para_lr-*/; do
            # para_lr-0.5
            para=$(basename $d)
            lr=${para/para_lr-/}

            dout="${dresult}/${model}/${phe}/score.cv${cvi}/"
            mv ${d}/boosting_n.score ${dout}/boosting_lr-${lr}_n.score
            mv ${d}/boosting_n.scorecov ${dout}/boosting_lr-${lr}_n.scorecov

            gzip -f ${dout}/boosting_lr-${lr}_n.score
            gzip -f ${dout}/boosting_lr-${lr}_n.scorecov

        done

    done
done
