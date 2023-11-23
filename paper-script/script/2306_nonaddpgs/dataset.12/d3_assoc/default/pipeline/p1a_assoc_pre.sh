#!/bin/bash

set -eux

dataname="dataset.12"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh

ddata="./data/${dataname}/"
dresult="./result/nonaddpgs/${dataname}/assoc/"
mkdir -p $dresult

batchn=1000

fgenot="${ddata}/ukb_imp"
fbatch="${dresult}/snvs.batch"

echo $fbatch

zstdcat ${fgenot}.pvar.zst | tail -n +2 | awk -v batchn=$batchn 'BEGIN{OFS="\t"}{print $0,int(NR/batchn)}' >$fbatch
