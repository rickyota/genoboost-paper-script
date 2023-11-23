#!/bin/bash

set -eux

kind="ex-chr6-gwasmethods"
program="ldpred"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh


cvn="$1"
phes="$2"

# copy ss from gwasmethods
ddata="./data/${dataname}/"
dtmp="${ddata}/tmp/${program}/${kind}/"
mkdir -p $dtmp

dgwasmethods="/hpgwork2/ricky/gwas/ldpred.v1.0.11/ss_firth/"

# copy ss from gwasmethods
for phe in $phes; do
    for cvi in $(seq 0 $((cvn - 1))); do
        cp ${dgwasmethods}/dataset.${phe}.12.ukb_imp_com.obs.cv${cvi}.tr.unrelated.add.firth.ss.format.LDpredformat ${dtmp}/${phe}.cv${cvi}.ss.LDpred_ori
    done
done

# exclude chrom6
for phe in $phes; do
    for cvi in $(seq 0 $((cvn - 1))); do
        cat ${dtmp}/${phe}.cv${cvi}.ss.LDpred_ori | awk '(NR==1 || $2!=6){print}' >${dtmp}/${phe}.cv${cvi}.ss.LDpred
    done
done
