#!/bin/bash

# rename and format wgt

set -eux

kind="simu"
program="ldpred"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh

#dataname="$1"
simu_type="$1"
phes="$2"
#cvn="$2"
#job_id="$4"

cvi=0

#rhos="0.1 0.01 0.001"
rhos_e=(1.0000e-01 1.0000e-02 1.0000e-03)
rhos=(0.1 0.01 0.001)

dresult="./result/nonaddpgs/${dataname}/${program}/${kind}/${simu_type}/"

format_wgt() {
    fin="$1"
    fout="$2"

    if [ -e $fin ]; then
        cat $fin |
            tr -s ' ' '\t' | sed -e 's/chrom_//' |
            awk '{ if (NR==1){ print "chrom\tpos\tsid\tA1\tA2\tbeta\twgt" }else{ print $0 }}' \
                >$fout
    fi
}

echo $phes
#for cvi in $(seq 0 $((cvn - 1))); do
for phe in $phes; do
    echo $phe

    #dout="${dresult}/${phe}/train.cv${cvi}/raw/"
    dout="${dresult}/${phe}/train.cv${cvi}/"
    dout_raw="${dout}/raw/"

    format_wgt ${dout_raw}/wgt_LDpred-inf.txt ${dout}/${program}_rho-1.0.wgt

    for ((ti = 0; ti < ${#rhos_e[@]}; ti++)); do
        rho_e=${rhos_e[ti]}
        rho=${rhos[ti]}

        format_wgt ${dout_raw}/wgt_LDpred_p${rho_e}.txt ${dout}/${program}_rho-${rho}.wgt

    done

done
#done
