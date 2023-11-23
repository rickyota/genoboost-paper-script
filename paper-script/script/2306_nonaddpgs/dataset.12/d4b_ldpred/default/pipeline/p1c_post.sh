#!/bin/bash

# rename and format wgt

set -eux

program="ldpred"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh

echo "NOT TESTED"

phes="$1"
cvn="$2"
#job_id="$4"

#rhos="0.1 0.01 0.001"
rhos_e=(1.0000e-01 1.0000e-02 1.0000e-03)
rhos=(0.1 0.01 0.001)

dresult="./result/nonaddpgs/${dataname}/${program}/"

for cvi in $(seq 0 $((cvn - 1))); do
    for phe in $phes; do
        dout="${dresult}/${phe}/train.cv${cvi}/raw/"
        dout_raw="${dout}/raw/"

        cat ${dout_raw}/wgt_LDpred-inf.txt |
            tr -s ' ' '\t' | sed -e 's/chrom_//' |
            awk '{ if (NR==1){ print "chrom\tpos\tsid\tA1\tA2\tbeta\twgt" }else{ print $0 }}' \
                >"${dout}/${program}_rho-1.0.wgt"

        for ((ti = 0; ti < ${#rhos_e[@]}; ti++)); do
            rho_e=${rhos_e[ti]}
            rho=${rhos[ti]}

            f="${dout_raw}/wgt_LDpred_p${rho_e}.txt"
            fout="${dout}/${program}_rho-${rho}.wgt"

            if [ ! -e $f ]; then
                continue
            fi

            cat $f | tr -s ' ' '\t' |
                sed -e 's/chrom_//' |
                awk '{ if (NR==1){ print "chrom\tpos\tsid\tA1\tA2\tbeta\twgt" }else{ print $0 }}' \
                    >"$fout"

        done

    done
done
