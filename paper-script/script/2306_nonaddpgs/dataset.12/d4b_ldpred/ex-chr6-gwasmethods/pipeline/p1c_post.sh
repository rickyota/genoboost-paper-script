#!/bin/bash

# rename and format wgt

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

#rhos="0.1 0.01 0.001"
rhos_e=(1.0000e-01 1.0000e-02 1.0000e-03)
rhos=(0.1 0.01 0.001)
rhos_all="1.0 0.1 0.01 0.001"

dresult="./result/nonaddpgs/${dataname}/${program}/${kind}/"

for cvi in $(seq 0 $((cvn - 1))); do
    for phe in $phes; do
        dout="${dresult}/${phe}/train.cv${cvi}/"
        dout_raw="${dout}/raw/"

        cat ${dout_raw}/wgt_LDpred-inf.txt |
            tr -s ' ' '\t' | sed -e 's/chrom_//' |
            awk '{ if (NR==1){ print "chrom\tpos\tsid\tA1\tA2\tbeta\twgt" }else{ print $0 }}' \
                >"${dout_raw}/${program}_rho-1.0.wgt"

        for ((ti = 0; ti < ${#rhos_e[@]}; ti++)); do
            rho_e=${rhos_e[ti]}
            rho=${rhos[ti]}

            f="${dout_raw}/wgt_LDpred_p${rho_e}.txt"
            fout="${dout_raw}/${program}_rho-${rho}.wgt"

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

# sida -> rs
ddata="./data/${dataname}/"

for cvi in $(seq 0 $((cvn - 1))); do
    for phe in $phes; do

        dout="${dresult}/${phe}/train.cv${cvi}/"
        dout_raw="${dout}/raw/"

        for rho in $rhos_all; do
            fwgt="${dout_raw}/${program}_rho-${rho}.wgt"
            fout="${dout}/${program}_rho-${rho}.wgt"

            python -m projects_py.pgs.ldpred_sida_to_rs \
                --fout $fout \
                --fwgt $fwgt \
                --fconvert ${ddata}/ukb_imp.sida_rs.tsv

        done

    done
done

source ./script/lib/system/conda_deactivate.sh
