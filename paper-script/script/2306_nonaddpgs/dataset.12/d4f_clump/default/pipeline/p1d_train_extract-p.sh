#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 16
#$ -q all.q
#$ -l hostname=a02|a03|z01|z02|g01|g02

set -eux

program="clump"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh

phes="$1"
cvn="$2"
clump_p_max="$3"
clump_ps_remain="$4"
clump_r2s="$5"

#ddata="./data/${dataname}/"
dresult="./result/nonaddpgs/${dataname}/${program}/"
mkdir -p $dresult

for cvi in $(seq 0 $((cvn - 1))); do
    for phe in $phes; do
        dout="${dresult}/${phe}/train.cv${cvi}/"

        for clump_p in $clump_ps_remain; do
            for clump_r2 in $clump_r2s; do
                echo "cvi, phe, clump_p,clump_r2:$cvi , $phe , $clump_p , $clump_r2"

                fclump_ori="${dout}/clump.p${clump_p_max}.r${clump_r2}.clumped"
                fclump="${dout}/clump.p${clump_p}.r${clump_r2}.clumped"

                cat $fclump_ori | awk -v clumpp="$clump_p" 'NR==1 || $7 < clumpp {print}' >$fclump
                # somehow cannot send clump_p into awk so directly access env
                #cat $fclump_ori | awk 'NR==1 || $7 < "'$clump_p'" {print}' >$fclump

            done
        done

    done
done
