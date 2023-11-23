#!/bin/bash

set -eux

dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh
source ./script/lib/system/module.sh load plink2

dimpinfo="$1"
info_thre="$2"

ddata="./data/${dataname}/"
ddata_tmp="${ddata}/tmp/"
mkdir -p $ddata_tmp

source ./script/lib/system/module.sh load csvkit

snvs_valid="${ddata}ukb_imp.valid.snvs"

for chrom in {1..22}; do
    #snv_valid="${ddata}ukb_imp_chr${chrom}.valid.snv"
    info="${dimpinfo}ukb_mfi_chr${chrom}_v3.txt"
    snv_valid_info_chrom="${ddata_tmp}ukb_imp_chr${chrom}.valid_info.snvs"

    # -t : sep=tab
    # -H : no header
    # -I : no inference
    # -c3,3: join on 3rd column on left and 3rd on right
    # warning here using one-column tsv
    # but no problem
    #csvjoin -c3,1 -t -H -I -v <(cat $snvs_valid | head -n 100) <(cut -f2,8 $info | head -n 100) | head
    csvjoin -c3,1 -t -H -I -v $snvs_valid <(cut -f2,8 $info) | csvformat -T | tail -n +2 |
        awk -v info_thre=$info_thre '{
        if ($6 > info_thre){
            print $0
        }
    }' | cut -f1-5 >$snv_valid_info_chrom

done

snv_valid_info="${ddata}ukb_imp.valid_info.snvs"
cat ${ddata_tmp}ukb_imp_chr{1..22}.valid_info.snvs >$snv_valid_info

source ./script/lib/system/module.sh unload csvkit
