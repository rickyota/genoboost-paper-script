#!/bin/bash

# create snv file containg
# - one character: a, c, g, t
# - not ambiguous: ref/alt=a/c, a/g, c/a, c/t, ...
# - not duplicated: select snv with largest maf if several alleles are on the same position

set -eux

dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh
source ./script/lib/system/module.sh load plink2


dgenot="$1"

ddata="./data/${dataname}/"
# files not passed to the next pipeline should be in tmp
#ddata_tmp="${ddata}/tmp_$(basename $(basename $0) .sh)/"
ddata_tmp="${ddata}/tmp/"
mkdir -p $ddata_tmp

# TODO: concat chroms

## concat maf
fori=${ddata_tmp}ukb_imp.maf.snvs

zstdcat ${ddata_tmp}ukb_imp_chr1.maf.pvar.zst | head -n 1 >$fori
for chrom in {1..22}; do
    zstdcat ${ddata_tmp}ukb_imp_chr${chrom}.maf.pvar.zst | tail -n +2 >>$fori
done

## not ambiguous
snvs_noamb="${ddata_tmp}ukb_imp.acgt_namb.snvs"

cat $fori | awk '
{
    if ((($4 == "A") && ($5 == "C")) ||
        (($4 == "A") && ($5 == "G")) ||
        (($4 == "T") && ($5 == "C")) ||
        (($4 == "T") && ($5 == "G")) ||
        (($4 == "C") && ($5 == "A")) ||
        (($4 == "C") && ($5 == "T")) ||
        (($4 == "G") && ($5 == "A")) ||
        (($4 == "G") && ($5 == "T"))){
        print $0
    }
}' >${snvs_noamb}

# duplicated snvs
# judged by position (id is not same when rsid is not registered)
snvs_dup="${ddata_tmp}ukb_imp.acgt_namb.dup.snvs"

#cat $snvs_noamb | awk 'NR>1{
cat $snvs_noamb | awk '{
if ($2 == prev_pos){
    count++;
    if (count ==2){
        print prev_line
    }
    if (count >=2 ){
        print $0
    }
}else{
    count=1
}
prev_line=$0
prev_pos=$2


}' >${snvs_dup}

# maf
# maf of not ambigous snps with qc samples
for chrom in {1..22}; do
    #snvs_dup="${ddata_tmp}ukb_imp_chr${chrom}.acgt_namb.dup.snvs"
    fgenot_dup_chrom="${ddata_tmp}ukb_imp_chr${chrom}.acgt_namb.dup"

    cat $snvs_dup | cut -f3 |
        plink2 --freq \
            cols='chrom,pos,ref,alt,altfreq,nobs' \
            --out $fgenot_dup_chrom \
            --threads 24 \
            --pfile ${dgenot}ukb_imp_chr${chrom} \
            --keep ${ddata}/qc.sample \
            --extract /dev/stdin
done

# concat chroms
fafreq_dup="${ddata_tmp}ukb_imp.acgt_namb.dup.afreq"

# not tested
awk 'FNR>1|NR==1' ${ddata_tmp}ukb_imp_chr${chrom}.acgt_namb.dup.afreq >$fafreq_dup

#cat ${ddata_tmp}ukb_imp_chr${chrom}.acgt_namb.dup.afreq | head -n 1 >$fafreq_dup
#for chrom in {1..22}; do
#    cat ${ddata_tmp}ukb_imp_chr${chrom}.acgt_namb.dup.afreq | tail -n +2 >>$fafreq_dup
#done

# exclude list containing duplicated snvs with lower maf
snv_dup_exclude_freq="${ddata_tmp}ukb_imp.acgt_namb.dup.exclude.afreq"

cat ${fafreq_dup} | awk 'NR>1{
if ($2 == prev_pos){
    if ($6 > prev_max_maf){
        # prev is not max
        print prev_max_line
        prev_max_maf=$6
        prev_max_line=$0
    }else{
        print $0
    }
}else{
    prev_pos=$2
    prev_max_maf=$6
    prev_max_line=$0
}
}' >${snv_dup_exclude_freq}

snv_dup="${ddata_tmp}ukb_imp.acgt_namb.dup.exclude.snv"

cat $snv_dup_exclude_freq | tail -n +2 | cut -f3 >$snv_dup

source ./script/lib/system/module.sh load csvkit

# valid snvs without duplicated snvs
snv_valid_sort="${ddata_tmp}ukb_imp.valid_sort.snv"
# var id only
snvs_valid="${ddata}ukb_imp.valid.snvs"

if [ -s <(comm -13 --nocheck-order <(cut -f3 $snvs_noamb | sort) <(sort $snv_dup)) ]; then
    echo "some snvs only in dup.exclude.snv"
    exit 1
fi

# TODO: how to get ordered set difference in one line?
cat <(cut -f3 $snvs_noamb) <(cat $snv_dup) | sort | uniq -u >$snv_valid_sort
csvjoin -c3,1 -t -H -I -v $snvs_noamb $snv_valid_sort | csvformat -T | tail -n +2 >$snvs_valid

source ./script/lib/system/module.sh unload csvkit
