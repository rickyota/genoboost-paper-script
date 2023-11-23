#!/bin/bash

# Run default/p1a_train_pre.sh first

set -eux

#kind="monitor"
program="snpboost"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh


cvn="$1"

#flog="./log/nonaddpgs/boosting/$(basename $0).$(date +'%y%m%d-%H%M%S').log"
#mkdir -p "$(dirname $flog)"

ddata="./data/${dataname}/"
dtmp="${ddata}/tmp/${program}/"
mkdir -p $dtmp

fphe="${ddata}ukb_imp.phe"
fphe2="${dtmp}ukb_imp.phe"
fphe3="${dtmp}ukb_imp.cv.phe"
#fphe2_cov_standard="${dtmp}ukb_imp_cov_standard.phe"

#ln -sf ../../ukb_imp.pvar.zst ${dtmp}/ukb_imp.pvar.zst
#ln -sf ../../ukb_imp.pgen ${dtmp}/ukb_imp.pgen

# create .psam with #FID
#cat ${ddata}/ukb_imp.psam | awk 'BEGIN {OFS = "\t"} {print $1,$0}' | sed -e "1 s/#IID/#FID/" | sed -e "1 s/#IID/IID/" >${dtmp}/ukb_imp.psam

# add #FID to fin_phe
#cat $fphe | awk 'BEGIN {OFS = "\t"} {print $1,$0}' | sed -e "1 s/IID/#FID/" >$fphe2

## standardize cov (age, PC)
#source ./script/lib/system/module.sh load plink2
#plink2 --write-covar \
#    --out ${fphe2_cov_standard} \
#    --pfile ${dtmp}/ukb_imp vzs \
#    --covar ${fphe2} \
#    --covar-name sex,age,PC_1-PC_10 \
#    --covar-variance-standardize age PC_1 PC_2 PC_3 PC_4 PC_5 PC_6 PC_7 PC_8 PC_9 PC_10
#
#paste ${fphe2_cov_standard}.cov <(cat ${fphe2} | cut -f15-) > ${fphe2_cov_standard}

fcv="${ddata}ukb_imp.cv"

# for monitor

paste <(cat $fphe2) <(cat $fcv | cut -f3- | sed "s/tr/train/g" | sed "s/va/val/g"  ) >$fphe3


for cvi in $(seq 0 $((cvn - 1))); do
    # for both sex
    cat $fcv | awk -v cv=${cvi} '($(3+cv)=="tr" || $(3+cv)=="va"){ print $1,$1 }' >${dtmp}/ukb_imp.cv${cvi}.monitor.samples
    #cat $fcv | awk -v cv=${cvi} '(NR ==1||$(3+cv) == "tr"){ print $0 }' | cut -f1,$((3 + cvi)) | sed "s/tr/train/" | sed "s/va/val/" | awk '{print $1,$0}' | sed "1s/IID/#FID/" >${dtmp}/ukb_imp.cv${cvi}.monitor.samples
    #cat $fcv | awk -v cv=${cvi} '($(3+cv) == "tr"){ print $1,$1 }' >${dtmp}/ukb_imp.cv${cvi}.samples
    #cat $fcv | awk -v cv=${cvi} '($(3+cv) == "va"){ print $1,$1 }' >${dtmp}/ukb_imp.cv${cvi}.va.samples
    # for female
    cat $fcv | awk -v cv=${cvi} '($2==0 && ($(3+cv)=="tr" || $(3+cv)=="va")){ print $1,$1 }' >${dtmp}/ukb_imp.cv${cvi}.female.monitor.samples
    #cat $fcv | awk -v cv=${cvi} '(NR ==1||($2==0 && $(3+cv) == "tr")){ print $0 }' | cut -f1,$((3 + cvi)) | sed "s/tr/train/" | sed "s/va/val/" | awk '{print $1,$0}' | sed "1s/IID/#FID/" >${dtmp}/ukb_imp.cv${cvi}.female.monitor.samples
    #cat $fcv | awk -v cv=${cvi} '($2==0 && $(3+cv) == "tr"){ print $1,$1 }' >${dtmp}/ukb_imp.cv${cvi}.female.samples
    #cat $fcv | awk -v cv=${cvi} '($2==0 && $(3+cv) == "va"){ print $1,$1 }' >${dtmp}/ukb_imp.cv${cvi}.va.female.samples
done
