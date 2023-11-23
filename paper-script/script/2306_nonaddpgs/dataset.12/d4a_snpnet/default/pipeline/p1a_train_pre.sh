#!/bin/bash

set -eux

program="snpnet"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh

cvn="$1"

ddata="./data/${dataname}/"
dtmp="${ddata}/tmp/${program}/"
mkdir -p $dtmp

fphe="${ddata}ukb_imp.phe"
fphe2="${dtmp}ukb_imp.phe"

# add #FID to fin_phe
cat $fphe | awk 'BEGIN {OFS = "\t"} {print $1,$0}' | sed -e "1 s/IID/#FID/" >$fphe2

# create .psam with #FID
cat ${ddata}/ukb_imp.psam | awk 'BEGIN {OFS = "\t"} {print $1,$0}' | sed -e "1 s/#IID/#FID/" | sed -e "1 s/#IID/IID/" >${dtmp}/ukb_imp.psam

# ng: fid=0
#source ./script/lib/system/module.sh load plink2
#plink2 --make-just-psam \
#    cols=fid,iid,sex \
#    --out ${dtmp}/ukb_imp \
#    --pfile ${ddata}/ukb_imp vzs

ln -sf ../../ukb_imp.pvar.zst ${dtmp}/ukb_imp.pvar.zst
ln -sf ../../ukb_imp.pgen ${dtmp}/ukb_imp.pgen

fcv="${ddata}ukb_imp.cv"

# snpnet cannot accept stdin as fsamples
for cvi in $(seq 0 $((cvn - 1))); do
    # for both sex
    cat $fcv | awk -v cv=${cvi} '($(3+cv) == "tr"){ print $1,$1 }' >${dtmp}/ukb_imp.cv${cvi}.samples
    # for female
    cat $fcv | awk -v cv=${cvi} '($2==0 && $(3+cv) == "tr"){ print $1,$1 }' >${dtmp}/ukb_imp.cv${cvi}.female.samples
done
