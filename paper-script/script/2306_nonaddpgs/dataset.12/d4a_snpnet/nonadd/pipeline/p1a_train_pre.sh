#!/bin/bash

set -eux

program="snpnet"
kind="nonadd"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

threads=110

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh
source ./script/lib/system/module.sh load plink2

cvn="$1"

ddata="./data/${dataname}/"
dtmp="${ddata}/tmp/${program}/${kind}/"
mkdir -p $dtmp

# create genotype file with x=0 if g=0,2; and x=1 if g=1

# use --export Av & create hetonly & --import-dosage
# https://zzz.bwh.harvard.edu/plink/dosage.shtml

plink2 --export Av \
    --out ${dtmp}/ukb_imp \
    --pfile ${ddata}/ukb_imp vzs \
    --threads ${threads}

fori="${dtmp}/ukb_imp.traw"
fnonadd="${dtmp}/ukb_imp_nonadd.traw"
# save all original var
cp $fori $fnonadd

# add hetonly
# .traw is count of **ref**, so make 2,1,0 -> 2,1,2
paste <(cat $fori | tail -n +2 | cut -f1-6 | awk 'BEGIN{OFS="\t"}{print $1,$2"-dom",$3,$4,$5,$6}') <(cat $fori | tail -n +2 | cut -f7- | sed -e "s/0/2/g") >>$fnonadd

# https://www.cog-genomics.org/plink/2.0/input#import_dosage
plink2 \
    --import-dosage $fnonadd \
    skip0=1 skip1=2 skip2=0 id-delim='_' \
    chr-col-num=1 pos-col-num=4 ref-first \
    --out ${dtmp}/ukb_imp_nonadd \
    --make-pgen vzs \
    --sort-vars \
    --psam ${ddata}/ukb_imp.psam \
    --threads 110

# plink2
ln -sf ./ukb_imp_nonadd.pvar.zst ${dtmp}/ukb_imp.pvar.zst
ln -sf ./ukb_imp_nonadd.pgen ${dtmp}/ukb_imp.pgen

fphe="${ddata}ukb_imp.phe"
fphe2="${dtmp}ukb_imp.phe"

# add #FID to fin_phe
cat $fphe | awk 'BEGIN {OFS = "\t"} {print $1,$0}' | sed -e "1 s/IID/#FID/" >$fphe2

# create .psam with #FID
cat ${dtmp}/ukb_imp_nonadd.psam | awk 'BEGIN {OFS = "\t"} {print $1,$0}' | sed -e "1 s/#IID/#FID/" | sed -e "1 s/#IID/IID/" >${dtmp}/ukb_imp.psam

fcv="${ddata}ukb_imp.cv"

# snpnet cannot accept stdin as fsamples
for cvi in $(seq 0 $((cvn - 1))); do
    # for both sex
    cat $fcv | awk -v cv=${cvi} '($(3+cv) == "tr"){ print $1,$1 }' >${dtmp}/ukb_imp.cv${cvi}.samples
    # for female
    cat $fcv | awk -v cv=${cvi} '($2==0 && $(3+cv) == "tr"){ print $1,$1 }' >${dtmp}/ukb_imp.cv${cvi}.female.samples
done

# /grid2/
if [ $HPCL_HOST = "pg" ]; then
    nodes="z01 z02 a02 a03"

    for node in $nodes; do
        echo $node
        ssh $node <<EOF
        mkdir -p /grid2/ricky/data/${dataname}/snpnet/nonadd/
EOF
    done

    #cd /nfs/data06/ricky/code/genetics

    for node in $nodes; do
        echo $node
        ddata="${node}:/grid2/ricky/data/${dataname}/snpnet/nonadd/"
        rsync -ahvz --copy-links ./${dtmp}/ukb_imp.{pvar.zst,psam,pgen} $ddata
    done

fi
