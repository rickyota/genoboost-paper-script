#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -j y
##$ -t 1:(batchn) # above
#$ -pe smp 1
##$ -q all.q
## TMP
#$ -q 2x.q
#$ -l hostname=a02|a03|z01|z02|g01|g02

# o smp1=30m?-40m/18phe
# smp4=11m30s/18phe
# -> almost the same
# It gets slower for later batches

set -eux

dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh

if [ -v SGE_TASK_ID ]; then
    batch_id=$((SGE_TASK_ID - 1))
    threads=$NSLOTS
else
    batch_id=$6
    threads=16
fi

#dataname="$1"
cvi="$1"
sex="$2" # male, female, both
phes="$3"
gmodels="$4"

ddata="./data/${dataname}/"
dresult="./result/nonaddpgs/${dataname}/assoc/"
dresult_tmp="${dresult}/tmp/"
mkdir -p $dresult_tmp

fbatch="${dresult}/snvs.batch"
if [ $HPCL_HOST = "pg" ]; then
    fgenot="/grid2/ricky/data/${dataname}/ukb_imp"
else
    fgenot="${ddata}/ukb_imp"
fi
fphe="${ddata}/ukb_imp.phe"
fcv="${ddata}/ukb_imp.cv"

for gmodel in $gmodels; do
    echo "gmodel $gmodel"

    bash ./script/lib/genetics/assoc.sh \
        $batch_id $gmodel $cvi $dresult_tmp \
        $fgenot $fphe $fcv $sex $fbatch "$phes" $threads

done
