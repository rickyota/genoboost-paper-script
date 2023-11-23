#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 28
#$ -q all.q
#$ -l hostname=a02|a03

set -eux

kind="ex-chr6-ts"
program="boosting"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh

if [ -v SGE_TASK_ID ]; then
    job_id=$((SGE_TASK_ID - 1))
    threads=$NSLOTS
else
    job_id=$6
    threads=16
fi

#dataname="$1"
#sex="$2" # male, female, both
phes="$1"
cvn="$2"
model="$3"
lrs="$4"
iters="$5"

phei=$((job_id / cvn))
cvi=$((job_id % cvn))

phes_ar=(${phes})
phe=${phes_ar[$phei]}

ddata="./data/${dataname}/"
dresult="./result/nonaddpgs/${dataname}/${program}/${kind}/"
mkdir -p $dresult

#dout="${dresult}/${model}/${phe}/score.cv${cvi}/raw/"
#mkdir -p $dout
#dtr="${dresult}/${model}/${phe}/train.cv${cvi}/raw/"

if [ $HPCL_HOST = "pg" ]; then
    fgenot="/grid2/ricky/data/${dataname}/ukb_imp"
else
    fgenot="${ddata}/ukb_imp"
fi
fphe="${ddata}ukb_imp_pcname.phe"
#fcv="${ddata}ukb_imp.cv"

# TODO: simu: 'ex-chr6' -> 'ex-chr6-230507'
#dtmp="${ddata}/tmp/boosting/ex-chr6-230507/"

flog="./log/nonaddpgs/${dataname}/${program}/${kind}/${model}/${phe}/score.cv${cvi}/score.$(date +'%y%m%d-%H%M%S').log"
mkdir -p "$(dirname $flog)"

#bash ./script/lib/boosting/score_230507.sh \
#    $phe "$lrs" "$iters" $sex $dout $dtr $fgenot $dtmp |&
#    tee $flog

# use newer score
dout="${dresult}/${model}/${phe}/score.cv${cvi}/raw/"
mkdir -p $dout
dtr="${dresult}/${model}/${phe}/train.cv${cvi}/raw/"

bash ./script/lib/boosting/score_no-monitor_use-snv-pos.sh \
    $dout $dtr $fgenot $fphe "$lrs" "$iters" $threads |& tee $flog

# in post
#for f in ${dout}/*.score; do
#    [ ! -f $f ] && continue
#    gzip -f $f
#done
