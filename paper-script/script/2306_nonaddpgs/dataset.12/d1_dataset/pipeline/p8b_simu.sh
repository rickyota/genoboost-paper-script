#!/bin/bash

# first time trial is in ./dataset.12/simu_save/
# at this time, I only used np.random.seed(seed) for the first time,
# and there is no reproducibility.
# Now, I recreate all simu data for reproducibility.

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 28
#$ -q all.q
#$ -l hostname=a02|a03
## run in centos8
##$ -l hostname=a02|a03|z01|z02|g01|g02

set -eux

dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh
source ./script/lib/system/conda_activate.sh genetics

if [ -v SGE_TASK_ID ]; then
	#job_id=$((SGE_TASK_ID - 1))
	threads=$NSLOTS
else
	#job_id=$2
	threads=16
fi

d="$1"

ddata="./data/${dataname}/"
ddata_tmp="${ddata}/simu/"
mkdir -p $ddata_tmp

flog="./log/nonaddpgs/dataset/${dataname}/$(basename $0).$(date +'%y%m%d-%H%M%S').log"
mkdir -p "$(dirname $flog)"


dwgt="${d}/true-wgt/"
dscore="${d}/score-genot/raw/"
#dscore="${d}/score-genot/raw/"
mkdir -p $dscore

# calculate genotype liability
export RUST_BACKTRACE=full
./projects_rust/target/release/genetics_res \
	score \
	--verbose \
	--threads $threads \
	--dir-score ${dscore} \
	--file-genot ${ddata}/ukb_imp \
	--genot-format "plink2-vzs" \
	--dir-wgt ${dwgt} \
	--nonadd

